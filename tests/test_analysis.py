from io import BytesIO
import tarfile

import numpy as np

from tcga_analysis import (
    analyze_gene,
    available_genes,
    make_example,
    merge_data,
    normalize_tcga_id,
    parse_timer2_clinical_archive,
    plot_pan_cancer_expression,
    plot_survival,
    search_pan_cancer_expression,
    search_timer2_gene_expression,
)


def test_example_analysis_end_to_end():
    clinical, expression = make_example(size=48)
    merged = merge_data(clinical, expression)
    assert len(merged) == 48
    assert available_genes(merged) == ["TP53", "BRCA1", "EGFR"]

    result = analyze_gene(merged, "TP53")
    assert result.sample_count == 48
    assert result.high_count + result.low_count == 48
    assert np.isfinite(result.p_value)
    assert 0 <= result.p_value <= 1
    assert plot_survival(result).startswith(b"\x89PNG")


def test_no_shared_samples_raises():
    clinical, expression = make_example(size=8)
    expression["sample_id"] = "OTHER-" + expression["sample_id"]
    try:
        merge_data(clinical, expression)
    except ValueError as exc:
        assert "공통" in str(exc)
    else:
        raise AssertionError("Expected validation error")


def test_timer2_firehose_archive_conversion():
    matrix = (
        "Hybridization REF\ttcga-aa-0001\ttcga-aa-0002\n"
        "Composite Element REF\tvalue\tvalue\n"
        "vital_status\t1\t0\n"
        "days_to_death\t120\tNA\n"
        "days_to_last_followup\tNA\t450\n"
        "years_to_birth\t61\t55\n"
    ).encode()
    archive_bytes = BytesIO()
    with tarfile.open(fileobj=archive_bytes, mode="w:gz") as archive:
        info = tarfile.TarInfo("TEST.clin.merged.picked.txt")
        info.size = len(matrix)
        archive.addfile(info, BytesIO(matrix))

    clinical = parse_timer2_clinical_archive(archive_bytes.getvalue(), "LUAD")
    assert clinical["sample_id"].tolist() == ["TCGA-AA-0001", "TCGA-AA-0002"]
    assert clinical["overall_survival_days"].tolist() == [120, 450]
    assert clinical["cancer_type"].eq("LUAD").all()
    assert normalize_tcga_id("tcga-aa-0001-01a") == "TCGA-AA-0001"


def test_gene_search_response_conversion(monkeypatch):
    class Response:
        def raise_for_status(self):
            return None

        def json(self):
            return {
                "mRNASeq": [
                    {"tcga_participant_barcode": f"TCGA-AA-{index:04d}", "expression_log2": index / 10, "gene": "TP53"}
                    for index in range(1, 7)
                ]
            }

    monkeypatch.setattr("tcga_analysis.httpx.get", lambda *args, **kwargs: Response())
    expression = search_timer2_gene_expression("LUAD", "tp53")
    assert expression.columns.tolist() == ["sample_id", "TP53"]
    assert len(expression) == 6


def test_pan_cancer_tumor_normal_plot(monkeypatch):
    class Response:
        def __init__(self, cohort):
            self.cohort = cohort

        def raise_for_status(self):
            return None

        def json(self):
            rows = []
            for index in range(1, 9):
                rows.append({
                    "tcga_participant_barcode": f"TCGA-{self.cohort[:2]}-{index:04d}",
                    "expression_log2": index / 2 + (1 if self.cohort == "LUAD" else 0),
                    "gene": "TP53",
                    "cohort": self.cohort,
                    "sample_type": "TP" if index <= 5 else "NT",
                })
            return {"mRNASeq": rows}

    def fake_get(*args, **kwargs):
        return Response(kwargs["params"]["cohort"])

    monkeypatch.setattr("tcga_analysis.httpx.get", fake_get)
    expression = search_pan_cancer_expression(["LUAD", "BRCA"], "TP53")
    image, stats = plot_pan_cancer_expression(expression, "TP53")
    assert image.startswith(b"\x89PNG")
    assert stats["cohort"].tolist() == ["LUAD", "BRCA"]
    assert stats["tumor_n"].tolist() == [5, 5]
    assert stats["normal_n"].tolist() == [3, 3]
