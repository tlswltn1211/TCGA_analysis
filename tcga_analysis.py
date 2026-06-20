"""Core data validation and survival-analysis functions for the TCGA app."""

from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path
import re
import tarfile

import httpx
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, mannwhitneyu


CLINICAL_COLUMNS = {"sample_id", "overall_survival_days", "vital_status"}
TIMER2_COHORTS = (
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
    "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
    "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
    "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM",
)
TIMER2_FIREHOSE_ROOT = (
    "https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data"
)
FIREBROWSE_MRNASEQ_URL = "http://firebrowse.org/api/v1/Samples/mRNASeq"


@dataclass(frozen=True)
class SurvivalResult:
    gene: str
    cutoff: float
    sample_count: int
    event_count: int
    high_count: int
    low_count: int
    logrank_chi2: float
    p_value: float
    table: pd.DataFrame


def read_csv(path: str | Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="cp949")


def normalize_tcga_id(value: object) -> str:
    """Normalize sample IDs to the 12-character TCGA participant barcode."""
    sample_id = str(value).strip().upper().replace("_", "-")
    tcga_barcode = re.match(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", sample_id)
    return tcga_barcode.group(1) if tcga_barcode else sample_id


def timer2_clinical_url(cohort: str) -> str:
    cohort = cohort.upper()
    if cohort not in TIMER2_COHORTS:
        raise ValueError(f"지원하지 않는 TCGA 암종입니다: {cohort}")
    archive = (
        f"gdac.broadinstitute.org_{cohort}.Clinical_Pick_Tier1."
        "Level_4.2016012800.0.0.tar.gz"
    )
    return f"{TIMER2_FIREHOSE_ROOT}/{cohort}/20160128/{archive}"


def parse_timer2_clinical_archive(content: bytes, cohort: str) -> pd.DataFrame:
    """Convert the GDAC Firehose matrix used by TIMER2.0 to one row per patient."""
    try:
        with tarfile.open(fileobj=BytesIO(content), mode="r:gz") as archive:
            members = [
                member
                for member in archive.getmembers()
                if member.isfile() and member.name.endswith(".clin.merged.picked.txt")
            ]
            if len(members) != 1:
                raise ValueError("TIMER2 임상 행렬을 아카이브에서 찾을 수 없습니다.")
            extracted = archive.extractfile(members[0])
            if extracted is None:
                raise ValueError("TIMER2 임상 행렬을 읽을 수 없습니다.")
            matrix = pd.read_csv(extracted, sep="\t", index_col=0, dtype=str)
    except (tarfile.TarError, OSError) as exc:
        raise ValueError("유효한 TIMER2/GDAC 임상 아카이브가 아닙니다.") from exc

    clinical = matrix.T
    clinical = clinical.drop(columns=["Composite Element REF"], errors="ignore")
    clinical.index.name = "sample_id"
    clinical = clinical.reset_index()
    required = {"sample_id", "vital_status", "days_to_death", "days_to_last_followup"}
    missing = required.difference(clinical.columns)
    if missing:
        raise ValueError(f"TIMER2 임상 필수 변수가 없습니다: {', '.join(sorted(missing))}")

    clinical["sample_id"] = clinical["sample_id"].map(normalize_tcga_id)
    clinical["vital_status"] = pd.to_numeric(clinical["vital_status"], errors="coerce")
    death = pd.to_numeric(clinical["days_to_death"], errors="coerce")
    followup = pd.to_numeric(clinical["days_to_last_followup"], errors="coerce")
    clinical["overall_survival_days"] = death.where(
        clinical["vital_status"].eq(1), followup
    ).fillna(death).fillna(followup)
    clinical["cancer_type"] = cohort.upper()
    if "years_to_birth" in clinical:
        clinical = clinical.rename(columns={"years_to_birth": "age_at_diagnosis"})
    return clinical


def download_timer2_clinical(cohort: str) -> pd.DataFrame:
    """Download the official GDAC Firehose clinical archive used by TIMER2.0."""
    url = timer2_clinical_url(cohort)
    try:
        # This archived Broad endpoint has an incomplete legacy certificate chain on
        # some managed Windows machines. The host and path are fixed above.
        response = httpx.get(url, follow_redirects=True, timeout=60, verify=False)
        response.raise_for_status()
    except httpx.HTTPError as exc:
        raise ValueError(f"TIMER2 임상 데이터 다운로드 실패: {cohort}") from exc
    return parse_timer2_clinical_archive(response.content, cohort)


def search_timer2_gene_expression(cohort: str, gene: str) -> pd.DataFrame:
    """Fetch primary-tumor RSEM expression used by the TIMER2 data pipeline."""
    cohort = cohort.upper().strip()
    gene = gene.upper().strip()
    if cohort not in TIMER2_COHORTS:
        raise ValueError(f"지원하지 않는 TCGA 암종입니다: {cohort}")
    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{0,39}", gene):
        raise ValueError("유전자 기호를 확인하세요. 예: TP53, EGFR, BRCA1")
    try:
        response = httpx.get(
            FIREBROWSE_MRNASEQ_URL,
            params={
                "format": "json",
                "gene": gene,
                "cohort": cohort,
                "sample_type": "TP",
                "page_size": 2000,
            },
            timeout=60,
        )
        response.raise_for_status()
        records = response.json().get("mRNASeq", [])
    except (httpx.HTTPError, ValueError) as exc:
        raise ValueError(f"유전자 발현 데이터 검색 실패: {cohort} / {gene}") from exc
    if not records:
        raise ValueError(f"{cohort}에서 유전자 {gene}의 발현 데이터를 찾지 못했습니다.")

    expression = pd.DataFrame.from_records(records)
    required = {"tcga_participant_barcode", "expression_log2", "gene"}
    if not required.issubset(expression.columns):
        raise ValueError("Firebrowse 발현 응답 형식이 올바르지 않습니다.")
    canonical_gene = str(expression["gene"].dropna().iloc[0]).upper()
    expression = expression.rename(
        columns={"tcga_participant_barcode": "sample_id", "expression_log2": canonical_gene}
    )
    expression["sample_id"] = expression["sample_id"].map(normalize_tcga_id)
    expression[canonical_gene] = pd.to_numeric(expression[canonical_gene], errors="coerce")
    expression = (
        expression[["sample_id", canonical_gene]]
        .dropna()
        .groupby("sample_id", as_index=False)[canonical_gene]
        .mean()
    )
    if len(expression) < 4:
        raise ValueError(f"{cohort} / {canonical_gene}의 분석 가능 샘플이 4개 미만입니다.")
    return expression


def search_pan_cancer_expression(cohorts: list[str], gene: str) -> pd.DataFrame:
    """Fetch tumor and adjacent-normal RSEM expression for selected TCGA cohorts."""
    selected = list(dict.fromkeys(cohort.upper().strip() for cohort in cohorts))
    invalid = [cohort for cohort in selected if cohort not in TIMER2_COHORTS]
    if not selected:
        raise ValueError("하나 이상의 암종을 선택하세요.")
    if invalid:
        raise ValueError(f"지원하지 않는 TCGA 암종입니다: {', '.join(invalid)}")
    gene = gene.upper().strip()
    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{0,39}", gene):
        raise ValueError("유전자 기호를 확인하세요. 예: TP53, EGFR, BRCA1")

    def fetch(cohort: str) -> list[dict]:
        try:
            response = httpx.get(
                FIREBROWSE_MRNASEQ_URL,
                params={
                    "format": "json",
                    "gene": gene,
                    "cohort": cohort,
                    "sample_type": "TP,NT",
                    "page_size": 2000,
                },
                timeout=60,
            )
            response.raise_for_status()
            return response.json().get("mRNASeq", [])
        except (httpx.HTTPError, ValueError) as exc:
            raise ValueError(f"{cohort} 발현 데이터 검색 실패") from exc

    with ThreadPoolExecutor(max_workers=min(8, len(selected))) as executor:
        batches = list(executor.map(fetch, selected))
    records = [record for batch in batches for record in batch]
    if not records:
        raise ValueError(f"선택한 암종에서 유전자 {gene}의 발현 데이터를 찾지 못했습니다.")

    frame = pd.DataFrame.from_records(records)
    required = {"tcga_participant_barcode", "expression_log2", "gene", "cohort", "sample_type"}
    if not required.issubset(frame.columns):
        raise ValueError("Firebrowse 발현 응답 형식이 올바르지 않습니다.")
    frame = frame.rename(
        columns={"tcga_participant_barcode": "sample_id", "expression_log2": "expression"}
    )
    frame["sample_id"] = frame["sample_id"].map(normalize_tcga_id)
    frame["expression"] = pd.to_numeric(frame["expression"], errors="coerce")
    frame = frame[frame["sample_type"].isin(["TP", "NT"])].dropna(subset=["expression"])
    frame = frame[["sample_id", "cohort", "sample_type", "expression", "gene"]]
    frame = frame.drop_duplicates(["sample_id", "cohort", "sample_type"])
    if not frame["sample_type"].eq("TP").any():
        raise ValueError("선택한 암종의 종양 발현 데이터가 없습니다.")
    return frame


def plot_pan_cancer_expression(
    expression: pd.DataFrame, gene: str
) -> tuple[bytes, pd.DataFrame]:
    """Draw selected-cohort tumor/normal boxes and return cohort statistics."""
    cohorts = list(dict.fromkeys(expression["cohort"].astype(str)))
    figure_width = max(7.2, min(18.0, 0.82 * len(cohorts)))
    fig, ax = plt.subplots(figsize=(figure_width, 6.2), dpi=140)
    rng = np.random.default_rng(42)
    colors = {"TP": "#FF1744", "NT": "#2979FF"}
    stats: list[dict] = []
    all_box_tops: list[float] = []
    star_labels: list[tuple[int, str]] = []

    for index, cohort in enumerate(cohorts, start=1):
        cohort_rows = expression[expression["cohort"] == cohort]
        tumor = cohort_rows.loc[cohort_rows["sample_type"] == "TP", "expression"].to_numpy(float)
        normal = cohort_rows.loc[cohort_rows["sample_type"] == "NT", "expression"].to_numpy(float)
        p_value = float(mannwhitneyu(tumor, normal).pvalue) if len(tumor) and len(normal) else np.nan
        stats.append({
            "cohort": cohort,
            "tumor_n": len(tumor),
            "normal_n": len(normal),
            "tumor_median": float(np.median(tumor)) if len(tumor) else np.nan,
            "normal_median": float(np.median(normal)) if len(normal) else np.nan,
            "p_value": p_value,
        })
        for values, sample_type, position in ((tumor, "TP", index - 0.18), (normal, "NT", index + 0.18)):
            if not len(values):
                continue
            box = ax.boxplot(
                values,
                positions=[position],
                widths=0.28,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "white", "linewidth": 1.4},
            )
            box["boxes"][0].set(facecolor=colors[sample_type], alpha=0.98)
            for element in box["whiskers"] + box["caps"]:
                element.set(color=colors[sample_type])
            box_top = float(np.max(box["caps"][1].get_ydata()))
            all_box_tops.append(box_top)
            sample = values if len(values) <= 180 else rng.choice(values, 180, replace=False)
            jitter = rng.normal(position, 0.045, len(sample))
            ax.scatter(jitter, sample, s=7, alpha=0.38, color=colors[sample_type], linewidths=0)
        if np.isfinite(p_value):
            stars = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            star_labels.append((index, stars))

    star_y: float | None = None
    if star_labels:
        plot_values = expression["expression"].to_numpy(float)
        star_offset = min(0.04, max(0.01, float(np.ptp(plot_values)) * 0.005))
        star_y = max(all_box_tops) + star_offset
        for index, stars in star_labels:
            ax.text(index, star_y, stars, ha="center", va="bottom", fontsize=20, fontweight="bold")

    from matplotlib.patches import Patch

    ax.set_title(f"{gene.upper()} expression across selected TCGA cancers", fontweight="bold")
    ax.set_ylabel("RSEM expression (log2)", fontsize=15)
    ax.set_xticks(range(1, len(cohorts) + 1), cohorts, rotation=45, ha="right", fontsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.set_xlim(0.55, len(cohorts) + 0.45)
    ax.grid(axis="y", alpha=0.18)
    if star_y is not None:
        y_bottom, y_top = ax.get_ylim()
        text_padding = (y_top - y_bottom) * 0.12
        ax.set_ylim(y_bottom, max(y_top, star_y + text_padding))
    ax.legend(
        handles=[Patch(facecolor=colors["TP"], label="Tumor"), Patch(facecolor=colors["NT"], label="Normal")],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
    )
    fig.tight_layout()
    output = BytesIO()
    fig.savefig(output, format="png", bbox_inches="tight")
    plt.close(fig)
    return output.getvalue(), pd.DataFrame(stats)


def validate_clinical(df: pd.DataFrame) -> pd.DataFrame:
    missing = CLINICAL_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(f"임상 CSV 필수 열이 없습니다: {', '.join(sorted(missing))}")
    out = df.copy()
    out["sample_id"] = out["sample_id"].map(normalize_tcga_id)
    out["overall_survival_days"] = pd.to_numeric(
        out["overall_survival_days"], errors="coerce"
    )
    status = out["vital_status"].astype(str).str.strip().str.lower()
    status_map = {
        "dead": 1,
        "deceased": 1,
        "1": 1,
        "true": 1,
        "alive": 0,
        "living": 0,
        "0": 0,
        "false": 0,
    }
    out["event"] = status.map(status_map)
    out = out.dropna(subset=["sample_id", "overall_survival_days", "event"])
    out = out[out["overall_survival_days"] >= 0].drop_duplicates("sample_id")
    out["event"] = out["event"].astype(int)
    if out.empty:
        raise ValueError("분석 가능한 임상 행이 없습니다.")
    return out


def validate_expression(df: pd.DataFrame) -> pd.DataFrame:
    if "sample_id" not in df.columns:
        raise ValueError("발현 CSV에는 sample_id 열이 필요합니다.")
    if len(df.columns) < 2:
        raise ValueError("발현 CSV에는 하나 이상의 유전자 열이 필요합니다.")
    out = df.copy()
    out["sample_id"] = out["sample_id"].map(normalize_tcga_id)
    gene_columns = [column for column in out.columns if column != "sample_id"]
    for column in gene_columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    out = out.drop_duplicates("sample_id")
    valid_genes = [column for column in gene_columns if out[column].notna().sum() >= 4]
    if not valid_genes:
        raise ValueError("숫자 발현값이 4개 이상인 유전자 열이 없습니다.")
    return out[["sample_id", *valid_genes]]


def merge_data(clinical: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    clinical_valid = validate_clinical(clinical)[
        ["sample_id", "overall_survival_days", "vital_status", "event"]
    ]
    merged = clinical_valid.merge(validate_expression(expression), on="sample_id", how="inner")
    if len(merged) < 4:
        raise ValueError("두 CSV에 공통으로 존재하는 샘플이 4개 미만입니다.")
    return merged


def available_genes(merged: pd.DataFrame) -> list[str]:
    excluded = CLINICAL_COLUMNS | {"event"}
    return [column for column in merged.columns if column not in excluded]


def _km_curve(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    event_times = np.sort(np.unique(time[event == 1]))
    x = [0.0]
    y = [1.0]
    survival = 1.0
    for current in event_times:
        at_risk = int(np.sum(time >= current))
        deaths = int(np.sum((time == current) & (event == 1)))
        if at_risk:
            x.extend([float(current), float(current)])
            y.extend([survival, survival * (1 - deaths / at_risk)])
            survival = y[-1]
    if len(time):
        x.append(float(np.max(time)))
        y.append(survival)
    return np.asarray(x), np.asarray(y)


def _logrank(time: np.ndarray, event: np.ndarray, high: np.ndarray) -> tuple[float, float]:
    observed = expected = variance = 0.0
    for current in np.sort(np.unique(time[event == 1])):
        risk_high = np.sum((time >= current) & high)
        risk_low = np.sum((time >= current) & ~high)
        deaths_high = np.sum((time == current) & (event == 1) & high)
        deaths_total = np.sum((time == current) & (event == 1))
        risk_total = risk_high + risk_low
        if risk_total <= 1:
            continue
        observed += deaths_high
        expected += deaths_total * risk_high / risk_total
        variance += (
            risk_high
            * risk_low
            * deaths_total
            * (risk_total - deaths_total)
            / (risk_total**2 * (risk_total - 1))
        )
    statistic = float((observed - expected) ** 2 / variance) if variance > 0 else 0.0
    return statistic, float(chi2.sf(statistic, 1))


def analyze_gene(merged: pd.DataFrame, gene: str) -> SurvivalResult:
    if gene not in available_genes(merged):
        raise ValueError(f"유전자 열을 찾을 수 없습니다: {gene}")
    table = merged[["sample_id", "overall_survival_days", "event", gene]].dropna()
    if len(table) < 4:
        raise ValueError("해당 유전자의 분석 가능 샘플이 4개 미만입니다.")
    cutoff = float(table[gene].median())
    table = table.assign(group=np.where(table[gene] >= cutoff, "High", "Low"))
    high = table["group"].eq("High").to_numpy()
    if high.all() or (~high).all():
        raise ValueError("중앙값 기준으로 두 발현군을 만들 수 없습니다.")
    statistic, p_value = _logrank(
        table["overall_survival_days"].to_numpy(float),
        table["event"].to_numpy(int),
        high,
    )
    return SurvivalResult(
        gene=gene,
        cutoff=cutoff,
        sample_count=len(table),
        event_count=int(table["event"].sum()),
        high_count=int(high.sum()),
        low_count=int((~high).sum()),
        logrank_chi2=statistic,
        p_value=p_value,
        table=table,
    )


def plot_survival(result: SurvivalResult) -> bytes:
    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=140)
    colors = {"High": "#ef476f", "Low": "#277da1"}
    for group in ("High", "Low"):
        rows = result.table[result.table["group"] == group]
        x, y = _km_curve(
            rows["overall_survival_days"].to_numpy(float),
            rows["event"].to_numpy(int),
        )
        ax.step(x, y, where="post", linewidth=2.4, color=colors[group], label=f"{group} (n={len(rows)})")
        censored = rows[rows["event"] == 0]
        if not censored.empty:
            censor_y = [y[np.searchsorted(x, t, side="right") - 1] for t in censored["overall_survival_days"]]
            ax.scatter(censored["overall_survival_days"], censor_y, marker="|", s=75, color=colors[group])
    ax.set(title=f"{result.gene} expression and overall survival", xlabel="Overall survival (days)", ylabel="Survival probability")
    ax.set_ylim(-0.03, 1.04)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    ax.text(0.99, 0.04, f"Log-rank p = {result.p_value:.3g}", transform=ax.transAxes, ha="right")
    fig.tight_layout()
    output = BytesIO()
    fig.savefig(output, format="png", bbox_inches="tight")
    plt.close(fig)
    return output.getvalue()


def make_example(seed: int = 1211, size: int = 96) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    sample_ids = [f"TCGA-DEMO-{index:04d}" for index in range(1, size + 1)]
    tp53 = rng.normal(8.0, 1.5, size)
    brca1 = rng.normal(6.5, 1.1, size)
    egfr = rng.normal(7.2, 1.7, size)
    hazard = np.exp((tp53 - np.median(tp53)) * 0.28)
    event_time = rng.exponential(1200 / hazard)
    censor_time = rng.exponential(1800, size)
    observed = np.minimum(event_time, censor_time)
    event = event_time <= censor_time
    clinical = pd.DataFrame({
        "sample_id": sample_ids,
        "overall_survival_days": observed.round(0).astype(int),
        "vital_status": np.where(event, "Dead", "Alive"),
        "age_at_diagnosis": rng.integers(35, 86, size),
    })
    expression = pd.DataFrame({"sample_id": sample_ids, "TP53": tp53.round(4), "BRCA1": brca1.round(4), "EGFR": egfr.round(4)})
    return clinical, expression
