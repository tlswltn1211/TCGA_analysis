"""Flet desktop UI for exploratory TCGA survival analysis."""

from __future__ import annotations

import asyncio
import os
from pathlib import Path

import certifi

# Some managed Windows installations expose an invalid CA chain to urllib.
# Flet Desktop downloads its client once on first launch, so use the bundled CA set.
os.environ.setdefault("SSL_CERT_FILE", certifi.where())

import flet as ft
import pandas as pd

from tcga_analysis import (
    SurvivalResult,
    analyze_gene,
    available_genes,
    download_timer2_clinical,
    merge_data,
    plot_pan_cancer_expression,
    plot_survival,
    search_pan_cancer_expression,
    search_timer2_gene_expression,
    TIMER2_COHORTS,
)


ACCENT = "#4F46E5"
SURFACE = "#F8FAFC"


def main(page: ft.Page) -> None:
    page.title = "TCGA Explorer_JS"
    page.theme_mode = ft.ThemeMode.LIGHT
    page.bgcolor = SURFACE
    page.padding = 0
    page.window.width = 1180
    page.window.height = 820
    page.window.min_width = 900
    page.window.min_height = 650

    picker = ft.FilePicker()
    page.services.append(picker)
    state: dict[str, object | None] = {
        "clinical": None,
        "expression": None,
        "merged": None,
        "result": None,
        "export_table": None,
        "export_name": None,
        "plot_bytes": None,
        "plot_name": None,
    }
    selected_cohorts = {"BRCA", "COAD", "KIRC", "LIHC", "LUAD", "LUSC"}

    clinical_name = ft.Text("선택되지 않음", color=ft.Colors.GREY_600)
    expression_name = ft.Text("검색되지 않음", color=ft.Colors.GREY_600)
    status = ft.Text("암종과 유전자를 선택해 검색하세요.", color=ft.Colors.GREY_700)
    gene_query = ft.TextField(
        label="유전자 기호",
        hint_text="예: TP53, EGFR, BRCA1",
        capitalization=ft.TextCapitalization.CHARACTERS,
        width=220,
    )
    cancer_select = ft.Dropdown(
        label="TCGA 암종",
        value="LUAD",
        options=[
            ft.DropdownOption(key=cohort, text=cohort)
            for cohort in TIMER2_COHORTS
            if cohort in selected_cohorts
        ],
        enable_search=True,
        width=180,
    )
    timer_button = ft.Button("TIMER2.0 임상 데이터 받기", icon=ft.Icons.CLOUD_DOWNLOAD)
    search_button = ft.Button("Tumor / Normal 발현 비교", icon=ft.Icons.COMPARE_ARROWS)
    survival_button = ft.OutlinedButton("선택 암종 생존 분석", icon=ft.Icons.QUERY_STATS)
    export_button = ft.OutlinedButton("결과 CSV 저장", icon=ft.Icons.DOWNLOAD, disabled=True)
    save_plot_button = ft.OutlinedButton("그래프 PNG 저장", icon=ft.Icons.IMAGE_OUTLINED, disabled=True)
    plot_area = ft.Container(
        content=ft.Column(
            [ft.Icon(ft.Icons.INSERT_CHART_OUTLINED, size=64, color=ft.Colors.GREY_300), ft.Text("분석 결과가 여기에 표시됩니다.", color=ft.Colors.GREY_500)],
            horizontal_alignment=ft.CrossAxisAlignment.CENTER,
            alignment=ft.MainAxisAlignment.CENTER,
        ),
        bgcolor=ft.Colors.WHITE,
        border=ft.Border.all(1, ft.Colors.GREY_200),
        border_radius=16,
        padding=20,
        expand=True,
    )
    metrics = ft.Row(spacing=12, wrap=True)

    def snack(message: str, error: bool = False) -> None:
        page.show_dialog(ft.SnackBar(ft.Text(message), bgcolor=ft.Colors.RED_700 if error else ft.Colors.GREEN_700))

    def metric_card(label: str, value: str) -> ft.Container:
        return ft.Container(
            content=ft.Column([ft.Text(label, size=12, color=ft.Colors.GREY_600), ft.Text(value, size=22, weight=ft.FontWeight.BOLD)], spacing=3),
            bgcolor=ft.Colors.WHITE,
            border=ft.Border.all(1, ft.Colors.GREY_200),
            border_radius=12,
            padding=14,
            width=150,
        )

    def prepare_data() -> None:
        clinical = state["clinical"]
        expression = state["expression"]
        if not isinstance(clinical, pd.DataFrame) or not isinstance(expression, pd.DataFrame):
            return
        merged = merge_data(clinical, expression)
        state["merged"] = merged
        genes = available_genes(merged)
        status.value = f"공통 샘플 {len(merged):,}개 · 분석 가능한 유전자 {len(genes):,}개"
        status.color = ACCENT
        metrics.controls = [metric_card("공통 샘플", f"{len(merged):,}"), metric_card("유전자", f"{len(genes):,}"), metric_card("사망 이벤트", f"{int(merged['event'].sum()):,}")]
        page.update()

    async def load_timer2_clinical(_=None) -> None:
        cohort = cancer_select.value or "LUAD"
        timer_button.disabled = True
        cancer_select.disabled = True
        status.value = f"TIMER2.0 {cohort} 임상 아카이브를 받는 중입니다…"
        status.color = ACCENT
        page.update()
        try:
            clinical = await asyncio.to_thread(download_timer2_clinical, cohort)
            state["clinical"] = clinical
            state["expression"] = None
            state["merged"] = None
            state["result"] = None
            expression_name.value = "검색되지 않음"
            export_button.disabled = True
            save_plot_button.disabled = True
            clinical_name.value = f"TIMER2.0 / GDAC {cohort} · {len(clinical):,}명"
            clinical_name.color = ft.Colors.GREY_900
            status.value = f"{cohort} 임상 데이터 {len(clinical):,}명 로드 완료"
            snack(f"TIMER2.0 {cohort} 임상 데이터를 불러왔습니다.")
        except Exception as exc:
            snack(str(exc), error=True)
            status.value = "TIMER2.0 임상 데이터를 불러오지 못했습니다."
            status.color = ft.Colors.RED_700
        finally:
            timer_button.disabled = False
            cancer_select.disabled = False
            page.update()

    timer_button.on_click = load_timer2_clinical

    def run_analysis() -> None:
        merged = state["merged"]
        if not isinstance(merged, pd.DataFrame):
            return
        try:
            genes = available_genes(merged)
            result = analyze_gene(merged, genes[0])
            state["result"] = result
            state["export_table"] = result.table
            state["export_name"] = f"{result.gene}_survival_result.csv"
            image = plot_survival(result)
            state["plot_bytes"] = image
            state["plot_name"] = f"{result.gene}_survival_plot.png"
            plot_area.content = ft.Image(src=image, fit=ft.BoxFit.CONTAIN, border_radius=10)
            metrics.controls = [
                metric_card("분석 샘플", f"{result.sample_count:,}"),
                metric_card("High / Low", f"{result.high_count} / {result.low_count}"),
                metric_card("중앙값 cutoff", f"{result.cutoff:.3f}"),
                metric_card("Log-rank p", f"{result.p_value:.3g}"),
            ]
            export_button.disabled = False
            save_plot_button.disabled = False
            status.value = f"{result.gene} 분석 완료 · χ²={result.logrank_chi2:.3f}"
            page.update()
        except Exception as exc:
            snack(str(exc), error=True)

    async def search_survival(_=None) -> None:
        cohort = cancer_select.value or "LUAD"
        gene = (gene_query.value or "").strip().upper()
        survival_button.disabled = True
        timer_button.disabled = True
        cancer_select.disabled = True
        gene_query.disabled = True
        status.value = f"{cohort} / {gene or '유전자'} 발현 데이터를 검색 중입니다…"
        status.color = ACCENT
        page.update()
        try:
            clinical = state["clinical"]
            loaded_cohort = None
            if isinstance(clinical, pd.DataFrame) and "cancer_type" in clinical:
                loaded_cohort = str(clinical["cancer_type"].iloc[0])
            if loaded_cohort != cohort:
                clinical = await asyncio.to_thread(download_timer2_clinical, cohort)
                state["clinical"] = clinical
                clinical_name.value = f"TIMER2.0 / GDAC {cohort} · {len(clinical):,}명"
                clinical_name.color = ft.Colors.GREY_900

            expression = await asyncio.to_thread(search_timer2_gene_expression, cohort, gene)
            state["expression"] = expression
            canonical_gene = next(column for column in expression.columns if column != "sample_id")
            gene_query.value = canonical_gene
            expression_name.value = f"Firebrowse RSEM {canonical_gene} · {len(expression):,}개 종양 샘플"
            expression_name.color = ft.Colors.GREY_900
            prepare_data()
            run_analysis()
            snack(f"{cohort} / {canonical_gene} 검색과 분석을 완료했습니다.")
        except Exception as exc:
            snack(str(exc), error=True)
            status.value = str(exc)
            status.color = ft.Colors.RED_700
        finally:
            survival_button.disabled = False
            timer_button.disabled = False
            cancer_select.disabled = False
            gene_query.disabled = False
            page.update()

    async def export_result(_=None) -> None:
        table = state["export_table"]
        if not isinstance(table, pd.DataFrame):
            return
        file_name = str(state["export_name"] or "tcga_result.csv")
        csv_bytes = table.to_csv(index=False).encode("utf-8-sig")
        path = await picker.save_file(
            dialog_title="분석 결과 저장",
            file_name=file_name,
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["csv"],
            src_bytes=csv_bytes,
        )
        if path:
            Path(path).write_bytes(csv_bytes)
            snack(f"저장 완료: {path}")

    async def save_plot(_=None) -> None:
        image = state["plot_bytes"]
        if not isinstance(image, bytes):
            return
        file_name = str(state["plot_name"] or "tcga_plot.png")
        path = await picker.save_file(
            dialog_title="분석 그래프 저장",
            file_name=file_name,
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["png"],
            src_bytes=image,
        )
        if path:
            Path(path).write_bytes(image)
            snack(f"그래프 저장 완료: {path}")

    async def search_pan_cancer(_=None) -> None:
        gene = (gene_query.value or "").strip().upper()
        cohorts = [cohort for cohort in TIMER2_COHORTS if cohort in selected_cohorts]
        search_button.disabled = True
        survival_button.disabled = True
        gene_query.disabled = True
        status.value = f"{len(cohorts)}개 암종에서 {gene or '유전자'} Tumor/Normal 발현을 검색 중입니다…"
        status.color = ACCENT
        page.update()
        try:
            expression = await asyncio.to_thread(search_pan_cancer_expression, cohorts, gene)
            canonical_gene = str(expression["gene"].iloc[0]).upper()
            image, stats = await asyncio.to_thread(plot_pan_cancer_expression, expression, canonical_gene)
            gene_query.value = canonical_gene
            state["export_table"] = expression
            state["export_name"] = f"{canonical_gene}_pan_cancer_expression.csv"
            state["result"] = None
            state["plot_bytes"] = image
            state["plot_name"] = f"{canonical_gene}_pan_cancer_expression.png"
            plot_area.content = ft.Image(src=image, fit=ft.BoxFit.CONTAIN, border_radius=10)
            tumor_n = int(expression["sample_type"].eq("TP").sum())
            normal_n = int(expression["sample_type"].eq("NT").sum())
            significant = int(stats["p_value"].lt(0.05).sum())
            metrics.controls = [
                metric_card("선택 암종", f"{len(cohorts)}"),
                metric_card("Tumor", f"{tumor_n:,}"),
                metric_card("Normal", f"{normal_n:,}"),
                metric_card("p < 0.05", f"{significant}"),
            ]
            expression_name.value = f"Firebrowse RSEM {canonical_gene} · 총 {len(expression):,}개 샘플"
            expression_name.color = ft.Colors.GREY_900
            export_button.disabled = False
            save_plot_button.disabled = False
            status.value = (
                f"{canonical_gene} pan-cancer 발현 비교 완료 · Mann–Whitney "
                "* p<0.05, ** p<0.005, *** p<0.0005"
            )
            snack(f"{canonical_gene}의 {len(cohorts)}개 암종 발현 비교를 완료했습니다.")
        except Exception as exc:
            snack(str(exc), error=True)
            status.value = str(exc)
            status.color = ft.Colors.RED_700
        finally:
            search_button.disabled = False
            survival_button.disabled = False
            gene_query.disabled = False
            page.update()

    search_button.on_click = search_pan_cancer
    gene_query.on_submit = search_pan_cancer
    survival_button.on_click = search_survival
    export_button.on_click = export_result
    save_plot_button.on_click = save_plot

    cohort_summary = ft.Text(
        ", ".join(cohort for cohort in TIMER2_COHORTS if cohort in selected_cohorts),
        size=12,
        color=ft.Colors.GREY_700,
    )
    cohort_checks = {
        cohort: ft.Checkbox(label=cohort, value=cohort in selected_cohorts, col=3)
        for cohort in TIMER2_COHORTS
    }

    def select_all_cohorts(_=None) -> None:
        for checkbox in cohort_checks.values():
            checkbox.value = True
        page.update()

    def clear_cohorts(_=None) -> None:
        for checkbox in cohort_checks.values():
            checkbox.value = False
        page.update()

    def apply_cohorts(_=None) -> None:
        chosen = {cohort for cohort, checkbox in cohort_checks.items() if checkbox.value}
        if not chosen:
            snack("하나 이상의 암종을 선택하세요.", error=True)
            return
        selected_cohorts.clear()
        selected_cohorts.update(chosen)
        ordered = [cohort for cohort in TIMER2_COHORTS if cohort in selected_cohorts]
        cohort_summary.value = f"{len(ordered)}개 선택: " + ", ".join(ordered)
        cancer_select.options = [ft.DropdownOption(key=cohort, text=cohort) for cohort in ordered]
        if cancer_select.value not in selected_cohorts:
            cancer_select.value = ordered[0]
        page.pop_dialog()
        page.update()

    cohort_dialog = ft.AlertDialog(
        modal=True,
        title="분석할 TCGA 암종 선택",
        content=ft.Container(
            content=ft.ResponsiveRow(list(cohort_checks.values()), spacing=2, run_spacing=2),
            width=680,
            height=390,
        ),
        actions=[
            ft.TextButton("전체 선택", on_click=select_all_cohorts),
            ft.TextButton("모두 해제", on_click=clear_cohorts),
            ft.Button("적용", on_click=apply_cohorts),
        ],
        scrollable=True,
    )

    def open_cohort_dialog(_=None) -> None:
        page.show_dialog(cohort_dialog)

    timer2_card = ft.Container(
        content=ft.Column([
            ft.Row([ft.Icon(ft.Icons.DNS_OUTLINED, color=ACCENT), ft.Text("TIMER2.0 임상 데이터", size=17, weight=ft.FontWeight.BOLD)]),
            ft.Text("TIMER2.0이 사용한 TCGA GDAC Firehose 2016-01-28 임상 자료", size=12, color=ft.Colors.GREY_600),
            ft.Row([cancer_select, timer_button], wrap=True),
            survival_button,
            clinical_name,
            ft.Text("OS = 사망 시 days_to_death, 생존 시 days_to_last_followup", size=11, color=ft.Colors.GREY_500),
        ], spacing=10),
        bgcolor=ft.Colors.WHITE,
        border=ft.Border.all(1, ft.Colors.GREY_200),
        border_radius=14,
        padding=18,
    )

    gene_search_card = ft.Container(
        content=ft.Column([
            ft.Row([ft.Icon(ft.Icons.MANAGE_SEARCH, color=ACCENT), ft.Text("유전자 발현 검색", size=17, weight=ft.FontWeight.BOLD)]),
            ft.Text("선택한 암종의 Tumor와 Normal RSEM 발현값을 동시에 비교", size=12, color=ft.Colors.GREY_600),
            gene_query,
            ft.OutlinedButton("암종 선택", icon=ft.Icons.CHECKLIST, on_click=open_cohort_dialog),
            cohort_summary,
            search_button,
            expression_name,
        ], spacing=10),
        bgcolor=ft.Colors.WHITE,
        border=ft.Border.all(1, ft.Colors.GREY_200),
        border_radius=14,
        padding=18,
    )

    header = ft.Container(
        content=ft.Row([
            ft.Column([ft.Text("TCGA Explorer_JS", size=26, weight=ft.FontWeight.BOLD, color=ft.Colors.WHITE), ft.Text("유전자 발현 기반 생존 분석", color="#C7D2FE")], spacing=2),
            ft.Container(expand=True),
            ft.TextButton("TIMER2.0 논문", icon=ft.Icons.OPEN_IN_NEW, style=ft.ButtonStyle(color=ft.Colors.WHITE), url="https://doi.org/10.1093/nar/gkaa407"),
        ]),
        bgcolor=ACCENT,
        padding=ft.Padding.symmetric(horizontal=28, vertical=18),
    )

    sidebar = ft.Container(
        content=ft.Column([
            ft.Text("1. 데이터", size=18, weight=ft.FontWeight.BOLD),
            timer2_card,
            ft.Divider(),
            ft.Text("2. 유전자 검색", size=18, weight=ft.FontWeight.BOLD),
            gene_search_card,
            ft.Row([export_button, save_plot_button], wrap=True),
            ft.Container(content=status, bgcolor="#EEF2FF", border_radius=10, padding=12),
        ], spacing=14, scroll=ft.ScrollMode.AUTO),
        width=360,
        padding=22,
        bgcolor="#F1F5F9",
    )

    content = ft.Column([
        ft.Text("TCGA 발현·생존 분석", size=22, weight=ft.FontWeight.BOLD),
        ft.Text("선택 암종의 Tumor/Normal 발현 비교 또는 개별 암종 생존 분석을 수행합니다.", color=ft.Colors.GREY_600),
        metrics,
        plot_area,
        ft.Text("탐색적 연구용 도구이며 임상적 의사결정에 사용할 수 없습니다.", size=11, color=ft.Colors.GREY_500),
    ], spacing=14, expand=True)

    page.add(header, ft.Row([sidebar, ft.Container(content=content, padding=24, expand=True)], spacing=0, expand=True))


if __name__ == "__main__":
    ft.run(main)
