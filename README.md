# TCGA Explorer_JS

Flet으로 만든 로컬 데스크톱 TCGA 생존 분석 앱입니다. TIMER2.0이 사용한 공식 TCGA GDAC Firehose 임상 자료와 Firebrowse RSEM 발현값을 검색·결합하고, 선택한 유전자의 중앙값을 기준으로 발현군을 나눈 뒤 Kaplan–Meier 곡선과 log-rank 검정 결과를 제공합니다.

## Windows 실행 파일

Python이나 Git 없이 사용하려면 [GitHub Releases](https://github.com/tlswltn1211/TCGA_analysis/releases/latest)에서 `TCGA-Explorer_JS.exe`를 내려받아 더블클릭하세요. Windows SmartScreen이 표시되면 파일 출처를 확인한 후 **추가 정보 → 실행**을 선택할 수 있습니다.

## 실행

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe app.py
```

현재 PC에서는 다음 명령으로 바로 실행할 수 있습니다.

```powershell
C:\Users\user\TCGA_analysis\.venv\Scripts\python.exe C:\Users\user\TCGA_analysis\app.py
```

## 입력 데이터

임상 데이터는 앱에서 TCGA 암종을 선택한 후 **TIMER2.0 임상 데이터 받기**를 누르면 가져옵니다. TIMER2.0 논문에 명시된 GDAC Firehose 2016-01-28 `Clinical_Pick_Tier1` 아카이브를 사용합니다.

- 환자 ID: `bcr_patient_barcode`
- 생존 사건: `vital_status`
- 사망 환자 생존일: `days_to_death`
- 생존 환자 추적일: `days_to_last_followup`

## 사용 방법

1. TCGA 암종을 선택합니다.
2. **암종 선택**에서 비교할 암종만 체크합니다.
3. 유전자 기호(예: `TP53`, `EGFR`)를 입력합니다.
4. **Tumor / Normal 발현 비교**를 누릅니다.

앱은 선택한 모든 암종의 primary tumor(`TP`)와 adjacent normal(`NT`) RSEM 발현값을 동시에 표시하고 Mann–Whitney 검정을 수행합니다. 개별 암종은 **선택 암종 생존 분석**으로 Kaplan–Meier 분석할 수 있습니다. TCGA 샘플 바코드는 환자 단위 12자리 바코드로 자동 정규화됩니다.

분석 결과는 CSV로, 현재 그래프는 PNG 파일로 저장할 수 있습니다.

데이터 출처: [TIMER2.0 논문](https://doi.org/10.1093/nar/gkaa407), [Broad GDAC Firehose](https://gdac.broadinstitute.org/)

## 검증

```powershell
.\.venv\Scripts\python.exe -m pytest
```

이 앱은 탐색적 연구용이며 의료 진단이나 임상적 의사결정용이 아닙니다.
