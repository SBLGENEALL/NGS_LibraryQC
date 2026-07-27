# NGS Library QC

Reference-defined amplicon library를 사내 Linux 워크스테이션에서 분석하는
오프라인 CLI 파이프라인입니다. Python 표준 라이브러리로 FASTQ 분석을
수행하고, R로 최종 figure를 만듭니다.

## 최종 프로젝트 구조

설치 후 `/data/user/MCET03/03_NGS/01_5UTR_Plasmid`에는 필요한 다섯 폴더만
남습니다.

```text
01_5UTR_Plasmid/
├── raw_data/
│   ├── 1pct/
│   └── full/
├── reference/
├── scripts/
├── config/
└── results/
    ├── 1pct/
    └── full/
```

기존의 `pipeline`, `subset_1pct`, 날짜·버전별 중복 결과와 기타 항목은
덮어쓰거나 삭제하지 않고 `/data/user/MCET03/03_NGS/archive`로 이동합니다.

## 설치

GitHub ZIP을 `/data/user/MCET03/03_NGS`에 풀고 다음 명령을 실행합니다.

```bash
cd /data/user/MCET03/03_NGS
bash NGS_LibraryQC-final-clean/install.sh
```

질문에는 `yes`를 입력합니다. 분석이 실행 중이면 설치가 중단되므로 기존
분석이 끝난 뒤 실행해야 합니다.

## 현재 1% 결과 다시 그리기

기존 분석은 다시 돌리지 않고 R figure만 새로 만듭니다.

```bash
bash /data/user/MCET03/03_NGS/01_5UTR_Plasmid/scripts/run.sh plot
```

기존 `ALL_ASSIGNED_variant_counts.tsv`와 새 CSV 결과를 모두 입력으로 읽을 수
있습니다. 새로 만들어지는 요약 표는 CSV입니다.

```text
results/1pct/figures/
├── 01_coverage_classes.png
├── 02_rank_abundance.png
├── 03_count_distribution.png
├── 04_length_gc_vs_count.png
├── coverage_class_summary.csv
├── plot_statistics.csv
└── library_qc_figures.pdf
```

Coverage plot은 모든 bar label을 plot 안에 표시할 공간을 고정적으로
확보하므로 가장 높은 `≥100 reads` 값도 잘리지 않습니다.

## 최종 파이프라인을 1% FASTQ로 검증

현재 1% 결과를 archive로 옮기고 동일한 1% FASTQ를 처음부터 다시
분석합니다.

```bash
bash /data/user/MCET03/03_NGS/01_5UTR_Plasmid/scripts/run.sh 1pct --replace
```

정상 완료 조건:

- 로그가 `[7/7]`까지 완료
- `resource_usage.txt`의 `Exit status`가 `0`
- `results/1pct/combined/ALL_ASSIGNED_variant_counts.csv` 생성
- `results/1pct/figures`의 PNG 4개와 PDF 생성

전체 FASTQ 분석은 필요할 때만 다음 명령으로 실행합니다.

```bash
bash /data/user/MCET03/03_NGS/01_5UTR_Plasmid/scripts/run.sh full
```

### 대용량 FASTQ 병렬 분석

`amplicon_qc_parallel.py`는 원본 분석 로직과 결과 형식을 유지하면서 read
pair batch를 여러 프로세스에서 병렬 처리합니다. 기본 worker 수는 서버에서
보이는 logical CPU의 절반이며 최대 128입니다. 256-thread 서버에서는
자동으로 128 worker를 사용합니다.

```bash
python3 amplicon_qc_parallel.py \
  --config config/libraryqc.ini \
  --input-dir raw_data/full \
  --outdir results/full_parallel
```

worker 수를 직접 지정할 수도 있습니다.

```bash
python3 amplicon_qc_parallel.py --workers 128 \
  --config config/libraryqc.ini \
  --input-dir raw_data/full \
  --outdir results/full_parallel
```

## 입력

### FASTQ

`fastq`, `fastq.gz`, `fq`, `fq.gz`를 지원하며 R1/R2/I1/I2, lane, chunk,
Undetermined sample을 자동 인식합니다.

### Reference

CSV, TSV, FASTA를 입력으로 지원합니다. ID와 sequence 열은 실제 header를
읽어 자동 탐지합니다. 권장 header는 다음과 같습니다.

```text
Variant_ID
5UTR candidate sequence
```

채팅에 옮겨 적은 파일명이나 header를 코드 판단에 사용하지 않습니다.
실제 reference 파일의 header만 읽습니다.

### 선택 입력

- Illumina section형 `SampleSheet.csv`
- `Top_Unknown_Barcodes.csv`
- I1/I2 FASTQ

SampleSheet와 Top Unknown Barcode는 선택 입력이며 사용하지 않아도 분석이
중단되지 않습니다.

## 결과 형식

사용자가 제공하는 reference는 TSV여도 되지만 파이프라인이 생성하는 모든
표 형식 결과는 CSV로 통일합니다. 대용량 observed sequence table만 gzip으로
압축된 CSV입니다.

```text
results/1pct/
├── run_qc_summary.csv
├── per_cycle_qc.csv
├── design_sequence_duplicates.csv
├── observed_target_sequences.csv.gz
├── index_qc/
│   ├── index_summary.csv
│   ├── index_position_substitutions.csv
│   ├── nearest_expected_barcodes.csv
│   └── top_observed_barcodes.csv
├── samples/
│   └── <sample>/
│       ├── <sample>_library_metrics.csv
│       ├── <sample>_variant_counts.csv
│       ├── <sample>_dropout_variants.csv
│       ├── <sample>_top_variants.csv
│       └── <sample>_classification.csv
├── combined/
│   ├── ALL_ASSIGNED_variant_counts.csv
│   ├── ALL_WITH_UNDETERMINED_variant_counts.csv
│   ├── all_sample_metrics.csv
│   ├── combined_metrics.csv
│   ├── variant_count_matrix.csv
│   └── variant_rpm_matrix.csv
├── figures/
├── RESULTS_TO_SHARE.txt
├── report.html
├── manifest.json
├── resource_usage.txt
└── run.log
```

## 분석 범위

- read/base/cycle Q20·Q30 및 염기 조성
- primer anchor 기반 target 복원
- paired-end quality-aware consensus
- exact 및 configurable near-match assignment
- coverage, dropout, CV, Gini, P90/P10, top read share
- Shannon effective library size
- 길이·GC와 abundance의 Spearman correlation
- SampleSheet/I1/I2/header/unknown barcode 기반 index QC
- PhiX 직접 분류 및 non-PhiX 기준 mapping·Undetermined 비율
- 실행시간, peak memory, 결과 크기, 입력/reference provenance 기록

## 테스트

```bash
bash tests/run_test.sh
```

테스트는 synthetic FASTQ end-to-end 분석, 올바른 reference header 자동 인식,
Illumina SampleSheet 처리, CSV-only output을 확인합니다. R이 설치된
환경에서는 실제 1% coverage 분포와 같은 `11 / 113 / 864 / 1,013` fixture로
label 잘림 없이 figure가 생성되는지도 검증합니다.
