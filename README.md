# Amplicon Library QC Pipeline

외부 인터넷이나 별도 Python 패키지 없이 사내 Linux 서버에서 실행하는
설정 파일 기반 amplicon-library 분석 파이프라인입니다. 5′UTR에 한정되지
않고 promoter, enhancer, barcode, regulatory element, 짧은 CDS fragment 등
reference-defined amplicon library에 재사용할 수 있습니다.

**Current version: v1.3 (2026-07-26)**

v1.3은 Illumina section형 SampleSheet의 surplus-column 오류를 수정하고,
오프라인 워크스테이션의 복잡한 폴더를 자동 탐색·정리하며 1%/full 분석을
한 명령으로 실행하도록 만든 버전입니다.

버전은 두 자리만 사용합니다. 세부 수정은 `v1.3 → v1.4`, 큰 구조 변경은
`v2.0 → v3.0`으로 진행합니다.

## GitHub ZIP으로 설치하고 기존 폴더 정리

GitHub에서 `v1.3` 브랜치의 ZIP을 내려받아 사내 워크스테이션의
`/data/user/MCET03/03_NGS` 바로 아래에 풉니다. 압축을 풀면 일반적으로
`NGS_LibraryQC-1.3` 폴더가 생깁니다.

먼저 읽기 전용 자동 탐색을 실행합니다.

```bash
cd /data/user/MCET03/03_NGS
bash NGS_LibraryQC-1.3/cleanup_project.sh
```

출력에는 다음 후보가 자동으로 표시됩니다.

- 원본 FASTQ 그룹
- 추출된 1% FASTQ 그룹
- reference 후보
- SampleSheet 및 index metadata 후보
- 핵심 combined 파일이 존재하는 기존 분석 결과
- 기존 figure/report

파일명을 입력할 필요가 없습니다. 실제 정리를 시작하려면:

```bash
bash NGS_LibraryQC-1.3/cleanup_project.sh apply
```

각 후보를 보존할지 `yes` 또는 `no`로만 선택합니다. 최종 구조는 다음처럼
번호 중복 없이 정리됩니다.

```text
/data/user/MCET03/03_NGS/
└── 01_5UTR_Plasmid/
    ├── raw_data/
    ├── reference/
    ├── subset_1pct/
    ├── results/
    ├── config/
    └── pipeline/
```

`00_Tools`, `99_ProjectTemplate`, `create_ngs_project.sh`, 중복 pipeline과
선택하지 않은 파일은 즉시 삭제하지 않고
`/data/user/MCET03/03_NGS_cleanup_quarantine_<timestamp>`로 이동합니다.
새 1% 분석이 `[7/7]`까지 완료된 뒤에만 다음 명령으로 각 quarantine을
`yes/no` 확인 후 영구 삭제합니다.

```bash
bash /data/user/MCET03/03_NGS/01_5UTR_Plasmid/pipeline/purge_quarantine.sh
```

## 정리 후 분석 실행

1% 분석:

```bash
bash /data/user/MCET03/03_NGS/01_5UTR_Plasmid/pipeline/run_project.sh 1pct
```

최초 실행에서는 자동으로 reference와 SampleSheet 후보를 보여주며,
역시 `yes/no`로만 선택합니다. 실행이 끝나면 다음을 자동으로 수행합니다.

- `[7/7]` 완료 검증
- Python 표준 라이브러리 기반 실제 시간과 peak memory 기록
- 결과 폴더 크기 기록 및 full FASTQ 대비 전체 실행시간·공간 추정
- R이 설치되어 있으면 coverage class, rank-abundance, count distribution,
  length/GC plot 생성

1% 결과를 확인한 뒤 전체 분석:

```bash
bash /data/user/MCET03/03_NGS/01_5UTR_Plasmid/pipeline/run_project.sh full
```

분석 날짜는 실행한 날짜를 사용하며 결과는
`results/YYYYMMDD_1pct_v1.3` 또는 `results/YYYYMMDD_full_v1.3`에 저장됩니다.
같은 날 다시 실행하면 `_rerun2`, `_rerun3`이 자동으로 붙어 기존 결과를
덮어쓰지 않습니다.

## 저장소 폴더 구조

```text
.
├── amplicon_qc.py
├── run_pipeline.sh
├── README.md
├── VERSION
├── CHANGELOG.md
├── config/
├── docs/
├── references/
└── tests/
```

분석에서 생성되는 FASTQ, 결과 폴더, HTML/TSV 보고서는 저장소에 포함하지
않습니다.

## 분석 범위

- R1/R2/I1/I2 FASTQ 자동 탐색 및 sample/lane/chunk 묶기
- FASTQ별 read 수, 길이, mean Q, Q20/Q30, N 비율
- cycle별 Q20/Q30, A/C/G/T/N 조성
- primer anchor 기반 target 서열 복원
- paired-end overlap의 quality-aware consensus
- reference에 exact match 및 configurable near match
- 검출률, dropout, depth, CV, Gini, P90/P10, top 1%/10% 점유율
- Shannon effective library size
- 길이·GC 함량과 read count의 Spearman correlation
- SampleSheet, I1/I2 FASTQ, FASTQ header, Top Unknown Barcodes를 이용한
  index 오류 진단
- `Undetermined`에서 reference-matching read를 별도 계산
- bundled PhiX174 `NC_001422.1` reference에 R1/R2를 직접 분류
- 전체 및 sample별 PhiX 검출률, non-PhiX:PhiX ratio, PhiX 제외 mapping률
- PhiX를 제외한 실제 `Undetermined` non-PhiX 비율 계산
- HTML 보고서, TSV 결과, 실행 manifest, 복사용 요약 생성

Python 3.9 이상과 표준 라이브러리만 사용합니다.

## 1. 입력

### FASTQ

일반적인 Illumina 파일명을 자동 인식합니다.

```text
fastq/
├── SampleA_S1_L001_R1_001.fastq.gz
├── SampleA_S1_L001_R2_001.fastq.gz
├── SampleB_S2_L001_R1_001.fastq.gz
├── SampleB_S2_L001_R2_001.fastq.gz
├── Undetermined_S0_L001_R1_001.fastq.gz
└── Undetermined_S0_L001_R2_001.fastq.gz
```

Single-end R1도 처리할 수 있지만, target 전체가 한 read 안에 있지 않으면
reference assignment가 감소합니다.

### Reference

CSV, TSV 또는 FASTA를 지원합니다.

```csv
Variant_ID,target_sequence
VAR_0001,ACGTACGTACGT
VAR_0002,TGCATGCATGCA
```

- ID는 고유해야 합니다.
- 같은 서열에 ID가 여러 개면 sequencing으로 구분할 수 없으므로 하나의
  sequence group으로 집계하고 경고합니다.
- `reference_mode = target`은 target 서열만 있는 reference입니다.
- `reference_mode = amplicon`은 primer가 포함된 full amplicon reference입니다.
- `reference_mode = auto`는 두 형식을 자동 판정합니다.
- FASTA를 사용할 때는 record name이 variant ID가 됩니다.

### 선택 입력

- `SampleSheet.csv`: expected i7/i5 index
- `Top_Unknown_Barcodes.csv`: demultiplexing 결과의 unknown index와 count
- I1/I2 FASTQ: raw index cycle 품질 및 염기 조성

## 2. 설정

템플릿을 복사합니다.

```bash
cp config/config_template.ini config.ini
```

`config.ini`에서 프로젝트마다 주로 다음 값만 바꿉니다.

```ini
[project]
target_name = promoter

[input]
input_dir = /path/to/fastq
reference = /path/to/reference.csv
outdir = /path/to/results
sample_sheet = /path/to/SampleSheet.csv
top_unknown = /path/to/Top_Unknown_Barcodes.csv

[reference]
id_column = Variant_ID
sequence_column = promoter_sequence
reference_mode = target

[amplicon]
left_primer = FORWARD_PRIMER_SEQUENCE
right_primer = REVERSE_PRIMER_SEQUENCE

[control_qc]
phix_enabled = yes
expected_phix_percent = 10
```

`right_primer`에는 reverse-complement된 amplicon 서열이 아니라 실제 주문한
reverse primer oligo 서열을 입력합니다. Stagger primer의 5′ random N은
제외하고 template에 결합하는 고정 부분만 입력합니다.

현재 5′UTR 실험 조건은 `config/config_5utr_example.ini`에 예제로
포함되어 있습니다.

## 3. 실행

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh config.ini
```

또는 Python 명령으로 직접 실행할 수 있습니다.

```bash
python3 amplicon_qc.py --config config.ini
```

설정 일부만 일시적으로 덮어쓸 수도 있습니다.

```bash
python3 amplicon_qc.py \
  --config config.ini \
  --outdir /path/to/test_results \
  --max-reads 100000 \
  --overwrite
```

첫 실행은 `max_reads = 100000`으로 파일 인식과 primer-anchor rate를 확인한
뒤, 새 output 경로에서 `max_reads`를 비우고 전체 분석하는 것을 권장합니다.

설정 파일 없이 모든 옵션을 CLI로 입력하는 방식도 지원합니다.

```bash
python3 amplicon_qc.py \
  --input-dir /path/to/fastq \
  --reference /path/to/reference.csv \
  --outdir /path/to/results \
  --left-primer ACGT... \
  --right-primer TGCA... \
  --id-column Variant_ID \
  --sequence-column target_sequence
```

## 4. 결과 구조

```text
results/
├── RESULTS_TO_SHARE.txt
├── report.html
├── manifest.json
├── run_qc_summary.tsv
├── per_cycle_qc.tsv
├── design_sequence_duplicates.tsv
├── observed_target_sequences.tsv.gz
├── index_qc/
│   ├── index_summary.tsv
│   ├── index_position_substitutions.tsv
│   ├── nearest_unknown_barcodes.tsv
│   └── top_observed_barcodes.tsv
├── samples/
│   └── <sample>/
│       ├── <sample>_library_metrics.tsv
│       ├── <sample>_variant_counts.tsv
│       ├── <sample>_dropout_variants.tsv
│       ├── <sample>_top_variants.tsv
│       └── <sample>_classification.tsv
└── combined/
    ├── ALL_ASSIGNED_variant_counts.tsv
    ├── ALL_WITH_UNDETERMINED_variant_counts.tsv
    ├── all_sample_metrics.tsv
    ├── variant_count_matrix.tsv
    ├── variant_rpm_matrix.tsv
    └── combined_metrics.tsv
```

`ALL_ASSIGNED`는 정상 demultiplexing된 sample을 합친 결과입니다.
`ALL_WITH_UNDETERMINED`는 unknown index로 빠진 read까지 같은 reference에
매칭한 진단용 결과입니다. 한 run에 같은 reference를 공유하는 여러 library가
있다면 `ALL_WITH_UNDETERMINED`를 공식 sample count로 사용하면 안 됩니다.
`variant_count_matrix.tsv`와 `variant_rpm_matrix.tsv`는 sample을 열로 배치한
downstream 통계·시각화용 matrix입니다.

## 5. 주요 지표

| 지표 | 의미 |
|---|---|
| `q30_base_percent` | base-call 정확도 Q30 이상 비율 |
| `mapped_fraction_percent` | 전체 read pair 중 reference에 배정된 비율 |
| `phix_fraction_percent` | 전체 read pair 중 PhiX174로 직접 분류된 비율 |
| `mapped_fraction_non_phix_percent` | PhiX를 분모에서 제외한 target mapping률 |
| `undetermined_non_phix_fraction_percent` | PhiX 제거 후 남은 실제 index-failure 후보 비율 |
| `detected_1x_percent` | 1 read 이상 검출된 설계 서열 비율 |
| `dropout_sequences` | 한 번도 검출되지 않은 unique reference 수 |
| `coverage_cv` | coverage 표준편차/평균 |
| `gini` | 0에 가까울수록 균등, 1에 가까울수록 편중 |
| `p90_over_p10` | 낮을수록 균등; P10이 0이면 계산 불가 |
| `top_10_percent_read_share_percent` | 상위 10% variant의 read 점유율 |
| `effective_library_size` | Shannon entropy 기반 실질 다양성 |
| `spearman_length_logcount` | 길이와 abundance의 순위상관 |
| `spearman_gc_logcount` | GC와 abundance의 순위상관 |
| `both_anchors_percent` | 양쪽 primer evidence가 모두 확인된 read 비율 |

Near match는 sequencing error rescue 용도입니다. mutation discovery나
variant calling을 목적으로 할 때는 reference assignment 결과만으로 새
mutation을 해석하지 말고 `observed_target_sequences.tsv.gz`를 별도로
검토해야 합니다.

PhiX는 bundled `references/NC_001422.1_phiX174.fasta`의 circular genome
27-mer를 이용해 R1/R2에서 직접 분류합니다. `expected_phix_percent`를
입력하면 `RESULTS_TO_SHARE.txt`에서 관측값과 투입값의 percentage-point
차이를 함께 표시합니다. `max_reads`를 sample마다 적용한 시험 실행에서는
sample과 Undetermined를 같은 수만큼 subsampling할 수 있으므로 전체 run의
PhiX 비율을 해석하지 말고, 제한을 제거한 full run 결과만 사용해야 합니다.

## 6. Index 문제 해석

가능한 정보원은 우선순위가 다릅니다.

1. I1/I2 FASTQ: raw index read의 cycle별 품질과 base composition
2. FASTQ header index: demultiplexed read에 기록된 observed index
3. Top Unknown Barcodes: demultiplexing 실패 read의 상위 index

`expected_G_to_N`과 `expected_G_to_other`는 expected index의 G 위치에서
나타난 오류를 집계합니다. Top Unknown Barcodes만 사용했을 경우 이는 전체
index error rate가 아니라 demultiplexing에 실패한 read 내부의 패턴입니다.

## 7. 다른 library에 재사용할 때

코드는 수정하지 않고 새 config와 reference만 만듭니다.

- 5′UTR: `target_name`, 경로, sequence column만 변경
- promoter library: promoter primer와 promoter reference 입력
- barcode library: barcode 양쪽 고정 primer와 barcode reference 입력
- enhancer library: enhancer-flanking primer와 reference 입력

Target 길이는 PE read로 양 끝에서 충분히 덮여야 정확한 full-sequence
assignment가 가능합니다. 현재 matcher는 reference 길이의 95% 이상에
evidence가 있을 때만 배정합니다.

## 8. 검증

포함된 synthetic paired-end dataset으로 전체 동작을 검사합니다.

```bash
bash tests/run_test.sh
```

성공하면 마지막에 다음 메시지가 출력됩니다.

```text
Synthetic integration test passed
```

## 9. 제한사항

- raw BCL demultiplexing 자체는 수행하지 않습니다.
- adapter/quality trimming된 새 FASTQ 파일을 생성하지 않습니다. Primer anchor
  사이의 target은 분석 중에 직접 추출하므로 reference quantification에는
  별도 primer trimming이 필수는 아닙니다.
- paired FASTQ 파일의 read 순서가 동기화되어 있다고 가정합니다.
- 긴 amplicon이 R1/R2 합산 범위를 벗어나면 full assignment가 감소합니다.
- indel이 많은 library에는 reference signature 방식보다 aligner 기반
  workflow가 적합합니다.
- exact 230–300 bp physical size selection 같은 wet-lab QC는 다루지 않습니다.
