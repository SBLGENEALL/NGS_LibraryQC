# NextSeq amplicon NGS QC와 trimming 가이드

## 1. 분석 단계

```text
Run folder/InterOp
  → BCL conversion + demultiplexing
  → raw FASTQ QC
  → PhiX 분리
  → 필요 시 adapter/primer/quality trimming
  → post-trim QC
  → amplicon target 추출 및 reference 정량
  → library representation/bias 평가
```

현재 파이프라인은 demultiplexing이 끝난 FASTQ부터 시작한다. BCL conversion
자체와 cleaned FASTQ 생성은 수행하지 않지만, primer anchor 사이의 target을
메모리에서 직접 복원하므로 designed amplicon library 정량에는 별도 primer
trimming이 필수는 아니다.

## 2. PhiX와 9:1 확인

전체 FASTQ의 read pair를 다음처럼 상호배타적으로 분류한다.

- PhiX
- designed reference exact match
- designed reference near match
- unresolved non-PhiX

계산식:

```text
Observed PhiX % = PhiX read pairs / all read pairs × 100
Observed library:PhiX = non-PhiX read pairs / PhiX read pairs
Adjusted Undetermined % =
  Undetermined non-PhiX read pairs / all non-PhiX read pairs × 100
```

90:10으로 투입했다면 기대 ratio는 9:1이다. 다만 wet-lab molarity 비율과
cluster/read 비율은 농도 측정 오차, fragment 특성, cluster generation 차이
때문에 완전히 같지 않을 수 있다. SAV의 `% Aligned`와 full-FASTQ direct
PhiX fraction을 함께 비교한다.

Unindexed PhiX는 보통 Undetermined FASTQ에 들어가므로 raw Undetermined
비율을 곧바로 index failure로 해석하면 안 된다.

## 3. 단계별 QC 항목

### A. Run-level: SAV 또는 InterOp

- total reads와 reads passing filter
- yield와 projected yield
- cluster density 또는 occupancy
- `% PF`
- R1/R2와 I1/I2의 `%Q30`
- cycle별 Q-score
- cycle별 A/C/G/T/N 조성
- phasing/prephasing
- PhiX aligned %
- PhiX error rate
- cycle 1 intensity와 low-diversity 이상

### B. Demultiplexing/index

- assigned와 Undetermined read 수 및 비율
- perfect index match
- one-mismatch와 two-mismatch index reads
- Top Unknown Barcodes
- expected index까지의 Hamming distance
- index position별 substitution
- expected G 위치의 `G→N`, `G→other`
- i5 orientation
- PhiX를 제거한 adjusted Undetermined 비율

### C. Raw FASTQ

- R1/R2 pair 수와 pair synchronization
- read length 분포
- mean Q, Q20, Q30
- cycle별 Q30
- N-base와 reads-with-N
- per-cycle base composition과 GC
- adapter content
- overrepresented sequence
- duplication

### D. Trimming 후

- retained read pairs와 retained base %
- adapter-trimmed %
- forward/reverse primer detected %
- too-short 또는 discarded read %
- length distribution
- R1/R2 pair synchronization
- post-trim Q30/N/adapter content

### E. Designed library

- target extraction과 both-primer-anchor %
- exact/near mapping
- PhiX 제외 mapping %
- unresolved non-PhiX %
- detected 1×, 10×, 100×, 500×
- dropout
- mean/median/min/max coverage
- CV, Gini, P90/P10, P95/P5
- top 1%/10% read share
- effective library size
- length 및 GC bias
- original-vector count와 contamination %
- assigned와 Undetermined variant-count correlation

## 4. Trimming 여부 결정

Trimming은 항상 수행하는 단계가 아니다.

### Adapter trimming이 필요한 경우

- FastQC/MultiQC에서 adapter content가 증가한다.
- insert가 read 길이보다 짧아 adapter read-through가 예상된다.
- downstream aligner 또는 merger가 adapter 때문에 실패한다.

### Primer trimming이 필요한 경우

- 새로운 mutation 또는 ASV를 찾는다.
- primer mismatch를 biological variant로 오인할 수 있다.
- primer를 제외한 insert FASTQ가 다른 도구의 입력으로 필요하다.

### Quality trimming이 필요한 경우

- read 말단에서 cycle별 품질이 뚜렷하게 하락한다.
- low-quality tail이 overlap consensus 또는 mapping을 방해한다.

Q30 미만 base를 무조건 모두 자르는 방식은 권장하지 않는다. 80–130 bp
amplicon의 2×150 데이터는 R1/R2 overlap이 충분하므로 quality-aware consensus가
오류를 교정할 수 있다. 지나친 trimming은 primer anchor 또는 짧은 variant를
잃게 만들 수 있다.

## 5. 널리 사용하는 도구

- FastQC: 개별 raw/post-trim FASTQ QC
- MultiQC: 여러 FastQC, Cutadapt, fastp 결과를 한 HTML로 통합
- Cutadapt: 알려진 adapter와 PCR primer를 정확하게 제거
- fastp: adapter/quality trimming과 기본 QC를 한 번에 수행
- nf-core/ampliseq: FastQC–Cutadapt–MultiQC–DADA2/QIIME2를 묶은
  재현 가능한 microbiome amplicon workflow

`nf-core/ampliseq`는 미지의 ASV와 taxonomy 분석에 적합하다. 현재처럼
주문한 2,001개 서열에 exact/near assignment하고 representation을 평가하는
실험에는 이 패키지의 reference-defined counting 단계가 더 직접적이다.

## 6. 권장 운영 방식

1. SAV/InterOp와 BCL Convert report를 보존한다.
2. raw FASTQ를 read-only 원본으로 보존한다.
3. raw FastQC/MultiQC를 만든다.
4. full FASTQ에서 PhiX를 직접 정량한다.
5. adapter/primer 문제가 확인될 때만 cleaned FASTQ를 별도 폴더에 만든다.
6. post-trim FastQC/MultiQC로 개선과 read loss를 동시에 확인한다.
7. reference-defined library QC를 실행한다.
8. 모든 config, 명령, 버전, report를 run별 결과 폴더에 보존한다.

`max_reads` 시험 실행은 primer와 reference 설정 확인용이다. Sample별로 같은
수의 read를 자르면 실제 assigned/Undetermined/PhiX 비율이 왜곡되므로 최종
run-level 비율은 반드시 제한 없는 full run에서 계산한다.
