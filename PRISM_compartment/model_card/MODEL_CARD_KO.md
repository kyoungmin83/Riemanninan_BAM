# PRISM integrated compartmental learned-rank — SV6 epoch 52

## 한 줄 정의

`compartmental e52`는 canonical Phase-I 병리 decoder와 Phase-II PRISM 분해를 하나의 연속 학습으로 결합하고, 논문에서 가져온 compartmental-threshold 계산 경로, train-only technical-zero bank, 자동 pathology-rank 학습, 자동 generator-cardinality 학습을 함께 사용한 SV6 checkpoint다.

이 모델은 이전의 `fixed-rank8 e52`, SV7 `graph-linear control`, frozen `personal-rank2 e20`과 서로 다른 checkpoint다. 이름의 `e52`는 epoch 번호이고 rank 수를 뜻하지 않는다.

## 정확한 위치

- 모델 checkpoint: `/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_freegen_recovery_e12_20260829/checkpoint_epoch_052.pt`
- 전체 run: `/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_freegen_recovery_e12_20260829`
- 학습 당시 source tree: `/home/kmlee/project_sv6/kmlee_bam_integrated_compartmental_pathrank_20260828`
- exact source cfg 사본: `cfg/source_config.json`
- runtime-resolved cfg 사본: `cfg/resolved_config.json`
- checkpoint manifest: `CHECKPOINT_MANIFEST.json`

## 무슨 모델인가

### 공통 골격

- 13,498 genes, 24 cell types, 414 candidate biological modules를 사용한다.
- ordinal gene-expression decoder와 cell-type reference anchoring을 사용한다.
- sex·technology nuisance control과 reference-centered clean latent를 사용한다.
- PRISM module activity를 shared pathology, personal baseline, candidate personal pathology response로 분리한다.
- personal response rank는 2로 고정돼 있다.

### 병리 rank

- canonical Phase I에서는 pathology rank 8을 사용했다.
- epoch 25부터 96개 후보 방향의 soft rank search를 시작했다.
- epoch 33부터 hard mask를 사용했고 epoch 37에 hard rank 8로 확정·동결했다.
- e52에서 네 pathology-correction axis 모두 95% energy rank가 6이고, forward hard rank는 8이다.

즉 이 모델은 처음부터 rank 8을 고정한 이전 모델과 결과 rank 수는 같지만, rank 후보를 학습한 뒤 8개를 선택했다는 점이 다르다.

### Generator 개수

- 후보 generator는 414개다.
- epoch 37부터 generator-cardinality 학습을 시작했다.
- 목표 개수와 최소 개수는 지정하지 않았다.
- blanket unique-gene-coverage 보호는 제거했고 singleton generator 10개만 보호했다.
- e52의 epoch-end 확정값은 hard 112개, expected 109.49개다.

### 논문 아이디어 적용부

- 참고 논문: `ref/Dendritic morphology and synaptic nonlinearities enhancefunctional complexity in human cortical neurons.pdf`
- model arm: `compartmental_threshold`
- module graph를 compartment처럼 사용하고 학습 가능한 threshold, slope, gain, inter-compartment mix를 적용한다.
- 이것은 dendrite나 NMDA를 직접 측정한 생물물리 모델이 아니라 논문 가설을 옮긴 계산적 비유다.
- 이 경로의 독립 기여는 아직 graph-linear arm과 동일 checkpoint stage에서 최종 비교되지 않았다. 따라서 전체 성능 향상을 비선형 경로 하나의 인과 효과라고 부르면 안 된다.

### Technical-zero

- 관측된 zero가 technical dropout인지 추정하는 train-only latent bank를 사용한다.
- cell type과 region을 조건화하고, 후보가 부족하면 zero로 fail-closed한다.
- sex-linked genes에는 별도 차단 규칙을 적용한다.
- official test donor를 이웃 bank에 넣지 않는다.

## 학습 일정

| Epoch | 주요 단계 |
|---:|---|
| 1–12 | Canonical Phase I biological warmup; pathology decoder rank 8 |
| 13–16 | PRISM common/personal 분해 cross-fade |
| 17–24 | Module-local 경로 개방; nonlinear path는 21–24 ramp |
| 25–32 | Pathology-rank soft search와 sparsity |
| 33–36 | Hard rank mask |
| 37 | Pathology rank 8 확정·동결; generator search 시작 |
| 37–40 | Generator shadow stage |
| 41–44 | Generator soft pruning |
| 45–55 | Generator hard pruning |

e52는 pathology rank가 이미 확정되고 generator hard pruning이 진행된 후기 checkpoint다.

## 왜 e52를 보존했는가

e52는 학습 중 자동 `checkpoint_best`가 아니다. 자동 primary loss 기준의 best epoch는 13이었다. e52는 test를 보기 전에 validation module recovery를 기준으로 `biological-recovery candidate` 역할로 별도 보존했다.

| Validation metric | e52 |
|---|---:|
| AD-module recovery | 0.8593 |
| Pooled module recovery | 0.8640 |
| Cell-centered recovery | 0.8359 |
| Disease-isolated mean | 0.7633 |
| Validation reconstruction NLL | 0.5576 |

이미 소비된 동일 test set에 대한 사후 재비교 결과는 다음과 같다.

| Consumed-test metric | rank2 e20 | e52 |
|---|---:|---:|
| AD-module recovery | 0.5266 | 0.7257 |
| Pooled module recovery | 0.5178 | 0.7016 |
| Cell-centered recovery | 0.5767 | 0.7851 |
| Disease-isolated mean | 0.5235 | 0.7393 |
| Reconstruction NLL | 0.5852 | 0.5572 |

e52는 test에서 rank2보다 18/24 cell types의 AD-module recovery가 높았다. 비신경세포 6종의 평균/최솟값은 0.716/0.528이었다.

## 반드시 함께 읽어야 할 한계

1. 공식 test 9명은 2026-08-28에 이미 개봉됐다. 현재 test 결과는 `CONSUMED_REUSE`이며 새로운 독립 holdout이 아니다.
2. test 결과로 e52를 다시 선택하거나 baseline 교체를 정당화할 수 없다. e52의 역할은 test 전에 validation으로 고정됐다.
3. clean latent의 test UMI Spearman은 0.927, detected-gene Spearman은 0.915다. sequencing-depth coupling이 매우 강해 생물학 복원 개선의 일부가 depth 구조를 이용했을 가능성을 배제하지 못한다.
4. sex leakage는 chance에 가깝지만 technology leakage는 e52에서 0.649로 높다.
5. candidate personal pathology-response의 독립 재현성은 확립되지 않았다.
6. compartmental-threshold 경로의 독립 기여는 같은 stage의 graph-linear 대조군과 parameter/function-class를 명확히 구분해 평가해야 한다.

## 현재 판정

- **용도:** 후기 biological-recovery 후보, compartmental-threshold 실험 checkpoint, consumed-test descriptive comparison.
- **강점:** 지금까지 가장 높은 validation 및 consumed-test module recovery, rank 8과 약 112 generators로 압축, 낮은 reconstruction NLL.
- **약점:** 심한 depth coupling, 소비된 test, personal-response 및 nonlinear-path 인과 기여 미확정.
- **공식 baseline:** 새로운 untouched cohort 또는 사전 고정 외부 holdout 전까지 frozen rank2 e20을 독립 baseline으로 유지한다.

## 관련 문서와 결과

- 설계: `docs/prism_integrated_compartmental_nonlinearity_techzero_v2_design_20260827.md`
- 구현 보고서: `docs/prism_compartmental_techzero_v2_implementation_report_20260827.md`
- 통합 curriculum: `docs/prism_integrated_celltype_ad_module_curriculum_v2_design_20260824.md`
- e12 recovery: `analysis_outputs/sv6_compartmental_recovery_20260829/RECOVERY_REPORT_KO.md`
- e52 validation: `analysis_outputs/prism_integrated_compartmental_pathrank_sv6_e52_interim_posthoc_20260901`
- e52 consumed-test: `analysis_outputs/prism_integrated_compartmental_pathrank_consumed_test_20260902/e52`
- rank2 및 이전 모델 비교: `analysis_outputs/prism_integrated_compartmental_pathrank_consumed_test_20260902/CONSUMED_TEST_COMPARISON_REPORT_KO.md`

## 무결성

- checkpoint SHA-256: `9fbb7df5119aba87a52d437fe443d656044e3c9c868fe13b1c39b73ac375944d`
- resolved cfg SHA-256: `0618ced5ac4fd7a096422ddddf3ac02e631d305c0e05596b6b04ad8de3e66846`
- source cfg SHA-256: `de65ca10d906eeb168294af25214dcf377a779296a8c6640a8603c668f5e494a`

모델 checkpoint는 627,238,305 bytes이며 model-card 폴더에는 중복 복사하지 않았다. cfg와 manifest만 보존하고 checkpoint는 위의 shared-storage 원본을 사용한다.
