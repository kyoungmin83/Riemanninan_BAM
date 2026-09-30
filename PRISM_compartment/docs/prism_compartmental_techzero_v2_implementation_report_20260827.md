# PRISM compartmental nonlinearity + technical-zero bank v2 구현 보고서

작성일: 2026-08-27  
수정일: 2026-08-28 — integrated learned pathology-rank 구현 추가  
상태: **코드 및 component regression 완료, full-data 학습 미실행, launch 잠금 유지**

## 1. 무엇을 구현했는가

이 작업은 `module-local residual v1`을 논문 실험이라고 이름만 바꾼 것이 아니다.

- v1은 paper feature 이전의 선형 reference로 그대로 보존한다.
- `compartmental_threshold`는 고정 module graph, 학습 가능한 bounded mix·threshold·slope·gain, threshold 부근의 supralinear curvature를 추가한 paper-inspired arm이다.
- `graph_linear_control`은 같은 graph·입력·curriculum을 사용하면서 `mix`와 `gain`만 학습하는 함수족 대조군이다. threshold·slope tensor는 checkpoint layout을 위해 보존하지만 동결하고 forward에서 사용하지 않는다.
- graph-linear arm은 raw tensor 41,400개를 저장하지만 실제 trainable raw parameter는 20,700개다. paper arm은 41,400개를 모두 학습하므로 이 비교를 parameter-count 또는 capacity matched 결과로 해석하지 않는다.
- graph-linear effective scale의 평균·RMS·최솟값·최댓값과 학습된 paper arm의 nonlinear-OFF NLL을 모두 기록한다.

scRNA가 실제 dendrite morphology나 NMDA spike를 관측한다는 주장은 하지 않는다. 논문의 구조·입력 상호작용과 conductance/slope 결합을 PRISM module 계산으로 옮긴 computational analogy다.

## 2. technical-zero 경로

기존 in-batch KNN을 바로 폐기하지 않고 다음처럼 전환한다.

- epoch 1~12: canonical in-batch `k=8` 계산을 그대로 유지해 Phase-I parity를 보존한다.
- epoch 13~16: train-only bank로 0.25, 0.50, 0.75, 1.00 비율로 전환한다.
- bank neighbor: 같은 cell type·같은 region, 자기 row 제외, 기본적으로 query donor 제외, donor마다 가장 가까운 세포 한 개만 투표한다.
- technology: neighbor 제한에 사용하지 않는다.
- sex: 봉인된 sex-linked whitelist에서만 same-sex 후보를 요구한다. 후보가 부족하면 bank 목적값은 0이며, ramp 종료 뒤 cross-sex fallback은 남지 않는다.
- 일반 gene도 distinct-donor 조건이 부족하면 기본 bank 목적값을 0으로 두며, `in_batch` fallback은 명시적인 legacy ablation에서만 가능하다.
- donor×cell-type에서 두 cell을 고를 때 두 region이 있으면 region별 한 개를 우선 선택한다.
- detection rate와 depth ECDF는 donor 안에서 먼저 평균한 뒤 donor 간 평균을 내므로 donor별 질량이 같다.
- epoch-frozen bank detection sufficient statistic·donor 수·depth ECDF만 refresh 시 cache하고, live step-wise global EMA shrinkage는 매 step 적용한다.
- proven/stable 진단은 batch 평균이 아니라 epoch 점수 합과 실제 위치 수로 계산하며 cell-type·region별 위치·세포·batch·distinct-donor 분모를 함께 남긴다.
- bank refresh는 encoder-only forward를 사용하며 thinning RNG와 model train/eval 상태를 복원한다.
- checkpoint의 bit-packed bank는 현재 train row 선택 및 row별 cell type·region·donor·sex metadata와 일치해야만 resume된다.

## 3. integrated learned pathology rank

기존 SV7 learned-rank 코드를 그대로 재사용하지 않았다. 기존 설정은 Phase II pathology contribution이 0으로 내려가기 전인 epoch 1~16에 96개 rank 후보와 generator-count 학습이 겹쳤고, canonical Phase-I fixed rank 8과도 같지 않았다. 새 코드 계약은 다음과 같다.

- 모델은 epoch 1 전에 algebraic capacity 96개를 한 번만 할당한다.
- epoch 1~24에는 앞의 8개 component만 정확히 1, 나머지 88개는 정확히 0인 live mask를 사용한다. gate gradient도 정확히 0이다.
- epoch 25 첫 step부터 96개 전체의 learned soft mask를 직접 활성화한다. warmup mask와의 cross-fade는 없다. 기존 8개의 keep probability는 0.995, 새 후보는 0.5이며 extra `U`는 0 초기화다.
- epoch 29~32에는 target/minimum 없는 rank cardinality penalty를 ramp한다.
- epoch 33~36에는 exact-hard forward와 straight-through gradient를 사용한다.
- epoch 37에 threshold 0.45·0.50·0.55 count spread와 near-threshold uncertainty가 설정 한도를 넘으면 mask 동결을 fail-closed로 거부한다.
- pathology route가 rank search 중 0이면 config load가 실패한다. 새 blueprint는 final epoch까지 scale 1.0을 요구한다.
- generator search가 rank freeze epoch보다 먼저 시작하면 config load가 실패한다. 제안값은 rank freeze epoch 37, generator shadow 37~40, soft 41~44, hard 45~55다.
- strict historical checkpoint 호환성을 위해 warmup mask는 immutable config에서 재구성하며 새 persistent checkpoint key를 추가하지 않았다.
- 로그에는 live expected/hard rank와 underlying search-logit expected/hard rank, mode, temperature, uncertainty, freeze readiness가 함께 남는다.
- generator 개수 선택은 rank mask가 동결된 뒤에만 시작한다. 반대 순서는 넓은 residual generator를 먼저 줄여 좁은 named-pathology rank가 일반 residual을 대신 맡게 할 위험이 있어 채택하지 않았다.

이 구조적 freeze gate는 validation 생물학 성능을 대신하지 않는다. epoch 36 validation reconstruction·세포형별 general/AD module recovery·tail·leakage noninferiority도 별도의 commit gate로 통과해야 한다.

### 3.1 새 모델 전용 실시간 로그

학습 중 tmux에서 구조 탐색 상태를 바로 볼 수 있도록 두 층의 로그를 추가했다.

- step log: `prism_integrated` console style에서 100 micro-batch마다 pathology route scale·실제 크기, live/search rank, rank mode·temperature·sparsity, freeze low/mid/high count·spread·uncertainty·ready, generator count·mode, rank/generator overlap 여부를 출력한다.
- 같은 step log에 `■ ref 논문 적용부` block을 별도로 출력한다. 현재 arm, computational-analogy 한계, nonlinear ramp, nonlinear/local RMS, threshold crossing, nonlinear-OFF branch/full NLL gain, mix·threshold·slope·gain을 명시한다.
- epoch log: validation과 generator safety audit 뒤 `[architecture-capacity]` dashboard에 `paper branch`와 `paper OFF audit`를 함께 출력한다.
- machine-readable: `architecture_capacity_history.jsonl`, `architecture_capacity_latest.json`, `architecture_capacity_epoch_NNN.json`을 즉시 기록하며 `paper_compartmental_branch` object를 포함한다.
- fail-closed: 두 cardinality objective가 같은 epoch에 활성화되면 runtime error로 중단한다.
- logging helper는 `torch.no_grad()`와 detached tensor만 사용하며 parameter·buffer·RNG가 바뀌지 않는 회귀 테스트를 추가했다.

## 4. 실제 데이터/registry 확인

실제 414-module registry로 graph builder와 runtime loader를 실행한 결과는 다음과 같다.

- shape: `(414, 414)`
- nonempty rows: `319`
- directed edges: `1081`
- top-k: `4`
- minimum Jaccard: `0.05`
- zero diagonal 및 nonempty-row normalization: loader 통과
- smoke artifact SHA-256: `845c5ed61f5557cf3a762c01c08773020eee97e9ec5d19c02bacfad3bfbd6e9d`

이 smoke artifact는 `/tmp` 검증물이며 아직 최종 sealed launch artifact가 아니다.

현재 13,498-gene decoder 순서에서 whitelist resolution도 다시 확인했다.

| Symbol | Ensembl ID | Decoder index |
|---|---|---:|
| KDM5D | ENSG00000012817 | 190 |
| ZFY | ENSG00000067646 | 606 |
| USP9Y | ENSG00000114374 | 2471 |
| RPS4Y1 | ENSG00000129824 | 3512 |
| EIF1AX | ENSG00000173674 | 7391 |
| UTY | ENSG00000183878 | 8153 |
| RPS4X | ENSG00000198034 | 8909 |
| EIF1AY | ENSG00000198692 | 8984 |
| DDX3X | ENSG00000215301 | 9429 |
| XIST | ENSG00000229807 | 10024 |

이름·Ensembl ID·index가 blueprint 계약과 다르면 runtime이 중단된다. 다만 최종 sex-linked biological annotation 범위는 launch 전에 별도 검토해야 한다.

## 5. 실행한 검증

SV7의 실제 PyTorch environment에서 저장소 전체 테스트를 실행했다. 이 workspace 환경에는 `pytest` package가 없어 동일한 `unittest.TestCase` 전체를 표준 discovery로 실행했다.

```bash
cd /home/kmlee/sv7/project_local/prism
PYTHONPATH=src:. /home/kmlee/sv7/miniconda3/envs/rie_bam/bin/python -m unittest discover -s tests -t .
```

2026-08-28 direct rank-first 전환과 paper-branch 실시간 logging까지 반영한 뒤 전체 결과는 `55 tests passed in 21.0s`다. rank·logging·compartmental targeted 회귀는 합계 `22 tests passed`다. 경고는 기존 Transformer nested-tensor 조건에 관한 PyTorch warning이며 실패는 없었다.

검증 범위에는 다음이 포함된다.

- nonlinear scale 0에서 v1과 `torch.equal` identity
- threshold innovation의 원점 값과 1차 기울기 0
- 두 neighbor drive의 supralinear threshold crossing
- graph-linear arm의 threshold·slope 동결, `mix·gain`만의 20,700 trainable raw parameter와 linear homogeneity
- stacked-input Jacobian의 graph-linear rank 48 대 threshold rank 96
- 414 modules × 24 cell types에서 paper arm parameter 수 41,400, trainable dense graph 없음
- graph artifact order/hash/split/normalization fail-closed
- same-cell-type·region/distinct-donor KNN 및 wrong-region 차단
- donor-balanced detection rate/depth ECDF
- 바뀌는 global EMA 여러 값에서 cached query와 reference loop의 float32 parity 및 cache 재사용
- sex-linked Phase-I parity와 full-bank fail-closed 전환
- 일반 gene의 bank-invalid 기본 zero fallback과 explicit legacy in-batch fallback
- canonical in-batch `k=8` parity와 bank `k=16` 분리
- bit-packed bank checkpoint round-trip 및 현재 train-row metadata 검증
- bank refresh의 encoder-only 실행, thinning state와 model mode 복원
- position-weighted epoch proven/stable 평균, support count, cell-type·region별 bank coverage
- 기존 v1 module-local 및 integrated config contract 회귀
- Phase-I fixed-rank8 live mask와 gate zero-gradient
- warmup rank8에서 epoch 25 full soft mask로의 직접 전환과 hard/frozen checkpoint round-trip
- threshold-sensitive rank mask의 fail-closed freeze
- Phase-II pathology route 0 설정 거부
- rank freeze 이전 generator search 설정 거부
- 기존 learned-rank config의 load compatibility
- architecture snapshot의 JSON 직렬화와 parameter·buffer 무변경
- rank/generator cardinality overlap 표시
- `prism_integrated` step log의 구조 용량 block 렌더링
- `prism_integrated` step log와 epoch JSON/dashboard의 explicit paper branch·OFF counterfactual 렌더링

변경 파일 Python byte-compile과 v2 blueprint JSON parsing도 통과했다.

## 6. 아직 통과하지 않은 launch gate

다음은 코드 unit/component 검증으로 대신할 수 없으며 아직 완료되지 않았다.

1. 최종 train-only module-local reliability artifact 생성 및 봉인
2. 최종 compartment graph artifact 경로와 SHA-256 봉인
3. canonical Phase-I epoch-12 full-system trajectory parity와 live rank 8 확인
4. 실제 full-data system-build와 epoch 1~55 merged curriculum traversal smoke
5. multi-GPU DDP forward/backward/checkpoint smoke
6. optimizer·scheduler·curriculum·rank gate·technical-zero bank를 포함한 end-to-end resume equivalence
7. exact optimizer-step budget audit
8. 기존에 사용된 9개 legacy test donor의 selection 배제와 새 untouched confirmatory holdout 봉인
9. epoch 36 rank threshold stability 및 validation noninferiority commit audit
10. rank freeze 전 generator objective zero-step audit
11. 세 function-class arm과 paper fixed-rank8 control의 같은 seed·split·step-budget 재학습
12. validation reconstruction, 세포별 general/AD module recovery, sex/tech/pathology leakage 비교
13. thinning proven-dropout AUROC/AP 또는 proven/stable gap 개선
14. sex-linked whitelist의 최종 biological annotation review
15. 동일 장비·batch에서 기준 0.57초/step 대비 bank 경로 처리량 회귀 5% 이내 확인
16. 명시적 launch 승인

따라서 blueprint의 `_launch_guard.launch_allowed`는 `false`다. 기존 SV6/SV7 학습, optimizer, checkpoint 및 tmux session은 변경하지 않았다.

## 7. 검토할 파일

- 설계: `docs/prism_integrated_compartmental_nonlinearity_techzero_v2_design_20260827.md`
- review blueprint: `configs/final/prism_integrated_compartmental_nonlinearity_techzero_v2_blueprint_20260827.json`
- nonlinear model: `src/kmlee_bam/model/precision_medicine.py`
- graph builder/loader: `scripts/training/build_prism_module_local_compartment_graph_20260827.py`, `src/kmlee_bam/data/module_local_reliability.py`
- technical-zero bank: `src/kmlee_bam/objectives/latent_knn_pi_tech.py`
- trainer/runtime wiring: `src/kmlee_bam/training/adaptive_subgroup_trainer.py`, `src/kmlee_bam/training/run_current.py`, `src/kmlee_bam/training/runner_base.py`
- learned rank: `src/kmlee_bam/training/learned_pathology_rank.py`, `src/kmlee_bam/training/learned_generator_count.py`
- live architecture logging: `src/kmlee_bam/training/architecture_capacity_logging.py`
- regression tests: `tests/test_module_local_compartmental_nonlinearity.py`, `tests/test_latent_knn_pi_tech_bank.py`, `tests/test_running_averages_sparse_metrics.py`, `tests/test_integrated_learned_pathology_rank_curriculum.py`, `tests/test_architecture_capacity_logging.py`
