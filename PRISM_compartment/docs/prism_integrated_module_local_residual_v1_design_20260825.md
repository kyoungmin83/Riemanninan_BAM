# 통합 PRISM module-local personal residual v1 설계

상태: **코드 구현 완료, 사용자 검증·smoke·학습 미승인**  
작성일: 2026-08-25  
목표: 기존 rank-2 personal code가 보존하는 donor 전역 특성은 유지하면서, target cell type을 전혀 보지 않고도 module별 donor 편차를 복원한다.

> 2026-08-27 addendum: v1의 정확한 no-op 경계와 target-celltype firewall을 유지한 채, 별도 default-OFF compartmental nonlinearity와 train-only technical-zero bank를 구현했다. 최신 설계와 검증 계약은 `docs/prism_integrated_compartmental_nonlinearity_techzero_v2_design_20260827.md`를 따른다. v1 문서는 당시 선형 residual의 독립적인 범위를 보존하기 위해 그대로 둔다.

범위: 이 문서의 v1은 논문 아이디어를 넣기 전에 만든 **선형 비교 baseline**이다. 따라서 v1 단독 경로에는 NMDA에서 착안한 학습 가능 threshold, gain–slope 결합, supralinear synergy가 없다. 현재 진행하는 논문 적용 실험은 v1을 대조군으로 보존한 채, 이 기능들을 추가한 `compartmental nonlinearity + technical-zero v2`이다. 즉 “논문 아이디어를 적용하지 않는다”는 뜻이 아니라, **v1 결과를 논문 기반 비선형 v2의 결과인 것처럼 해석하지 말자**는 뜻이다.

## 1. 확정 결정

- 최종 모델 형식은 하나의 architecture, 하나의 연속 run, 하나의 optimizer/checkpoint lineage를 갖는 통합 PRISM이다.
- module-local residual parameter는 epoch 1 이전에 모두 생성한다. 다만 기존 PRISM 경로가 안정되기 전에는 출력 scale을 0으로 유지하고 curriculum으로 연다.
- local adapter의 첫 검증 rank는 8로 고정한다. adaptive rank는 local residual 자체의 효과가 입증된 뒤 별도 요인으로 시험한다.
- 현재 SV6 fixed pathology-rank 8 run과 SV7 learned pathology-rank run은 변경하지 않는다. 두 run 중 validation 계약을 통과한 rank policy를 다음 통합 run의 matched comparator와 module-local arm에 똑같이 적용한다.
- 위 문장은 살아 있는 run을 변경하지 않는다는 뜻이다. 다음 module-local run에서 기존 승인 계약인 `hard pruning from epoch 25`를 늦추는 것은 별도 설계 변경이며 아직 승인되지 않았다.
- 공식 test donor 9명은 계속 봉인한다. rank policy 선택, module reliability 선택, hyperparameter 선택에 test 결과를 사용하지 않는다.

## 2. 해결할 실제 병목

현재 personal 경로는 target cell type의 모든 region context를 제외한 donor support를 Transformer로 통합한 뒤 rank-2 code로 압축하고, 이를 다시 414개 module coefficient로 복원한다. 이 경로는 donor의 전역적 개인차를 요약하는 데 적합하지만 특정 module만 donor별로 움직이는 국소 편차를 잃을 수 있다.

새 경로는 기존 rank-2 personal path를 대체하지 않는다. 두 경로의 역할은 다음과 같이 고정한다.

- rank-2 personal path: donor의 전역적인 개인 특성
- module-local residual: target-free support에서 보존되는 module별 국소 개인 편차

## 3. 입력과 leakage 계약

source context index를 `q`, source cell type을 `c(q)`, source region을 `r(q)`, target cell type을 `t`, module을 `m`이라고 한다. local residual의 입력 weight는 다음과 같다.

$$w_{d,t,q,m}=\mathbf{1}[\mathrm{observed}_{d,q}]\mathbf{1}[c(q)\ne t]r^{\mathrm{count}}_{d,q}r^{\mathrm{split}}_{q,m}$$

- `r_count`는 현재 context cell-count reliability다.
- `r_split`은 train-donor-only split-half module reliability다.
- target cell type의 DLPFC와 MTG context는 모두 weight가 정확히 0이어야 한다.
- donor ID는 support table을 찾는 index로만 쓰며 trainable donor embedding을 만들지 않는다.
- pathology label, age, target region, target nucleus latent, BAM uncertainty는 module-local input으로 사용하지 않는다.
- validation/test donor에서는 그 donor의 다른 cell type support를 사용할 수 있지만, module reliability prior와 모든 adapter parameter는 train donor에서만 학습한다.

### 3.1 Region-balanced module summary

cell 수가 많은 region이 local summary를 지배하지 않도록 먼저 region 안에서 pooling하고, 관측 가능한 region을 동일 질량으로 합친다.

$$S_{d,t,r,m}=\frac{\sum_{q:r(q)=r}w_{d,t,q,m}X_{d,q,m}}{\sum_{q:r(q)=r}w_{d,t,q,m}+\epsilon}$$

$$I_{d,t,r,m}=\mathbf{1}\!\left[\sum_{q:r(q)=r}w_{d,t,q,m}>0\right],\qquad \widetilde S_{d,t,m}=\frac{\sum_r I_{d,t,r,m}S_{d,t,r,m}}{\sum_r I_{d,t,r,m}+\epsilon}$$

`X`는 train-only pathology/nuisance residualizer를 통과해 context별 scale로 표준화된 기존 `source_module`이다. 입력은 outlier 하나가 low-rank path 전체를 흔들지 않도록 train-only 고정 한계로 clip한다.

### 3.2 새 reliability artifact

현재 `source_reliability`는 donor×context 스칼라이므로 module-local 경로에는 부족하다. 다음 artifact를 새로 만든다.

`module_local_context_reliability_train64_s42_v1.npz`

필수 내용:

- `context_module_reliability`: shape `[n_contexts, 414]`
- `context_module_reliable_mask`: shape `[n_contexts, 414]`
- deterministic cell split seed와 split rule
- train donor allow-list와 hash
- module/context name과 순서
- registry, activity dictionary, source context artifact SHA-256
- validation donor와 test donor를 reliability fit에 사용하지 않았다는 증거

각 donor×cell-type×region의 세포를 donor ID와 cell ID의 deterministic hash로 양분하고, train donor 사이의 두 half pseudobulk module activity 재현성을 사용한다. reliability 산출은 validation 성능을 보고 반복 조정하지 않는다.

## 4. Adapter architecture

module 수를 `M=414`, 고정 local rank를 `R=8`로 둔다. 읽기 행렬 `V`는 모든 target cell type이 공유하고, 쓰기 행렬 `U_t`와 diagonal scale `d_t`는 global parameter에 shrink된 target-cell-type deviation을 갖는다.

$$\ell_{d,t}=V\widetilde S_{d,t},\qquad h_{d,t}=d_t\odot\widetilde S_{d,t}+U_t\ell_{d,t}$$

$$d_t=d_0+\delta d_t,\qquad U_t=U_0+\Delta U_t,\qquad \sum_t\delta d_t=0,\qquad \sum_t\Delta U_t=0$$

여기서 `V`의 shape은 `[8,414]`, `U_t`의 shape은 `[414,8]`, `d_t`의 shape은 `[414]`이다. 완전한 414×414 dense transform은 금지한다.

bounded raw residual은 다음과 같다.

$$\widehat{\Delta a}_{d,t}^{\mathrm{local}}=a_{\max}\tanh\!\left(g_t\odot h_{d,t}\right),\qquad g_t=0\ \text{at initialization}$$

`a_max`는 임의 상수가 아니라 frozen rank-2 epoch-20 comparator의 **train-donor personal coefficient 절대값 99.5% 분위수**에서 사전 고정한 cap이다. `build_prism_module_local_output_cap_20260825.py`가 sealed checkpoint SHA-256, source-config SHA-256, rank 2, train-donor allow-list와 source-context SHA-256을 검증한 뒤 별도 artifact로 저장하며, runtime은 이 artifact 없이는 fail closed한다. 실제 e20/train64 calibration 결과는 `a_max=2.3915109324455255`, coefficient 수는 635,904개이고 artifact는 `provenance/module_local/rank2_e20_train64_output_cap_q995_20260825.npz`, SHA-256은 `76edbc17c78a29318302468fc0cf045c88cfe9b58b147f2ec57af673bcd8666d`다. validation/test coefficient로 cap을 조정하지 않는다. `g_t=0`이므로 새 parameter가 존재해도 초기 forward는 정확한 no-op이다.

adapter 본체의 parameter 수는 nuisance probe를 제외하면 약 106,398개다. 414×414 dense matrix 하나보다 작고, 모든 저차원 교환은 고정 rank 8을 통과한다.

### 4.1 rank-2 personal path와의 역할 분리

새 adapter가 기존 rank-2 personal basis를 복제해 rank-2 경로를 장식물로 만들지 못하도록 raw residual의 current personal-basis span 성분을 **ridge-stabilized 방식으로 억제한다**. target cell type `t`의 기존 personal basis를 `B_t`라고 한다.

$$B_t=\operatorname{stopgrad}(\mathrm{personal\_basis}_t),\qquad P_t=B_t^\top(B_tB_t^\top+\epsilon I_2)^{-1}B_t$$

$$\Delta a_{d,t}^{\mathrm{local}}=(I_{414}-P_t)\widehat{\Delta a}_{d,t}^{\mathrm{local}}$$

`stopgrad`는 local path가 기존 global personal subspace 자체를 움직여 projection을 우회하지 못하게 한다. ridge가 양수이므로 이는 exact orthogonal removal이 아니며, singular direction별 잔여 비율은 0이 아니라 ridge와 singular value에 의해 정해진다. personal basis의 수치 rank가 낮은 초기 구간에는 이 안정성이 필요하다.

구현에서는 `P_t` 또는 `I_414`를 실제 414×414 tensor로 만들지 않는다. residual을 `B_t`에 투영하고 2×2 linear solve를 거쳐 다시 module space로 보내는 동등한 rank-2 계산만 사용한다.

### 4.2 최종 personal coefficient

$$a_{d,t}^{\mathrm{personal,new}}=a_{d,t}^{\mathrm{personal,rank2}}+\lambda_{\mathrm{local}}(e)\Delta a_{d,t}^{\mathrm{local}}$$

나머지 normal-region, age, common pathology, donor-by-pathology response 항은 변경하지 않는다. 새 coefficient는 기존 fixed module activity dictionary로 gene score에 lift되며 full decoder와 target-latent-free branch NLL에 동일하게 들어간다.

## 5. 초기화와 checkpoint 호환성

- `g_t`: 정확히 0으로 초기화
- `d_0`: 1로 초기화
- `delta d_t`: 0으로 초기화
- `V`: deterministic semi-orthogonal initialization
- `U_0`: 작은 deterministic normal initialization
- `Delta U_t`: 0으로 초기화
- local curriculum scale: 0으로 초기화하고 checkpoint에 저장

`g_t=0`이면 첫 backward에서 주로 `g_t`가 먼저 학습되고 내부 `U/V/d` gradient는 0일 수 있다. 이는 의도된 exact no-op warm start다. `g_t`가 열린 다음 step부터 내부 parameter가 학습된다. graph와 reliability는 학습 전에 결정되어 있으므로 이 한-step 지연은 문제가 아니다.

기존 checkpoint를 새 architecture에 warm-start할 때 허용되는 missing key는 `precision_head.module_local_*` prefix로만 제한한다. 그 외 missing/unexpected key가 있으면 fail closed한다. optimizer, scheduler, PHU, EMA, generator mask, curriculum stage는 새로운 통합 run의 checkpoint에 함께 저장한다.

## 6. 하나의 연속 curriculum

새 parameter는 처음부터 존재하지만 loss/output strength만 순차적으로 연다.

| Epoch | module-local 상태 | 다른 핵심 상태 |
|---|---|---|
| 1–12 | instantiated, exact output 0 | 기존 biological warmup |
| 13–16 | output 0 | common/personal/response PRISM crossfade 안정화 |
| 17–24 | `lambda_local`을 0에서 1로 완만히 증가 | rank-2 personal path와 local role audit |
| 25–35 | full local residual, module rescue gradient 허용 | **미승인 제안:** generator는 shadow/soft search만 허용 |
| 36–43 | local consolidation | **미승인 제안:** PHU 및 protected hardening |
| 44–55 | local parameter를 유지하며 마지막 안정화 | **미승인 제안:** biological guard를 통과한 generator cardinality만 확정 |

이는 Phase I model과 Phase II model을 따로 학습하는 방식이 아니다. architecture, optimizer, scheduler, checkpoint history는 처음부터 끝까지 하나다. 다만 epoch 25 이후의 표는 기존 승인된 generator 계약을 변경하는 제안일 뿐 확정 schedule이 아니다. 다음 launch 전에 `epoch 25 hard pruning 유지`와 `local 안정화 후 지연` 중 하나를 사용자가 명시적으로 결정해야 하며, 총 epoch 55도 같은 결정에 포함된다.

## 7. Loss와 안전장치

기존 full reconstruction NLL과 target-latent-free branch NLL이 주 학습 신호다. 다음 보조항을 train split에서만 사용한다.

### 7.1 donor-balanced zero center

local path가 common disease 또는 normal-region offset을 흡수하지 않도록 target cell type별 donor 평균을 0 근처에 둔다.

$$L_{\mathrm{local,center}}=\frac{1}{T}\sum_t\left\|\frac{\sum_d\omega_d\Delta a_{d,t}^{\mathrm{local}}}{\sum_d\omega_d+\epsilon}\right\|_2^2$$

### 7.2 hierarchical shrinkage

$$L_{\mathrm{local,hier}}=\frac{1}{T}\sum_t\left(\|\delta d_t\|_2^2+\|\Delta U_t\|_F^2\right)$$

### 7.3 contribution cap

local residual이 기존 personal path 전체를 압도하면 role separation이 실패한 것이다. 초기에는 local-to-rank2 RMS ratio에 soft cap을 걸고, validation에서 cap에 계속 붙는 경우만 train-only calibration으로 완화한다.

$$L_{\mathrm{local,size}}=\left[\frac{\operatorname{RMS}(\lambda_{\mathrm{local}}\Delta a^{\mathrm{local}})}{\operatorname{RMS}(a^{\mathrm{personal,rank2}})+\epsilon}-\rho_{\max}\right]_+^2$$

### 7.4 nuisance leakage guard

`local_code=V*S_tilde`에 작은 gradient-reversal pathology/age probe를 둔다. 이 probe는 local baseline이 named pathology와 age를 다시 운반하지 못하게 할 뿐, target nucleus나 validation label을 입력으로 사용하지 않는다.

### 7.5 Cell-type/module-family graph prior — optional secondary arm

같은 donor의 target cell type 관측값은 절대 쓰지 않는다. training-only registry에서 정한 DLPFC–MTG module-family edge에 한해 parameter-level graph shrinkage를 optional secondary arm으로 둔다. e19 interim에서는 Astrocyte가 최저 recovery이자 rank-2 대비 최대 하락이고, Microglia-PVM과 Sst Chodl도 하락한 반면 Oligodendrocyte는 최대 개선이었다. 따라서 Oligodendrocyte 전용 prior를 선험적 primary 목표로 두지 않는다. 다음 matched comparator에서 deficit을 재확인한 뒤 graph prior의 대상 cell type을 고정한다.

$$L_{\mathrm{family,graph}}^{(t)}=\sum_{(m_1,m_2)\in E_t}w_{m_1m_2}\left\|\Theta_{t,m_1}-\Theta_{t,m_2}\right\|_2^2$$

primary local-residual ablation에서는 이 항을 끄고, local residual 자체가 유효하며 train-only deficit 근거가 재현된 경우에만 켜서 효과를 분리한다.

## 8. Optimizer와 module-rescue 계약

- 기존 `precision_head`와 local adapter를 별도 optimizer group으로 분리한다.
- 첫 제안값은 local adapter learning-rate multiplier 10이다. 실제 값은 smoke에서 finite-gradient와 output-open 속도만 보고 고정하며 validation 성능으로 반복 탐색하지 않는다.
- `prism_module_rescue` parameter allow-list에 local adapter를 명시적으로 추가한다.
- local residual이 완전히 열린 epoch 25 이후에만 module-rescue gradient가 `module_local_*` parameter로 들어갈 수 있다.
- generator hard-pruning 지연은 합리적 가설이지만 아직 승인되지 않았다. 기존 계약은 epoch 25 hard pruning이며, 변경하려면 별도 사용자 결정을 기록한다.

## 9. 필수 로그

매 epoch와 curriculum boundary에 다음을 저장한다.

- local residual RMS, rank-2 personal RMS, 두 값의 ratio
- local residual module coverage와 zero-coverage fraction
- local code의 singular spectrum, participation rank, 95% energy rank
- ridge 억제 후 local residual의 personal-basis row별 정규화 overlap 평균과 최댓값
- target cell type별 local residual RMS와 module recovery 기여
- local-on 대 local-off paired full NLL과 branch NLL
- local-on 대 local-off 일반/AD module Spearman
- pathology/age leakage probe 성능
- rank-2 personal effective dimension과 ablation contribution
- direct, Lie, affine, explicit PRISM, local residual의 RMS와 ablation contribution
- Astrocyte, Microglia-PVM, Sst Chodl recovery와 Oligodendrocyte 개선 보존 여부

RMS만으로 route bypass를 판단하지 않는다. paired ablation에 따른 NLL과 module recovery 변화도 함께 기록한다.

## 10. 필수 검증

### 10.1 Exact identity

이전 checkpoint를 새 architecture에 load하고 `g_t=0`, `lambda_local=0`일 때 score, probability, NLL이 기존 forward와 정확히 같아야 한다.

### 10.2 Target perturbation firewall

- source buffer에서 target cell type의 두 region 값을 임의로 크게 바꿔도 local output은 정확히 변하지 않아야 한다.
- non-target support를 바꾸면 local output은 변해야 한다.
- target nucleus latent와 BAM uncertainty를 바꿔도 local output은 변하지 않아야 한다.

### 10.3 Reliability provenance

- reliability artifact fit donor가 정확히 train donor allow-list와 같아야 한다.
- validation/test donor가 fit에 하나라도 포함되면 launch를 차단한다.
- module/context name과 registry 순서가 다르면 launch를 차단한다.

### 10.4 Gradient and resume

- output-open 첫 step에는 `g_t`에 finite nonzero gradient가 있어야 한다.
- 다음 step부터 `U/V/d`에 finite nonzero gradient가 있어야 한다.
- DDP rank별 parameter와 curriculum state가 일치해야 한다.
- checkpoint resume 후 local scale, optimizer moments, PHU/EMA, generator state가 정확히 복원되어야 한다.

## 11. 비교 실험

공식 primary 비교는 동일 seed, donor split, optimizer-update budget, rank policy를 사용한다.

1. matched integrated comparator: module-local OFF
2. diagonal-only adapter
3. fixed-rank8 low-rank-only adapter
4. proposed diagonal + fixed-rank8 adapter
5. 4번이 통과하고 train-only deficit이 재현된 뒤에만 optional cell-type graph prior

최종 test는 validation에서 architecture와 epoch가 고정된 뒤 한 번만 연다. 현재 진행 중인 SV6 fixed-rank8와 SV7 learned-rank 결과는 module-local architecture 선택에 사용하지 않고, 다음 matched pair에서 공통으로 사용할 pathology-rank policy만 선택한다.

## 12. 합격 조건

- 기존 global reconstruction과 leakage gate를 모두 통과한다.
- 전체 및 AD module recovery가 matched comparator보다 개선되거나 사전 정의된 non-inferiority 범위 안에 있다.
- cell-type 25분위와 reliable worst-cell regression guard를 통과한다.
- e19 최악/하락 cell type인 Astrocyte, Microglia-PVM, Sst Chodl은 명시적으로 보고하고, Oligodendrocyte의 기존 개선도 donor-level bootstrap confidence interval과 함께 보존 여부를 보고한다.
- cell-type leakage balanced accuracy는 matched integrated comparator보다 개선되어야 하며 frozen rank-2 epoch-20 대비 사전 정의 margin을 통과해야 한다. module recovery만 좋아지고 leakage가 남으면 replacement-ready로 판정하지 않는다.
- rank-2 personal path의 effective dimension과 ablation contribution이 붕괴하지 않는다.
- local path가 pathology/age를 baseline 경로로 우회 운반하지 않는다.
- direct route 증가만으로 성능이 오른 경우는 module-local 성공으로 인정하지 않는다.

## 13. 구현 touchpoint

- `src/kmlee_bam/model/precision_medicine.py`: config, target-free pooling, adapter, output, penalties
- `src/kmlee_bam/training/runner_base.py`: reliability/output-cap artifact load와 sealed-donor 검증
- `src/kmlee_bam/training/core_trainer.py`: local diagnostics와 paired ablation
- `src/kmlee_bam/training/run_current.py`: audited warm-start missing-key prefix와 resume state
- `src/kmlee_bam/training/prism_module_rescue_training.py`: local parameter allow-list와 stage gate
- `scripts/training/build_prism_module_local_reliability_20260825.py`: train-only split-half reliability artifact
- `scripts/training/build_prism_module_local_output_cap_20260825.py`: frozen rank-2 e20 train64 output-cap artifact
- unit/contract tests: identity, target firewall, gradient opening, normalized span overlap, artifact provenance, optimizer/rescue/warm-start 계약
- launch 전 별도 integration tests: 실제 DDP consistency와 checkpoint resume

## 14. 구현 및 launch 상태

2026-08-25에 다음 구현을 완료했다.

- default-off `PrecisionMedicineConfig`와 target-excluded module-local adapter
- train-only split-half reliability artifact builder와 fail-closed loader
- frozen rank-2 epoch-20 train-donor personal coefficient에서 output cap을 만드는 sealed builder와 fail-closed loader
- full/branch local-on 대 local-off paired NLL 진단
- local regularizer, effective-rank 진단, 별도 optimizer group
- epoch 17–24 curriculum state와 epoch 25 rescue allow-list gate
- 기존 precision checkpoint에서 `precision_head.module_local_*`만 허용하는 strict warm-start 감사
- support-query 음성 donor에도 동일한 module-local counterfactual을 적용하는 대칭성 보장

정적 compile과 격리된 SV7 CPU 환경의 14개 unit/contract test가 통과했고, 기존 통합 rank-8 설정의 module-local default-off load와 실제 frozen rank-2 e20/train64 output-cap calibration도 통과했다. 이는 blueprint의 DDP 및 실제 resume test가 통과했다는 뜻이 아니다. 실제 reliability artifact 생성, full-data system build, DDP smoke/resume, 기존 checkpoint exact score audit, posthoc module-recovery traversal과 generator curriculum 변경 승인은 사용자 검증 단계로 남아 있다. 이 검증과 별도 launch 승인 전에는 새 training launch를 허용하지 않는다.
