# 통합 PRISM v2: compartmental module-local 비선형성과 train-only technical-zero bank

작성일: 2026-08-27  
수정일: 2026-08-28 — 통합 learned pathology-rank curriculum 추가  
상태: 코드 구현 및 로컬/CPU 회귀 검증 완료, **학습 미실행·검토 전용**

## 1. 결론

기존 `module-local residual v1` 코드는 이미 존재했다. 이번 구현은 그것을 폐기하거나 별도 모델로 만드는 것이 아니라, 같은 통합 PRISM 안에 다음 두 가지를 기본값 OFF인 선택 기능으로 추가한다.

1. 모듈을 고정된 계산 구획(compartment)으로 보고, 인접 모듈의 공동 활성에 대해 threshold·slope·gain을 학습하는 비선형 module-local 가지
2. 현재 mini-batch KNN을 같은 세포형·donor-balanced·train-only latent bank로 점진적으로 교체하는 technical-zero 추정기

현재 실행 중인 학습의 모델·optimizer·scheduler·checkpoint는 변경하지 않았다. 이 설계는 다음 통합 학습 후보를 위한 것이다.

## 2. 논문에서 직접 가져온 사실과 PRISM으로 옮긴 해석

참조 논문은 Aizenbud 외, *Dendritic morphology and synaptic nonlinearities enhance functional complexity in human cortical neurons*, PNAS 2026, DOI 10.1073/pnas.2533168123이다.

논문이 직접 보인 핵심은 다음과 같다.

- 단순한 분지 개수보다 전체 수상돌기 면적과 분지 구조의 조합이 FCI를 더 잘 설명했다. 논문에서 단일 `total dendritic area`의 설명력은 약 0.74, `number of bifurcation branches`는 약 0.29였고, 면적과 가장 긴 bifurcation branch의 조합은 약 0.81이었다.
- 더 큰 NMDA conductance 또는 더 가파른 voltage-dependence 하나만 바꾼 hybrid보다 두 성질을 함께 가진 human-type synapse에서 비선형성과 FCI 차이가 더 뚜렷했다.
- 특정 모델에서는 약 35개의 동시 활성 synapse 부근에서 sublinear 합산이 supralinear 합산으로 전환됐지만, 이것은 그 생물물리 모델의 수치이지 scRNA 모듈의 보편 임계값이 아니다.

PRISM에서의 `module = dendritic compartment` 대응은 **계산적 비유**이다. scRNA 데이터가 실제 수상돌기 모양, NMDA spike 또는 synapse 수를 관측한다고 주장하지 않는다. 따라서 논문의 35라는 값을 복사하지 않고, 모듈별 임계값과 기울기를 제한된 범위 안에서 학습한다.

## 3. 기존 v1 경로를 보존하는 방법

기존 v1은 target 세포형의 관측을 완전히 제외한 같은 donor의 다른 세포형 정보로 모듈 요약을 만든다. donor를 $d$, target 세포형을 $c$, 모듈을 $m$이라 할 때 기존 선형 변환 결과를 다음처럼 둔다.

$$u_{dcm} = \operatorname{LinearLocal}_{cm}\!\left(s_{dm}^{(-c)}\right)$$

여기서 $s_{dm}^{(-c)}$에는 target 세포형의 두 region 값이 모두 들어가지 않는다. 기존 personal rank-2 span도 ridge projection으로 제거된다.

새 비선형 기능이 꺼져 있거나 curriculum 계수가 0이면 출력 계산은 기존 v1 경로를 그대로 사용한다.

$$\lambda_{NL} = 0 \Longrightarrow u'_{dcm} = u_{dcm}$$

이 조건은 단순 근사가 아니라 회귀 테스트에서 `torch.equal`로 확인한다. 기본 OFF 구성은 새 parameter와 persistent buffer를 만들지 않는다.

## 4. 고정 compartment graph

모듈 registry의 binary gene membership만 사용해 모듈 간 Jaccard overlap을 계산한다.

$$J_{mn} = \frac{|G_m \cap G_n|}{|G_m \cup G_n|}$$

각 모듈에서 overlap이 큰 이웃 최대 4개만 남기고, 최소 Jaccard는 0.05로 두며, 자기 연결은 제거한다. 남은 행은 합이 1이 되도록 고정한다.

$$A_{mm} = 0,\qquad A_{mn} \geq 0,\qquad \sum_n A_{mn} = 1$$

고립 모듈의 행은 전부 0이다. 이 graph는 donor·cell·pathology·validation·test 결과를 전혀 열지 않고 만들어지며, registry SHA-256과 정확한 모듈 순서를 loader가 검증한다. $A$는 학습되지 않으므로 자유로운 $M \times M$ 행렬이 생물학적 모듈을 재배치하지 못한다.

## 5. 비선형 compartment 경로

인접 모듈 drive와 기존 모듈 drive를 섞는다.

$$h_{dcm} = \sum_n A_{mn}u_{dcn},\qquad x_{dcm} = u_{dcm} + \mu_{cm}h_{dcm}$$

mix, threshold, slope, gain은 global 값과 평균 0인 cell-type delta의 합으로부터 bounded sigmoid로 만든다.

$$\mu_{cm} = \mu_{\max}\sigma(a^{\mu}_m+\Delta a^{\mu}_{cm})$$

$$\theta_{cm} = \theta_{\min}+(\theta_{\max}-\theta_{\min})\sigma(a^{\theta}_m+\Delta a^{\theta}_{cm})$$

$$\beta_{cm} = \beta_{\min}+(\beta_{\max}-\beta_{\min})\sigma(a^{\beta}_m+\Delta a^{\beta}_{cm})$$

$$g_{cm} = g_{\max}\sigma(a^{g}_m+\Delta a^{g}_{cm})$$

현재 review blueprint의 범위는 $\mu_{\max}=0.5$, $\theta_{\min}=0.5$, $\theta_{\max}=2.0$, $\beta_{\min}=1.0$, $\beta_{\max}=8.0$, $g_{\max}=1.0$이다.

비선형 innovation은 다음과 같다.

$$r_{dcm} = g_{cm}x_{dcm}\left[\sigma\!\left(\beta_{cm}(|x_{dcm}|-\theta_{cm})\right)-\sigma(-\beta_{cm}\theta_{cm})\right]$$

두 번째 sigmoid를 빼는 이유는 원점에서 값뿐 아니라 1차 기울기도 0으로 만들기 위해서다.

$$r_{dcm}(0) = 0,\qquad \left.\frac{\partial r_{dcm}}{\partial x_{dcm}}\right|_{x_{dcm}=0} = 0$$

따라서 새 가지는 기존 선형 adapter를 단순히 다른 gain으로 복제하기보다 threshold·curvature가 필요한 부분을 맡아야 한다. 최종 local 입력과 bounded 출력은 다음과 같다.

$$u'_{dcm} = u_{dcm}+\lambda_{NL}r_{dcm}$$

$$b_{dcm} = C\tanh(q_{cm}u'_{dcm})$$

그 뒤 기존과 동일하게 personal rank-2 span을 투영해 제거하고, 전체 module-local curriculum 계수 $\lambda_L$을 곱한다. 새 경로가 rank-2 또는 생성자의 일을 무제한으로 빼앗지 못하게 하는 장치는 다음과 같다.

- 고정 sparse graph
- bounded mix·threshold·slope·gain
- 원점 1차 기울기 0
- 기존 rank-2 span 제거
- v1보다 늦은 ramp
- 전체 local 크기 제한과 cell-type hierarchy penalty 유지
- nonlinear-OFF counterfactual NLL을 매 epoch 기록

모듈 수가 $M=414$, 세포형 수가 $C=24$이면 새 학습 parameter 수는 다음과 같다.

$$N_{NL} = 4M(C+1) = 41{,}400$$

학습 가능한 dense $414 \times 414$ 행렬은 추가하지 않는다.

### 5.1 graph-linear 함수족 대조군

이 비교는 parameter-count 또는 effective-capacity matched 실험이 아니다. 현재 선형 함수족에서 threshold·slope·gain을 하나의 scale로 합치면 여러 raw parameter가 같은 함수로 퇴화하고, AdamW와 hierarchy penalty가 paper arm과 다른 최적화 기하를 만든다. 따라서 그런 과잉파라미터화는 대조군으로 사용하지 않는다.

`graph_linear_control`은 paper arm과 같은 고정 graph·입력·curriculum을 사용하되 `mix`와 `gain`만 학습한다. checkpoint tensor layout을 일정하게 유지하기 위해 threshold·slope tensor는 저장하지만 `requires_grad = False`로 동결하고 forward와 hierarchy regularizer에서 사용하지 않는다.

$$r_{dcm}^{lin} = g_{cm}x_{dcm}$$

모듈 수가 $M=414$, 세포형 수가 $C=24$일 때 graph-linear arm은 raw tensor 41,400개를 저장하지만 실제 학습 raw parameter는 20,700개다. cell-type delta의 평균 0 제약을 반영한 generic effective-field 상한은 graph-linear arm이 19,872개이고 paper threshold arm이 39,744개다. 따라서 결론은 **같은 graph와 입력에서 선형 함수족과 threshold-curvature 함수족 중 무엇이 더 유용한가**로만 해석한다.

검증은 같은 seed·split·step budget·graph·ramp와 공유된 mix·gain 초기값으로 다음 세 arm을 각각 처음부터 연속 학습해 비교한다.

1. v1 linear reference: 논문 기능 이전 경로와의 총 차이
2. graph-linear function-class control: 고정 graph와 선형 혼합의 효과
3. paper compartmental threshold: bounded threshold·slope·gain과 curvature의 추가 효과

graph-linear arm에서는 effective linear scale의 평균·RMS·최솟값·최댓값을 매 epoch 기록한다. 작은 회귀 테스트는 여러 입력을 쌓은 Jacobian을 사용하며, 모듈 8개·세포형 3개·세포형당 최소 4개의 다양한 입력에서 graph-linear arm의 generic rank를 48, threshold arm의 generic rank를 96으로 확인한다. 이 Jacobian은 국소 식별 가능성 회귀 테스트이지 전체 표현력의 증명은 아니다.

학습된 paper arm에서 nonlinear contribution만 빼는 counterfactual NLL도 유지하지만, 이것은 별도로 재학습한 graph-linear 함수족 대조군을 대신하지 않는다.

## 6. technical-zero bank가 decoder의 sex·tech 항과 별도로 필요한 이유

decoder의 sex·tech 항은 **예측값을 분해**한다. 반면 $q_{ig}^{tech0}$는 관측된 0이라는 **학습 label을 얼마나 믿을지** 정한다. decoder에 nuisance 항이 있어도 잘못된 hard-zero gradient가 먼저 들어가면 label reliability 문제는 자동으로 해결되지 않는다. 그러므로 두 기능은 중복이 아니다.

다만 nuisance를 KNN 제한에 똑같이 적용해서도 안 된다.

- cell type은 유전자의 정상 ON/OFF 문맥을 크게 바꾸므로 반드시 같은 cell type으로 제한한다.
- brain region은 알려진 생물학적 문맥이고 PRISM decoder의 별도 항으로 설명되므로 DLPFC와 MTG도 섞지 않는다.
- technology는 서로 다른 assay에서의 검출 차이 자체가 technical zero의 증거가 될 수 있으므로 제한하지 않는다.
- sex는 전체 후보를 절반으로 줄이지 않고, 사전에 봉인한 sex-linked gene whitelist에만 같은-sex 이웃을 사용한다.

## 7. train-only same-celltype·region donor bank

각 epoch 경계에서 train split의 각 donor·cell-type 조합으로부터 최대 2개 세포를 고정 seed로 선택한다. 한 조합에 DLPFC와 MTG가 모두 있으면 먼저 region별로 한 개씩 선택한다. validation, 이미 사용된 legacy test donor, 새 confirmatory holdout은 bank 생성에 사용하지 않는다.

cell-type·region별 gene detection rate와 depth ECDF는 선택된 세포를 바로 평균하지 않는다. 먼저 donor 안에서 평균한 뒤 donor 평균을 다시 평균하여, 희귀 cell type에서 한 donor가 한 세포만 제공하고 다른 donor가 두 세포를 제공해도 donor당 질량이 정확히 같게 한다. bank 갱신 forward는 $z_i^{clean}$만 필요하므로 ordinal decoder와 auxiliary head는 실행하지 않는다.

query 세포 $i$에 대해 후보는 다음 조건을 모두 만족해야 한다.

1. 같은 cell type
2. 같은 brain region
3. 자기 자신과 다른 row
4. 기본값에서는 query donor와 다른 donor
5. donor마다 가장 가까운 세포 하나만 투표

그 뒤 가까운 donor $k$개의 검출 여부를 평균한다.

$$n_{ig}^{bank} = \frac{1}{K_i}\sum_{j\in\mathcal{N}_i^{bank}}\mathbf{1}(o_{jg}>0)$$

$K_i$가 최소 4개 donor보다 작으면 일반 gene과 sex-linked gene 모두 bank 목적값을 보수적으로 0으로 둔다. 희소 cell type·region에서 구 in-batch 경로가 조용히 되살아나는 것을 막는 fail-closed 기본값이다. `bank_invalid_fallback = in_batch`는 명시적인 legacy ablation에서만 허용한다. bank의 기본 $k$는 16~32 후보이며, review blueprint는 $k=16$을 제안한다. 다만 canonical Phase-I parity를 위해 epoch 1~12의 in-batch 값은 기존 $k=8$을 그대로 쓴다.

`depth_proxy = n_detected`일 때 low-depth percentile도 current mini-batch 전체가 아니라 같은 cell-type·region bank의 donor-balanced depth 분포에서 계산한다. gene detection rate도 같은 cell type·region의 bank rate를 global EMA 쪽으로 shrink한다.

$$D_{crg}^{shrunk} = \omega_{cr}D_{crg}^{bank}+(1-\omega_{cr})D_g^{EMA},\qquad \omega_{cr} = \frac{N_{cr}^{donor}}{N_{cr}^{donor}+8}$$

처리량을 위해 epoch bank에서만 결정되는 donor-balanced $D_{crg}^{bank}$, donor 수, depth ECDF를 bank refresh 때 cache한다. 최종 shrink 결과를 고정하지는 않는다. 매 training step에서 현재의 live global EMA인 $D_g^{EMA}$를 cached bank 통계와 다시 섞으므로 기존 step-wise EMA 의미가 유지된다. 서로 다른 EMA 벡터를 넣은 reference-loop 대 cached-vectorized parity test를 필수로 둔다.

sex-linked gene에서는 이웃 검출률과 detection rate를 모두 같은 sex·같은 cell type으로 계산한다. sex가 unknown이거나 같은-sex donor가 부족하면 해당 gene의 점수는 보수적으로 다음처럼 둔다.

review blueprint의 10개 후보 symbol은 현재 13,498-gene decoder의 `feature_name`과 Ensembl ID/index에 모두 exact resolution됐으며, 기대 Ensembl ID와 decoder index도 함께 봉인한다. decoder에 없는 후보는 실제 whitelist에서 제외했다. 다만 10개 후보의 최종 생물학적 annotation 범위는 launch 전에 한 번 더 검토해야 한다.

$$q_{ig}^{tech0} = 0$$

기존 in-batch KNN은 sex 조건을 모르지만 canonical Phase-I parity를 깨지 않기 위해 epoch 1~12에는 원래 값을 그대로 유지한다. epoch 13~16에는 whitelist 유전자를 기존 in-batch 값에서 same-sex bank 값으로 점진 전환하며, same-sex donor support가 부족하면 전환 목적값을 0으로 둔다. 따라서 ramp가 끝난 뒤에는 부적절한 cross-sex fallback이 남지 않는다.

technology별 제한 또는 technology별 detection rate는 사용하지 않는다.

## 8. 기존 추정기와의 점진적 연결

bank를 갑자기 전면 적용하지 않고 기존 mini-batch 점수와 섞는다.

$$q_{ig}^{tech0} = (1-\eta_e)q_{ig}^{batch}+\eta_eq_{ig}^{bank}$$

제안 일정에서는 epoch 13~16에 $\eta_e=0.25,0.50,0.75,1.00$으로 증가한다. bank 후보가 부족하면 모든 gene에서 $q_{ig}^{bank}=0$으로 둔다. 따라서 ramp 초반에는 기존 in-batch 기여가 점진적으로 남지만 ramp 종료 뒤에는 invalid query가 완전히 fail-closed 된다.

$q_{ig}^{tech0}$는 계속 detached 값이며, reconstruction weight·zero/nonzero BCE weight·soft target의 기존 세 경로만 조절한다. module-local 또는 비선형 path의 입력으로 넣지 않는다.

bank의 latent·row·cell type·region·donor·sex와 bit-packed detection mask는 checkpoint에 저장된다. epoch refresh는 dataset의 thinning RNG를 진행시키지 않도록 thinned view를 잠시 끄고 별도 generator를 쓴 뒤 상태를 원복한다.

resume 시에는 저장된 bank를 곧바로 신뢰하지 않는다. 현재 train dataset과 고정 seed로 다시 계산한 선택 row 집합이 정확히 같은지 확인하고, 각 row의 cell type·region·donor·sex metadata와 decoder gene 수가 모두 일치해야만 복구한다.

## 9. 하나의 통합 run 안에서의 curriculum

이 일정은 architecture를 갈아 끼우거나 optimizer를 재시작하는 Phase 전환이 아니다. 한 모델·한 optimizer/scheduler trajectory·한 run directory·한 checkpoint history 안의 curriculum이다. technical-zero와 module-local 경로에는 수치적 ramp가 남지만, rank 후보 mask에는 warmup-to-search cross-fade를 사용하지 않는다.

### 9.1 먼저 고쳐야 했던 구조적 문제

기존 SV7 learned-rank 설정은 epoch 13~16에 pathology decoder를 0으로 cross-fade하고 epoch 17에 rank mask를 동결했다. 따라서 후반에 보이는 rank 수는 최종 reconstruction에 계속 기여하는 경로의 용량이 아니었다. pathology 경로가 0인 뒤에는 task gradient가 rank 필요성을 검증할 수 없고 cardinality 압력만 남을 수 있다. 새 설계에서는 `phase1_pathology_scale = 1.0`과 `phase2_pathology_scale = 1.0`을 사용해 named-pathology decoder와 auxiliary tie를 마지막 epoch까지 활성 상태로 유지한다.

한편 epoch 1부터 96개 후보를 모두 학습하면 canonical Phase-I rank-8 warmup과 같지 않고, technical-zero 및 논문 적용 compartmental 경로가 아직 만들어지기 전에 rank가 그 미완성 residual을 대신 떠맡을 수 있다. 새 gate는 모델을 처음부터 최대 96개 후보로 한 번만 만들되 epoch 1~24에는 앞의 8개만 정확히 1, 나머지는 정확히 0인 mask를 사용한다. 이 기간 gate logit의 gradient도 정확히 0이다. checkpoint 이식이나 optimizer 재생성 없이 동일 모델이 fixed rank 8처럼 작동한다.

### 9.2 rank와 generator의 역할 및 순서

여기서 학습하는 rank는 personal rank 2나 Lie-generator 내부 rank가 아니다. four-axis pathology head가 만든 Braak·Thal·LATE·Lewy 좌표에 따라 세포형별 gene correction을 표현하는 `path_single_U`와 `path_V`의 공유 basis 개수다. 반면 generator는 전체 residual latent에 따른 일반 세포상태 변형을 담당한다. 두 경로 모두 gene score를 바꿀 수 있어 완전히 식별 가능하지 않으므로, 동시에 개수 벌점을 주면 한쪽이 줄어드는 동안 다른 쪽이 그 역할을 받아 성능은 유지되는 것처럼 보일 수 있다.

그래서 개수 학습은 직렬화한다.

1. 모든 generator를 켠 상태에서 pathology rank를 먼저 학습한다.
2. rank mask를 검증하고 완전히 동결한다.
3. 그 다음에만 generator shadow·soft·hard cardinality 학습을 시작한다.

두 cardinality penalty가 같은 epoch에 활성화되는 구성은 runtime config validation이 거부한다. generator search에는 `require_frozen_pathology_rank_before_search = true`를 요구한다.

반대 순서인 generator-first는 사용하지 않는다. generator는 전체 residual state를 담는 더 넓은 경로이므로 먼저 줄이면 named-pathology rank가 일반 residual까지 대신 담아 rank가 부풀 수 있다. 반대로 rank path는 Braak·Thal·LATE·Lewy 입력과 auxiliary supervision에 묶인 더 좁은 경로이므로 이를 먼저 확정한 뒤 넓은 generator를 압축하는 편이 역할 혼동이 작다. 여기서 “rank 먼저, generator 나중”은 **개수 선택**의 순서다. 한 통합 모델의 일반 reconstruction·encoder·decoder weight 학습은 계속 이어지며 optimizer나 checkpoint lineage를 끊지 않는다.

| Epoch | 활성화 내용 |
|---:|---|
| 1~12 | canonical CT64 biological warmup을 그대로 유지한다. four-axis first-order pathology decoder와 donor-balanced auxiliary head를 실제로 학습한다. 최대 96개 tensor는 이미 존재하지만 live mask는 정확히 rank 8이다. |
| 13~16 | rank 8을 유지한 채 train-only technical-zero bank를 0.25씩 ramp한다. |
| 17~24 | rank 8을 유지한 채 target-celltype-excluded module-local v1을 8 epoch 동안 ramp한다. compartmental nonlinearity는 epoch 21~24에 4 epoch ramp한다. |
| 25~28 | 모든 generator를 live all-on으로 유지한다. epoch 25 첫 step부터 rank 후보 96개의 learned soft mask를 **한 번에 직접 활성화**한다. warmup mask와의 cross-fade는 없다. task가 새 basis를 발견하게 하되 cardinality penalty는 아직 0이다. 기존 8개 logit은 keep probability 0.995, 새 88개는 중립적인 0.5에서 시작하고 extra `U`는 0으로 시작한다. 따라서 새 88개가 경계에서 출력을 튀게 하지는 않으면서 첫 step부터 gradient를 받는다. |
| 29~32 | rank soft mask를 유지하면서 target 없는 normalized cardinality penalty를 4 epoch 동안 ramp한다. 수동 target rank와 minimum rank는 두지 않는다. |
| 33~36 | rank는 exact hard forward와 straight-through gradient로 안정화한다. generator는 계속 all-on이며 generator-count objective는 0이다. |
| 37~40 | epoch 36 rank commit audit를 통과한 mask만 동결한다. 그 뒤 generator shadow search를 시작하지만 live decoder의 generator는 모두 켠다. |
| 41~44 | rank를 고정한 채 generator soft adaptation을 수행한다. |
| 45~55 | rank를 고정한 채 generator exact-hard constrained search와 consolidation을 수행한다. module rescue·AD-axis rescue·tail protection은 late stage 내내 유지한다. |

rank commit은 단순히 epoch 37이 됐다는 이유만으로 과학적 승인되지 않는다. epoch 36에서 threshold 0.45·0.50·0.55의 count spread가 4 이하이고, 0.45~0.55에 놓인 후보 비율이 0.15 이하이며, validation reconstruction·세포형별 general/AD module recovery·tail floor·sex/tech/depth/pathology leakage가 비열등해야 한다. 실패하면 generator search와 model release를 시작하지 않는다. 이는 rank 8을 보호하는 최소값이 아니라, 경계 부근의 불확실한 mask를 확정하지 않기 위한 안정성 gate다.

세 arm은 같은 schedule을 공유하므로 arm 간 최종 비교는 가능하지만, epoch 21~23의 nonlinear-OFF NLL은 아직 ramp 중인 v1 기저에 대한 값이다. 논문 가설의 주효과와 arm 간 성능 결론은 두 ramp가 모두 완료된 epoch 24 이후 값만 사용한다.

모든 paper-function arm은 같은 learned-rank 정책을 사용한다. 다만 비선형 효과와 rank-learning 효과를 분리하기 위해 `paper_compartmental_threshold + fixed pathology rank 8` 대조군을 별도로 한 번 더 학습한다. 이 대조군과 learned-rank primary는 seed·split·optimizer-step budget·graph·비선형 schedule이 같아야 한다.

## 10. 새 로그와 진단

비선형 경로는 다음을 기록한다.

- paper branch enabled 여부와 `compartmental_threshold`/`graph_linear_control` arm 이름
- nonlinear module RMS
- nonlinear/local RMS 비율
- threshold crossing fraction
- nonlinear ramp
- nonlinear-OFF branch NLL과 full NLL
- nonlinear-ON 대비 NLL gain
- 학습된 mix·threshold·slope·gain의 요약값
- graph-linear arm의 effective linear scale 평균·RMS·최솟값·최댓값

technical-zero bank는 다음을 기록한다.

- bank blend
- 유효 query 비율
- 평균 distinct-donor $K_i$
- bank epoch
- sex-linked query 유효 비율
- thinning proven-dropout과 stable-zero에서의 점수 및 gap
- cell-type·region별 bank valid fraction과 평균 distinct-donor 수
- cell-type·region별 proven/stable 점수 합·위치 수·세포 수·probe batch 수·distinct donor 수

pathology rank와 generator는 다음을 매 epoch 시작과 종료에 함께 기록한다.

- rank capacity, live expected rank, live hard rank, search-logit expected/hard rank
- fixed-rank warmup·soft·hard-ST·frozen mode
- rank keep-probability uncertainty와 temperature
- axis별 correction participation rank·95% energy rank·singular-value spectrum
- rank sparsity multiplier와 finalized 여부
- generator live hard count·expected count·uncertainty
- generator all-on·shadow·soft·hard mode

### 10.1 실시간 확인 형식

새 run은 `KMLEE_CONSOLE_LOG_STYLE = prism_integrated`, `train.log_every = 100`, `train.progress_log_every = 500`, `train.eval_log_every = 2000`을 사용한다. 따라서 tmux log를 보고 있으면 100 training micro-batch마다 다음 구조 용량 block이 나타난다.

```text
■ 구조 용량 자동학습 — pathology rank를 확정한 뒤 generator를 줄임
  Pathology rank live hard/E ....... 40/42.0 / 96
  rank gate 온도/희소화 ............ 0.800/0.50
  rank 확정 low/mid/high ........... 42/40/39 · spread 3 · 경계후보 3% · ready YES
  개수학습 분리 검사 ............... rank penalty ON · generator objective OFF · overlap NO
```

같은 주기로 논문 적용부도 별도 block으로 출력한다. 숫자만 남기지 않고 출처 가설, 현재 arm, computational analogy라는 한계, 활성 ramp, 비선형 크기, 문턱 통과율, OFF counterfactual NLL 이득을 함께 표시한다.

```text
■ ref 논문 적용부 — dendritic compartment·문턱·초선형성의 계산적 비유
  arm/활성도 ........................ paper_compartmental_threshold · ramp 1.00
  비선형/기본 local RMS ............. 0.01000/0.02000 · 비율 0.500 · 문턱통과 20.00%
  이 경로를 끄면 생기는 ΔNLL ....... branch/full +0.00300/+0.00200 (양수면 도움)
  학습된 mix/threshold/slope/gain .. 0.100/1.200/4.000/0.400
```

각 epoch validation과 generator safety audit가 끝나면 `[architecture-capacity epoch NNN/055]` dashboard를 출력한다. 동시에 다음 파일을 즉시 갱신한다.

- `architecture_capacity_history.jsonl`: epoch별 append-only 구조 용량 이력
- `architecture_capacity_latest.json`: 현재 최신 상태
- `architecture_capacity_epoch_NNN.json`: epoch별 고정 snapshot
- 기존 `pathology_correction_rank.jsonl`: axis별 participation rank·95% energy rank·spectrum
- 기존 `joint_generator_count_latest.json`과 epoch별 파일: generator mask와 SHA-256

로그는 모두 detached/no-grad 진단이다. parameter·optimizer·scheduler·RNG를 바꾸지 않는다. pathology route가 rank search 중 0이 되거나, rank freeze가 threshold-sensitive하거나, rank와 generator cardinality가 동시에 활성화되면 경고만 남기지 않고 fail-closed로 중단한다.

기존 `RunningAverages`의 sparse-key 수정은 metric이 없는 step에 의해 평균이 희석되는 문제만 해결한다. proven/stable 진단에는 그것을 사용하지 않는다. 각 GPU에서 epoch 동안 점수 합과 실제 cell-gene 위치 수를 누적한 뒤 DDP 합산하여 평균과 gap을 계산한다. 위치 수뿐 아니라 세포 수와 distinct donor 수도 함께 기록하며, 최소 위치 256개·세포 8개·donor 2명 중 하나라도 미달하면 subgroup 점수와 gap은 기록하지 않는다. 이때 누락은 0이 아니라 unsupported/NA를 뜻하고 raw 분모는 항상 남긴다.

## 11. 검증 gate

학습 launch 전에 다음이 모두 필요하다.

1. canonical Phase-I epoch-12 parity와 live pathology correction rank 8 확인
2. extra 88개 rank 후보가 Phase I forward와 gate gradient에 정확히 0으로 기여하고, epoch 25 첫 step에서 cross-fade 없이 전체 learned mask가 직접 적용되는지 확인. extra `U`는 zero-init이므로 첫 backward에서는 `U`가 먼저 gradient를 받고, gate logit은 `U`가 0에서 벗어난 뒤 task gradient를 받는 것이 정상이다.
3. rank gate·optimizer·frozen mask와 bank를 포함한 resume-equivalence 및 checkpoint round-trip
4. merged curriculum smoke traversal
5. exact optimizer-step budget audit
6. 기존 9개 test donor를 model selection에서 완전히 제외하고 새 untouched confirmatory holdout을 봉인했는지 확인
7. pathology route가 rank search 전체와 최종 epoch에서 0이 아닌지 확인
8. epoch 36 rank threshold-sensitivity·uncertainty·validation noninferiority commit gate 통과
9. rank mask가 동결되기 전에 generator objective가 한 step도 실행되지 않는지 확인
10. graph registry SHA·module order·zero diagonal·row normalization 검증
11. sex-linked whitelist의 decoder gene-order exact resolution
12. thinning proven-dropout 기준으로 bank가 기존 in-batch보다 AUROC·AP 또는 proven/stable gap을 개선하는지 확인
13. fixed rank-2 epoch-20 baseline 대비 reconstruction·세포별 module recovery·AD module recovery·sex/tech/pathology leakage 비열등성 확인
14. `graph_linear_control` 재학습 함수족 대조군 대비 validation NLL·module recovery의 이득이 실제로 양수인지 확인하고, parameter-count matched 결과로 해석하지 않기
15. 같은 paper arm 안에서 nonlinear-OFF counterfactual validation NLL gain도 양수인지 확인
16. learned-rank paper arm과 fixed-rank8 paper arm을 같은 budget으로 비교
17. nonlinear/local RMS와 threshold crossing이 0 또는 포화 상태로 붕괴하지 않는지 확인
18. stacked-input Jacobian 회귀에서 graph-linear와 threshold arm의 generic effective-field rank가 각각 48과 96인지 확인
19. cached bank 통계가 여러 live EMA 값에서 reference loop와 float32 허용오차 안에 일치하는지 확인
20. 기준 0.57초/step 대비 bank 경로의 wall-time 증가가 5% 이내인지 동일 장비·동일 batch 조건으로 확인

기존 9개 test donor는 이미 사후분석에 사용됐으므로 이 새 실험의 untouched test로 다시 부를 수 없다. 이 donor들은 참고용 legacy test로만 보고 model·epoch·rank·generator 선택에는 쓰지 않는다. 최종 확인에는 새로 봉인한 외부 donor holdout이 필요하다.

## 12. 구현 위치

- `src/kmlee_bam/model/precision_medicine.py`: bounded compartmental nonlinearity와 curriculum
- `src/kmlee_bam/data/module_local_reliability.py`: graph artifact loader와 provenance 검증
- `scripts/training/build_prism_module_local_compartment_graph_20260827.py`: registry-only graph builder
- `src/kmlee_bam/objectives/latent_knn_pi_tech.py`: train-only bank, same-celltype·region/distinct-donor KNN, sex-linked safety, checkpoint packing
- `src/kmlee_bam/data/ordinal_dataset.py`: runtime split identity를 보존해 bank가 train split인지 fail-closed 검증
- `src/kmlee_bam/training/adaptive_subgroup_trainer.py`: bank refresh·blend·checkpoint 연결
- `src/kmlee_bam/training/core_trainer.py`: nonlinear counterfactual 진단과 sparse-metric 분모 수정
- `src/kmlee_bam/training/run_current.py`: config·train-dataset·sex-linked whitelist 연결
- `src/kmlee_bam/training/prism_module_rescue_training.py`: nonlinear ramp 완료 후 rescue allowlist 연결
- `src/kmlee_bam/training/learned_pathology_rank.py`: Phase-I fixed-rank8 mask, 이후 soft/hard/frozen rank selection
- `src/kmlee_bam/training/runner_base.py`: active pathology route와 rank-before-generator ordering fail-closed 검증

이 문서와 blueprint는 launch 승인이 아니다. 실제 full config 생성·artifact 봉인·학습 실행은 별도 검토와 필수 launch gate 이후에만 수행한다.
