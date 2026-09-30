# PRISM 다음 학습 설계: projection-aware whitening과 module-rescue v2

- 작성일: 2026-08-04
- 상태: 코드 구현 및 독립 검토 반영 완료, 본학습은 아직 시작하지 않음
- 대상 모델: PRISM, Lie generator 414개 고정
- 목적: 생성자 수 선택과 섞지 않고, 잠재공간 수식을 바로잡고 donor-centered 모듈 복원도를 먼저 개선한다.

## 1. 먼저 확정할 결론

이번 단계에는 두 변화만 다룬다.

1. 성별 방향을 제거한 뒤 whitening 목표를 전체 32차원 단위행렬이 아니라 실제로 남은 부분공간에 맞춘다.
2. 현재 module-rescue가 최종 평가값인 donor-centered 모듈 복원도를 더 직접 배우도록 손실, 표본추출, 업데이트 범위를 고친다.

Lie generator는 414개를 전부 켜고 고정한다. 자동 generator 선택, pathology interaction, 새로운 decoder 용량 변경은 넣지 않는다. 모듈 복원과 생성자 수를 동시에 바꾸면 어떤 변화가 성능을 만들었는지 알 수 없기 때문이다.

## 2. 현재 구현을 감사해 확인한 사실

현재 module-rescue에는 이미 다음 기능이 있다.

- 24개 세포형을 빠짐없이 순환한다.
- 같은 세포형 안에서 여러 donor를 동시에 비교한다.
- 예측값과 관측값을 각각 donor-center한다.
- 관측 목표는 뽑힌 2개 세포가 아니라, train donor의 적격 세포 전체로 미리 만든 고정 pseudobulk다.
- target-cell latent를 쓰지 않는 PRISM branch view와 최종 full-decoder view를 함께 학습한다.
- 고정 module activity 행렬 $W^{act}$를 사용해 예측 유전자 상태와 관측 유전자 상태를 동일한 모듈 좌표로 옮긴다.

따라서 다음 학습에서 “donor-centering을 새로 추가한다”거나 “24개 세포형을 처음으로 모두 사용한다”고 설명하면 사실과 다르다. 개선 대상은 이미 존재하는 rescue 신호의 정확도와 학습 효율이다.

현재 sv6 설정은 GPU 6개에서 rank당 8번, 즉 전역적으로 epoch당 48개 세포형 contrast block을 처리한다. 24개 세포형은 각각 두 번 등장하며, 한 block에서는 donor마다 세포 2개를 예측에 사용한다. sv7 설정도 GPU 수에 맞춰 전역 48 block이 되도록 구성됐다.

## 3. 설계 1: projection-aware whitening

### 3.1 문제

정제 전 잠재상태를 $z_i^{\perp}\in\mathbb R^{32}$라 하고, 실제로 제거된 성별 방향의 직교기저를 $Q^{sex}\in\mathbb R^{32\times r_{sex}}$라 한다. 정제된 상태는 $z_i^{clean}=z_i^{\perp}P$이고, 여기서 $P=I_{32}-Q^{sex}(Q^{sex})^{\top}$이다.

실제 제거 rank가 $r_{sex}$이면 남은 유효차원은 $d_{eff}=32-r_{sex}$이다. 네 방향이 모두 유효하면 $d_{eff}=28$이다. $z^{clean}$은 제거된 네 방향에서 값이 0이어야 하므로, 공분산 $C=\operatorname{Cov}(z^{clean})$를 $I_{32}$에 맞추는 기존 목표에는 없앨 수 없는 상수 오차가 남는다.

### 3.2 새 목표

단일변수 비교 실험에서는 새 whitening loss를 $\mathcal L_{white}^{proj}=\|C-P\|_F^2/32^2$로 정의한다. 분모를 우선 기존과 같은 $32^2$로 유지해 loss 크기 변화가 결과를 교란하지 않게 한다. $d_{eff}^2$ 분모를 쓰는 버전은 별도 실험으로만 검토한다.

### 3.3 중요한 해석

$Q^{sex}$가 detach되어 있고 $z^{clean}$이 정확히 $P$ 공간에 있다면 기존 loss와 새 loss의 차이는 거의 상수 $r_{sex}/32^2$다. 같은 분모를 쓰면 gradient도 이론적으로 동일하다. 따라서 이 수정은 먼저 수식, 로그, checkpoint 비교를 정확하게 만드는 작업이며, 이것 하나가 모듈 복원도를 크게 높일 것이라고 주장하지 않는다.

### 3.4 코드 변경 계약

1. nuisance projector가 그 forward에서 실제 사용한 detached $Q^{sex}$, $P$, $r_{sex}$, $d_{eff}$를 제공한다.
2. whitening 함수는 선택적 인수 `whitening_target`을 받는다. 값이 없으면 기존처럼 $I_{32}$를 써 backward compatibility를 유지한다.
3. projection이 꺼졌거나 유효 방향이 없으면 $P=I_{32}$가 되어 기존 loss와 정확히 같아야 한다.
4. rescue 전용 forward는 nuisance EMA 통계를 갱신하지 않는다. 한 번에 한 세포형만 담긴 rescue batch가 성별 제거축을 편향시키지 않게 한다.
5. DDP의 모든 rank가 같은 $Q^{sex}$와 $P$를 사용하는지 hash와 최대 원소 차이로 확인한다.

### 3.5 추가 로그

- 실제 제거 rank $r_{sex}$
- 남은 유효차원 $d_{eff}=32-r_{sex}$
- active-space whitening error $\|C-P\|_F^2/32^2$
- removed-space leakage $\|(I-P)C(I-P)\|_F^2$
- active/removed cross leakage $\|PC(I-P)\|_F^2$
- 기존 identity-target loss와 새 projector-target loss의 차이

### 3.6 필수 시험

1. projection OFF에서 새 loss와 기존 loss 및 gradient가 일치해야 한다.
2. rank 4이고 $C=P$인 합성 자료에서 새 loss는 0에 가까워야 한다.
3. 같은 $32^2$ 분모와 detached $P$를 쓸 때 기존 loss와 새 loss의 gradient가 수치 허용오차 안에서 같아야 한다.
4. 유효 rank가 0, 1, 2, 3, 4로 바뀌어도 shape와 loss가 정상이어야 한다.
5. checkpoint 저장·재시작 뒤 $Q^{sex}$와 $P$가 재현되어야 한다.

## 4. 설계 2: module-rescue v2

### 4.1 현재 결과가 말하는 정확한 병목

같은 pooled val+test 비교에서 CT64와 PRISM recovery epoch 6은 다음 패턴을 보였다.

| 지표 | CT64 | PRISM recovery | 해석 |
|---|---:|---:|---|
| pooled module recovery | 0.8458 | 0.6935 | PRISM이 낮음 |
| donor-centered blind recovery | 0.8127 | 0.6875 | 가장 큰 개선 대상 |
| cross-disjoint recovery | 0.5693 | 0.6304 | PRISM이 오히려 높음 |
| blind raw recovery | 0.9793 | 0.9804 | 세포형 평균 구조는 이미 잘 맞음 |

기존 PRISM에서 module-rescue를 켜자 donor-centered recovery는 0.6611에서 0.6875로 실제 상승했다. 즉 방향은 맞지만 힘이 부족했다. raw recovery는 이미 CT64와 같고 cross-disjoint는 더 좋으므로, 모델 전체를 크게 바꾸기보다 donor별 미세한 모듈 차이를 더 안정적으로 배우게 해야 한다.

### 4.2 그대로 보존할 것

- train-only donor-context target artifact
- donor-balanced sampling
- 예측과 관측의 별도 donor-centering
- train-only robust scale
- branch/full dual view
- 24개 세포형 균등 노출
- singleton 10개를 포함한 414개 모듈 경로
- test donor를 학습이나 checkpoint 선택에 사용하지 않는 원칙

### 4.3 개선 A: 더 많은 고유 세포를 보되 GPU 메모리는 늘리지 않는다

현재 한 contrast block은 donor마다 세포 2개를 사용한다. 이것은 donor의 모든 세포를 버린다는 뜻이 아니다. 관측 target은 이미 그 donor의 적격 세포 전체로 만들지만, 그 target을 향해 gradient를 계산할 때 예측 쪽 표본이 2개라 흔들릴 수 있다는 뜻이다.

한 번에 4~8개를 올리면 과거 smoke에서 GPU 메모리가 위험했다. 따라서 한 forward의 크기는 donor당 2개로 유지하고, epoch 안에서 서로 겹치지 않는 독립 draw 수를 2회에서 4회로 늘린다.

- 전역 contrast block: 48개/epoch에서 96개/epoch로 증가
- 각 세포형: 2 draw/epoch에서 4 draw/epoch로 증가
- 각 draw: donor당 2개 유지
- 결과: 가능한 context에서는 donor당 epoch마다 최대 8개의 서로 다른 세포를 본다.
- 세포가 8개보다 적으면 모두 한 번씩 본 뒤에만 deterministic reshuffle을 허용한다.

이 방식은 한 번에 드는 GPU 메모리를 늘리지 않고 Monte Carlo 잡음을 줄인다. 구현은 서버의 GPU 수가 달라도 전역 block 수뿐 아니라 optimizer update 수도 같도록 한다. 12개 전역 block을 한 update로 묶으므로 GPU 6개에서는 rank당 2 block, GPU 4개에서는 rank당 3 block을 누적한 뒤 update한다. 두 서버 모두 epoch당 전역 96 block과 optimizer update 8회를 사용한다.

### 4.4 개선 B: 최종 평가와 같은 방향의 loss를 추가한다

현재 학습 loss는 donor-centered 예측과 관측 사이의 원소별 Huber 오차다. 이것은 값의 차이를 줄이는 데 좋지만, 최종 평가는 주로 “donor가 높고 낮은 순서를 같이 맞혔는가”를 Spearman 상관으로 본다. 현재 Pearson과 CCC는 계산되지만 detach되어 보고용일 뿐 학습에는 쓰이지 않는다.

view $v\in\{branch,full\}$, 세포형 $t$, donor $d$, 모듈 $m$의 raw donor-centered 값을 $\widehat Z^{v}_{dtm}$과 $Z^{obs}_{dtm}$라 한다. 기존 Huber는 train-scale standardized 값에서 유지한다. 추가 concordance 항은 실제 `blind_centered` 평가와 같은 배열 구조를 사용하기 위해 세포형별 donor×module 표를 펼쳐 한 개의 differentiable CCC를 계산한다. 즉 module별 CCC를 평균하지 않는다.

새 view loss 후보는 $\mathcal L_v=\mathcal L_{Huber,v}+\lambda_{CCC}(1-\operatorname{mean}_{t}\operatorname{CCC}^{v}_{t,flat})$다. exact Spearman은 미분하기 어려우므로 CCC를 진폭에도 민감한 학습 대리값으로 쓰고, checkpoint 선택은 계속 raw donor-centered exact Spearman으로 한다.

$\lambda_{CCC}$를 결과에 맞춰 임의로 고르지 않는다. train-only 고정 calibration block에서 Huber와 CCC의 decoder gradient norm을 재고 두 항의 초기 gradient 크기가 같아지는 계수를 계산한 뒤, 미리 정한 범위 $[0.1,1.0]$로 clip하고 본학습 동안 고정한다. 계산값과 hash를 manifest에 저장한다. validation이나 test 성능으로 이 값을 조정하지 않는다.

목표 분산이 사실상 0인 세포형×모듈은 CCC 분모가 불안정하므로 train-only variance mask로 CCC 항에서만 제외하고 Huber 항에는 남긴다. 이 mask는 prediction과 무관하게 고정되며 branch와 full이 정확히 같은 mask를 사용한다. 예측이 상수로 붕괴한 경우에는 제외하지 않고 $\operatorname{CCC}=0$, loss $=1$로 처리해 복원 gradient가 생기게 한다.

### 4.5 개선 C: full view를 우선 개선하되 PRISM branch를 보호한다

현재 loss는 branch 50%, full 50%다. 결과상 cross-disjoint branch 성질은 이미 CT64보다 좋고, 최종 full output의 donor-centered module recovery가 낮다. 따라서 v2의 주 학습 혼합은 branch 25%, full 75%로 두는 것이 합리적이다: $\mathcal L_{rescue}=0.25\mathcal L_{branch}+0.75\mathcal L_{full}$.

그러나 branch를 희생해 full만 올리는 checkpoint는 채택하지 않는다. validation에서 branch cross-disjoint가 기준 PRISM보다 0.02를 초과해 하락하면 탈락시킨다. 즉 full view는 주 개선 대상이고 branch view는 비열등성 안전장치다.

이 25:75 변화는 CCC 추가와 분리해 짧은 train-only gradient audit에서 다음 세 조건을 먼저 확인한다.

1. full-view gradient norm이 실제로 증가한다.
2. branch-view gradient가 0으로 사라지지 않는다.
3. Lie, translation, direct, explicit PRISM 경로 중 하나로 전부 우회하지 않는다.

### 4.6 개선 D: rescue step이 모델 전체를 흔들지 못하게 한다

현재 extra rescue optimizer step은 사실상 모델 전체에 gradient를 줄 수 있다. donor별 작은 모듈 오차를 고치려다 encoder, nuisance axis, cell-type baseline, ordinal threshold까지 함께 움직이면 이미 잘 잡힌 누출 제어와 기하가 흔들릴 수 있다.

ordinary ELBO/NLL step에서는 기존처럼 전체 모델을 학습한다. rescue 전용 step에서만 gradient를 다음 경로로 제한한다.

- Lie action의 $u$, $v$, coefficient head
- masked affine translation $a$
- dense direct-state head
- explicit PRISM common/personal/response module-to-gene head

score mixer는 base, technical, state gate를 하나의 softmax로 함께 조절해 state 전용 파라미터만 분리할 수 없으므로 rescue step에서는 고정한다.

rescue 전용 step에서는 다음을 고정한다.

- gene/module tokenizer와 transformer encoder
- posterior mean/scale head와 cell-type prior
- sex projection basis 및 관련 EMA 통계
- cell-type, assay, sex baseline
- ordinal threshold

같은 optimizer의 Adam 상태를 유지하되, rescue 대상이 아닌 파라미터의 gradient를 `None`으로 만든 뒤 step한다. rescue step에서는 AdamW weight decay를 일시적으로 0으로 두고 step-wise scheduler를 진행하지 않아 보조손실과 무관한 추가 decay 및 학습률 일정 가속을 막는다. 적용된 update가 예정된 8회보다 적거나 AMP/non-finite 문제로 하나라도 skip되면 v2 학습은 즉시 실패 처리한다.

grouped sampler가 이미 donor마다 같은 수의 세포를 공급하므로 rescue forward에서는 ordinary training의 inverse-frequency `donor_balance_weight`를 모두 1로 덮어쓴다. 이로써 full view에만 donor 균형이 이중 적용되는 것을 막고 branch:full $=0.25:0.75$가 실제 gradient에서도 같은 의미를 갖게 한다.

rescue forward는 posterior sample 대신 posterior mean을 사용하고, 가능한 stochastic BAM/Dropout을 명시적으로 끈 결정론적 보조 경로를 사용한다. ordinary training은 기존의 reparameterized sample과 stochastic regularization을 그대로 유지한다. 이는 평가가 posterior mean으로 수행된다는 조건과 rescue gradient를 맞추고 표본 잡음을 줄이기 위한 것이다.

### 4.7 개선 E: 세포형 균형을 검사 가능한 계약으로 만든다

24개 세포형 순환은 이미 있으므로 새 알고리즘을 만들지 않는다. 대신 다음을 assertion과 로그로 고정한다.

- 매 epoch 각 세포형의 전역 block 수가 정확히 4회인지 확인
- 각 세포형에서 본 고유 donor 수와 고유 세포 수 기록
- donor별 반복 표집률과 가장 적게 본 donor 기록
- rescue-ineligible 세포형이 생기면 즉시 경고하고 ordinary training support에서는 제외하지 않음
- per-celltype Huber, CCC, Pearson을 branch/full로 분리해 기록

특정 세포형이 어렵다는 이유만으로 validation 결과를 보고 학습 가중치를 올리지 않는다. 첫 본실험은 uniform cell-type weighting으로 한다. 이후에도 Oligodendrocyte, Microglia-PVM, Astrocyte 같은 특정 세포형이 반복적으로 낮다면 train-only 난이도 기반의 capped weighting을 별도 실험으로 분리한다.

## 5. best epoch와 early stopping

test donor는 마지막 한 번의 최종 평가에만 사용한다. checkpoint 선택은 validation donor 16명만 사용한다.

매 epoch에는 싸게 계산되는 학습 로그를 저장하고, 매 2 epoch마다 고정 donor-balanced validation panel에서 다음을 계산한다.

1. primary: cell-type별 donor-centered module Spearman의 중앙값
2. fairness: 위 Spearman의 25백분위수
3. branch safety: cross-disjoint recovery
4. reconstruction safety: ordinal NLL과 one-bin accuracy
5. leakage safety: sex/assay balanced accuracy
6. geometry safety: effective dimension, R_GF degeneracy 진단용 요약

후보 checkpoint는 다음 안전조건을 먼저 통과해야 한다.

- branch cross-disjoint가 기준 PRISM보다 0.02 초과 하락하지 않음
- validation NLL이 기준 PRISM보다 2% 초과 악화하지 않음
- sex/assay leakage가 기준 checkpoint보다 통계적으로 악화하지 않음
- effective dimension이 붕괴하지 않음

안전조건을 통과한 checkpoint 중 primary module Spearman이 가장 높은 epoch를 biological best로 고른다. 차이가 validation bootstrap 95% 신뢰구간 안이면 더 이른 epoch를 선택한다. 자동 `loss/total` best와 biological best를 혼동하지 않고 둘 다 기록한다.

## 6. 실험 순서와 인과 분리

### 단계 0: 단위시험과 bf16 producer smoke

- projection target/gradient 시험
- CCC의 finite-gradient 및 zero-variance mask 시험
- rescue parameter gradient allow-list 시험
- 24세포형×4 draw schedule 시험
- posterior-mean deterministic rescue 재현성 시험

### 단계 1: whitening 단일변수 대조

구현 결정(2026-08-04): 이 단계는 별도의 장기 본학습 arm으로 실행하지 않고 해석적 동등성 및 bf16 회귀 smoke를 통과하는 gate로 대체한다. 정사영 뒤의 공분산이 $C=PCP$이고 $P$를 detach하며 분모를 $32^2$로 유지하면 $\mathcal L_I=\lVert C-I\rVert_F^2/32^2=\lVert C-P\rVert_F^2/32^2+4/1024=\mathcal L_P+0.00390625$이고 따라서 $\nabla_\theta\mathcal L_I=\nabla_\theta\mathcal L_P$다. 이 관계와 fp32 projector 성질을 smoke에서 확인하며, 이후 성능 변화는 whitening 단독 효과가 아니라 module-rescue v2 패키지 효과로만 해석한다.

기존 recovery 설정을 모두 고정하고 identity target만 projector target으로 바꾼다. 같은 분모를 쓰므로 큰 성능 차이가 없어야 한다. 목적은 수식과 로그의 정확성 확인이다.

### 단계 2A: rescue update 범위와 표본 안정화

projection-aware whitening을 고정하고 다음만 바꾼다.

- rescue step의 parameter allow-list
- 24세포형당 4개의 non-overlapping draw
- deterministic posterior-mean rescue
- branch:full을 25:75로 변경

Huber-only로 먼저 돌려 CCC의 효과와 섞이지 않게 한다.

### 단계 2B: endpoint-aligned CCC 추가

2B 설정은 train-only gradient calibration 기록의 SHA256과 고정 $\lambda_{CCC}$가 있어야만 생성·실행할 수 있다. 현재 생성된 본학습 설정은 인과 분리를 위한 2A Huber-only이며, 2B 본학습 설정은 아직 만들지 않았다.

2A를 초기값이나 teacher로 물려받지 않고 동일한 시작 checkpoint에서 독립적으로 학습한다. 유일한 차이는 calibrated CCC 항이다. 2A와 비교해 CCC가 donor-centered Spearman을 실제로 올리는지 판정한다.

### 단계 3: 최종 사후분석

- cell type별 module recovery
- cell type×disease gene/module recovery
- sex/assay/cell-type leakage
- common/personal/response branch ablation
- R_GF normal-versus-disease shape test
- geodesic-versus-Euclidean comparison
- Hodge decomposition
- generator functional Gram

이 단계까지 통과한 뒤에만 generator 수 선택을 시작한다.

## 7. generator 수 학습을 지금 하지 않는 이유

현재 module-rescue가 CT64보다 낮은 상태에서 generator를 줄이면 “생성자가 불필요해서 꺼진 것”과 “모듈을 아직 잘 못 배워서 gradient가 없었던 것”을 구분할 수 없다. 따라서 1·2번 실험에서는 414개를 전부 고정한다.

module-rescue v2가 확정된 뒤에는 다음 순서로 별도 진행한다.

1. 고정 checkpoint에서 functional Gram으로 선형 중복 후보를 선별한다.
2. deterministic hard top-$k$ 제거로 validation Pareto 곡선을 측정한다.
3. 평평한 구간의 1~2개 $k$를 gate 없이 고정해 처음부터 재학습한다.
4. 고정-$k$ 재학습이 all-414의 module, disease, leakage, R_GF, geodesic 안전조건을 만족할 때만 최종 개수로 채택한다.

soft L1/continuous gate나 확률적 Bernoulli gate를 이번 단계에서 쓰지 않는다. 겹치는 모듈 사이의 우회와 불안정한 개수 선택을 다시 만들 위험이 크기 때문이다.

## 8. 성공·실패 판정

module-rescue v2의 1차 성공은 validation donor에서 판단한다.

- donor-centered module recovery가 기존 PRISM recovery보다 명확히 상승
- CT64와의 격차가 감소
- cross-disjoint disease recovery의 하락이 0.02 이내
- NLL 악화가 2% 이내
- sex/assay leakage 악화 없음
- common/personal/response 중 필수 branch가 제거검정에서 유지
- R_GF, geodesic, Hodge의 기존 PRISM 장점이 재현

평균만 오르고 일부 세포형이 붕괴하면 성공으로 보지 않는다. 중앙값과 25백분위수를 함께 보고한다. CCC를 넣어도 Huber-only 2A보다 validation Spearman이 좋아지지 않으면 CCC 항은 폐기하고 2A를 채택한다.

## 9. 쉬운 비유

현재 PRISM은 24개 반 학생을 모두 가르치고 있다. 문제는 한 학생의 성적표를 만들 때 연습문제 두 개만 보고 판단해 점수가 흔들리고, 학습 문제는 “정답과 숫자가 얼마나 가까운가”만 묻는데 최종 시험은 “학생들의 순서를 제대로 맞혔는가”도 묻는다는 점이다.

v2에서는 한 번에 드는 책상 수는 그대로 두되, 서로 다른 문제지를 더 여러 번 풀게 한다. 숫자 오차를 줄이는 Huber 문제와 순서·변동을 맞히는 CCC 문제를 함께 준다. 보충수업인 rescue 시간에는 학교 전체 규칙이나 학생 신분표를 바꾸지 못하게 하고, 실제 답을 만드는 decoder 부분만 고치게 한다.

성별 방향 네 개를 지운 뒤의 whitening은 32개 의자를 모두 채우라고 요구하지 않고, 실제로 남은 28개 의자만 고르게 쓰라고 요구하는 수정이다. 다만 기존 손실과 gradient가 거의 같으므로 이것은 주로 채점표를 수학적으로 정확하게 만드는 일이고, 성적 향상의 주역은 module-rescue v2가 되어야 한다.

## 10. 2026-08-04 구현·검토 상태

- projection-aware whitening: CUDA bf16 autocast 차단 수정 및 producer 회귀시험 PASS
- Huber-only module-rescue v2A: 수식·학습통합 PASS, precision 전체가 아니라 명시적 출력 dictionary만 허용하도록 allow-list 축소
- CCC-enabled v2B: 목적함수 시험 PASS, train-only gradient calibration 전 실행 NO-GO
- 이번 수정 관련 회귀시험 80개 PASS, CUDA bf16 대상 시험 PASS
- 실제 sv6 6-GPU·sv7 4-GPU topology smoke: 사후분석 종료 후 실행 예정
- 본학습: 시작하지 않음

본학습 시작 전에는 sv6·sv7의 warm-start checkpoint SHA256과 `precision_head` 존재 여부, module-rescue stats artifact의 byte 또는 semantic SHA256, 빈 output directory를 확인한다. 실제 4-GPU·6-GPU NCCL smoke에서 전역 12 block/update, allow-list hash, gradient checksum도 기록한다. 이 provenance가 일치하지 않으면 두 서버 결과를 paired experiment로 해석하지 않는다.

Rescue allow-list는 decoder coefficient/generator/direct-state head와 `precision_head`의 명시적 `common_*`, `personal_basis`, `response_basis` 및 활성화된 경우의 interaction 출력 dictionary만 포함한다. `context_encoder`, `personal_posterior`, `normal_region_delta`, `age_*`, adversary는 rescue-only step에서 고정한다. `source_config.json`은 원본 JSON으로 보존하고, 실제 환경변수·GPU topology·라이브러리 버전·git 상태·핵심 소스 SHA256은 별도 `runtime_manifest.json`에 기록한다.
