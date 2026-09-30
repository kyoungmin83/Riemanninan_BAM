# Compartment E52 아키텍처·학습 코드 검토 묶음

작성일: 2026-09-23. 논문 후보인 기존 compartment E52를 검토하기 위한 복사본이다. 최신 RLT/Dream 학습 실행본이 아니다. 원래 학습 파일은 변경하지 않았다.

## 먼저 읽을 순서

1. `model_card/MODEL_CARD_KO.md`: 모델과 평가 설명.
2. `model_card/CHECKPOINT_MANIFEST.json`: E52 checkpoint 식별자·당시 실행 경로·기록된 해시.
3. `model_card/cfg/source_config.json`, `resolved_config.json`: E52의 설정 기준. 루트 `configs/`의 다른 실험 설정을 대신 선택하지 말 것.
4. `docs/prism_integrated_module_local_residual_v1_design_20260825.md`: 모듈별 개인 경로.
5. `docs/prism_integrated_compartmental_nonlinearity_techzero_v2_design_20260827.md`: 구획 비선형 설계.
6. `docs/prism_compartmental_techzero_v2_implementation_report_20260827.md`: 구현 기록.
7. `review_notes/`: 이후 발견된 학습 문제와 사후분석. 설계 의도와 실제 검증 결과는 구분할 것.

## 코드 지도

- `src/kmlee_bam/model/`: encoder, latent nuisance projection, decoder, system 등 모델 코드.
- `src/kmlee_bam/training/`: 학습 루프, loss 결합, 생성자 선택, module rescue.
- `src/kmlee_bam/objectives/`: 목적함수.
- `scripts/`, `configs/`, `tests/`: 당시 실행 디렉터리의 보조 코드·설정·검사.
- `preprocessing_reference/`: 현재 PRISM 작업 폴더에서 가져온 전처리 보조 자료. 당시 E52와 파일별 동일성이 검증된 동결본이라는 뜻은 아니다.

쉽게 말하면 공통 병리 경로는 여러 사람에게 공통된 변화, 개인 기본 경로는 사람마다 원래 다른 특성, 개인 반응 경로는 같은 병리에 대한 서로 다른 반응을 표현한다. 모듈은 유전자 묶음을 읽는 좌표이고 생성자 gate는 모델 용량 선택 장치다. 생성자 제거를 해당 모듈 정보 전체의 삭제로 해석하면 안 된다. 구획이라는 이름도 실제 수상돌기의 인과적 측정을 뜻하지 않는다.

## 출처와 재현 범위

`src/`, `scripts/`, `configs/`, `docs/`, `tests/`는 SV6의 아래 실행 디렉터리에서 복사했다.

`/home/kmlee/project_sv6/kmlee_bam_integrated_compartmental_pathrank_20260828/`

이는 2026-09-23에 해당 위치에서 읽은 복사본이다. 디렉터리 이름만으로 모든 파일이 학습 당시와 bitwise 동일하다고 보증하지 않는다. `FILE_MANIFEST.json`은 이번 복사본의 파일 해시와 양 서버 간 동일성 확인용이며, 역사적 소스 해시 검증과 다르다.

E52 원본 run (SV6):

`/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_freegen_recovery_e12_20260829/`

주 checkpoint: `checkpoint_epoch_052.pt`.

환자 원시 데이터, 대형 Zarr/NPZ 입력, checkpoint, 학습 로그 전체는 이 묶음에 복사하지 않았다. 설정에 기록된 외부 절대경로는 그대로 보존했다. E12 recovery checkpoint 및 module graph/reliability/cap/rescue 통계도 정확한 재학습에 필요하다. 따라서 이 묶음은 코드 검토용이며, 독립 실행 가능한 완전한 재현 패키지라고 주장하지 않는다.

현재 PRISM 루트의 README/SOURCE_PROVENANCE/REQUIRED_ARTIFACTS는 과거 rank2 E20에 관한 내용이므로 E52 안내로 혼용하지 않는다. 보조 전처리 자료의 provenance 역시 그 적용 범위를 확인해야 한다.

이 폴더의 어떤 실행 스크립트도 이번 복사 작업에서 실행하지 않았다. 기존 SV6/SV7 학습 및 원본은 변경하지 않았다.
