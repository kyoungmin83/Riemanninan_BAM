import scanpy as sc
import torch 
import os
import sys

data_path = os.path.expanduser("~/project/riemann_bam/data/tabula_blood_normal/TS_Blood.h5ad")

print(f"--- 환경 및 데이터 점검 시작 ---")
print (data_path)

# [GPU 및 CUDA 점검]
print("\n[1] GPU 하드웨어 및 CUDA 버전 확인")
print(f"▶ PyTorch 버전: {torch.__version__}")
print(f"▶ CUDA 사용 가능 여부: {torch.cuda.is_available()}")

if torch.cuda.is_available():
    print(f"▶ 현재 CUDA 버전 (PyTorch 빌드): {torch.version.cuda}")
    print(f"▶ 감지된 GPU 수: {torch.cuda.device_count()}")
    print(f"▶ 현재 GPU 장치: {torch.cuda.get_device_name(0)}")
    
    # RTX A6000은 48GB여야 함
    vram_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
    print(f"▶ VRAM 용량: {vram_gb:.2f} GB")
    
    if "A6000" in torch.cuda.get_device_name(0) and vram_gb > 40:
        print(">> [성공] RTX A6000 (48GB)가 정상적으로 인식되었습니다.")
    else:
        print(">> [주의] 예상된 GPU(A6000)가 아니거나 메모리가 다릅니다.")
else:
    print(">> [경고] GPU가 인식되지 않습니다! 'nvidia-smi'를 확인하세요.")

# [데이터 로딩]
print(f"\n[2] 데이터 로딩 테스트: {data_path}")

if not os.path.exists(data_path):
    print(f">> [에러] 파일이 존재하지 않습니다: {data_path}")
    print(">> 마운트 경로가 정확한지, 파일명이 맞는지 확인해주세요.")
else:
    try:
        adata = sc.read_h5ad(data_path)
        print(f"▶ 성공적으로 로드됨: {adata.n_obs} cells x {adata.n_vars} genes")
        
        # [데이터 상태 확인]
        print("\n[3] 데이터 값(Count vs Normalized) 확인")
        # 희소 행렬(Sparse) 처리
        if hasattr(adata.X, "toarray"):
            sample_data = adata.X[:5, :5].toarray()
            is_sparse = True
        else:
            sample_data = adata.X[:5, :5]
            is_sparse = False
            
        print(f"▶ 데이터 포맷: {'Sparse Matrix (희소행렬)' if is_sparse else 'Dense Matrix'}")
        print("▶ 샘플 데이터 (5x5):")
        print(sample_data)
        
        # 정수인지 소수인지 판별
        is_integer = bool((sample_data % 1 == 0).all())
        if is_integer:
            print(">> [분석] 데이터가 '정수(Raw Count)'입니다. (Negative Binomial Loss 추천)")
        else:
            print(">> [분석] 데이터가 '실수(Normalized/Log-transformed)'입니다. (Gaussian/MSE Loss 추천)")

    except Exception as e:
        print(f">> [에러] 데이터 로드 중 문제 발생: {e}")

print (adata.obs.head())