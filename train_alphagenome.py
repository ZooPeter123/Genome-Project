import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm  # 진행률 표시줄 라이브러리

# --- 모델 정의 ---
class SimpleAlphaGenome(nn.Module):
    def __init__(self, input_length=1000, num_tracks=1):
        super().__init__()
        self.conv1 = nn.Conv1d(4, 32, kernel_size=7, padding=3)
        self.relu = nn.ReLU()
        self.conv2 = nn.Conv1d(32, 64, kernel_size=7, padding=3)
        self.final = nn.Conv1d(64, num_tracks, kernel_size=1)
        
    def forward(self, x):
        x = self.relu(self.conv1(x))
        x = self.relu(self.conv2(x))
        return self.final(x)

# --- 1. 가상의 유전체 데이터셋 생성 ---
class SyntheticDNADataset(Dataset):
    def __init__(self, num_samples=500, seq_len=1000): # 길이 1000bp로 수정됨
        self.num_samples = num_samples
        self.seq_len = seq_len
        self.motif = [2, 0, 3, 0] # G(2), A(0), T(3), A(0)
        
    def __len__(self):
        return self.num_samples
    
    def __getitem__(self, idx):
        # 1. 랜덤 DNA 생성
        seq_int = np.random.randint(0, 4, self.seq_len)
        
        # 2. Target Signal 생성 (기본 잡음)
        target = np.random.normal(0, 0.1, self.seq_len)
        
        # 3. 모티프 심기
        # 1000bp 길이 내에서 안전하게 모티프를 심을 공간 확보
        insert_loc = np.random.randint(100, self.seq_len - 100)
        seq_int[insert_loc:insert_loc+4] = self.motif
        target[insert_loc:insert_loc+4] += 5.0 # Signal Peak
        
        # 4. One-hot Encoding
        seq_onehot = np.zeros((4, self.seq_len), dtype=np.float32)
        for i, val in enumerate(seq_int):
            seq_onehot[val, i] = 1.0
            
        return torch.tensor(seq_onehot), torch.tensor(target, dtype=torch.float32).unsqueeze(0)

# --- 2. 학습 함수 (게이지바 추가됨) ---
def train_model(model, dataloader, epochs=5):
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    print(f"--- AlphaGenome 학습 시작 (총 {epochs} Epochs) ---")
    model.train()
    
    # 전체 Epoch 루프
    for epoch in range(epochs):
        # tqdm을 사용하여 진행률 게이지 생성
        with tqdm(dataloader, unit="batch") as tepoch:
            tepoch.set_description(f"Epoch {epoch+1}/{epochs}")
            
            epoch_loss = 0
            for seq, target in tepoch:
                optimizer.zero_grad()
                output = model(seq)
                loss = criterion(output, target)
                loss.backward()
                optimizer.step()
                
                epoch_loss += loss.item()
                # 게이지 바 옆에 현재 Loss 값을 실시간으로 표시
                tepoch.set_postfix(loss=loss.item())
                
    print("--- 학습 완료 ---")

# --- 3. 변이 효과 예측 (Peak 타격 로직 적용됨) ---
def simulate_variant_effect(model, original_seq_onehot):
    model.eval()
    with torch.no_grad():
        # A. 정상 서열 예측 (Wild Type)
        pred_wt = model(original_seq_onehot.unsqueeze(0)).squeeze().numpy()
        
        # B. 가장 신호가 높은 위치(Peak) 찾기 - '진짜 GATA' 위치 추적
        peak_idx = np.argmax(pred_wt) 
        print(f"\n[분석 결과]")
        print(f"- 가장 강력한 신호 위치 발견: {peak_idx}bp (Signal: {pred_wt[peak_idx]:.4f})")
        
        # C. Peak 위치 타격 (돌연변이 생성)
        mutant_seq = original_seq_onehot.clone()
        
        # Peak 위치의 염기를 강제로 변경 (0->1, 1->0 등으로 반전시켜 모티프 파괴)
        current_base_idx = torch.argmax(mutant_seq[:, peak_idx])
        new_base_idx = (current_base_idx + 1) % 4 
        
        mutant_seq[:, peak_idx] = 0
        mutant_seq[new_base_idx, peak_idx] = 1
        
        print(f"- 병원성 변이(Pathogenic Variant) 유도됨: {peak_idx}bp")
        
        # D. 변이 서열 예측
        pred_mut = model(mutant_seq.unsqueeze(0)).squeeze().numpy()
        
    return pred_wt, pred_mut

# --- 메인 실행부 ---
if __name__ == "__main__":
    # 데이터 준비 (1000bp로 넉넉하게 설정)
    dataset = SyntheticDNADataset(num_samples=1000, seq_len=1000)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
    
    # 모델 준비
    model = SimpleAlphaGenome(input_length=1000, num_tracks=1)
    
    # 학습 실행 (게이지가 차오르는 것을 확인하세요)
    train_model(model, dataloader, epochs=5)
    
    # 변이 효과 분석
    sample_seq, sample_target = dataset[0] 
    pred_wt, pred_mut = simulate_variant_effect(model, sample_seq)
    
    # 결과 시각화
    plt.figure(figsize=(10, 5))
    plt.plot(pred_wt, label='Wild Type (Original)', color='blue', alpha=0.7)
    plt.plot(pred_mut, label='Mutant (Pathogenic)', color='red', linestyle='--', alpha=0.7)
    plt.title("Effect of Pathogenic Mutation on Predicted Signal")
    plt.xlabel("Genomic Position (bp)")
    plt.ylabel("Predicted Signal")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()