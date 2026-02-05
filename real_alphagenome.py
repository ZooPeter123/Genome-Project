import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
from tqdm import tqdm

# --- 1. Numpy 데이터를 위한 데이터셋 클래스 (라이브러리 의존성 제거됨) ---
class NumpyGenomicDataset(Dataset):
    def __init__(self, seq_path, signal_path, length=8192, num_samples=2000):
        """
        이제 복잡한 파싱 없이 Numpy 배열을 직접 로드합니다.
        속도가 훨씬 빠르고, 윈도우/맥/리눅스 어디서든 돌아갑니다.
        """
        print(f"📂 데이터 로딩 중... ({seq_path})")
        # 메모리 맵(mmap) 모드를 사용하여 RAM을 아낍니다.
        self.seq_data = np.load(seq_path, mmap_mode='r') 
        self.signal_data = np.load(signal_path, mmap_mode='r')
        
        self.chrom_len = len(self.seq_data)
        self.length = length
        self.num_samples = num_samples
        print("✅ 데이터 로드 완료!")

    def __len__(self):
        return self.num_samples

    def __getitem__(self, idx):
        # 1. 랜덤 위치 선택
        start = np.random.randint(0, self.chrom_len - self.length)
        end = start + self.length
        
        # 2. Sequence 가져오기 (이미 숫자로 변환되어 있음)
        seq_int = self.seq_data[start:end]
        # One-hot Encoding (즉석 변환)
        # 0->A, 1->C, 2->G, 3->T, 4->N
        seq_tensor = torch.zeros(4, self.length, dtype=torch.float32)
        
        # 벡터화된 연산으로 고속 처리
        seq_tensor[0, seq_int == 0] = 1.0 # A
        seq_tensor[1, seq_int == 1] = 1.0 # C
        seq_tensor[2, seq_int == 2] = 1.0 # G
        seq_tensor[3, seq_int == 3] = 1.0 # T
        
        # 3. Signal 가져오기
        signal_val = self.signal_data[start:end]
        signal_tensor = torch.tensor(signal_val, dtype=torch.float32).unsqueeze(0)
        
        return seq_tensor, signal_tensor

# --- 2. 모델 아키텍처 (이전과 동일) ---
class ResBlock(nn.Module):
    def __init__(self, channels, kernel_size=5, dilation=1):
        super().__init__()
        self.conv1 = nn.Conv1d(channels, channels, kernel_size, padding='same', dilation=dilation)
        self.bn1 = nn.BatchNorm1d(channels)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size, padding='same', dilation=dilation)
        self.bn2 = nn.BatchNorm1d(channels)
        self.gelu = nn.GELU()

    def forward(self, x):
        residual = x
        out = self.gelu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        return self.gelu(out + residual)

class AlphaGenomeMedium(nn.Module):
    def __init__(self, input_len=32768, num_tracks=1):
        super().__init__()
        self.stem = nn.Sequential(nn.Conv1d(4, 128, kernel_size=15, padding=7), nn.GELU(), nn.BatchNorm1d(128))
        self.res_blocks = nn.Sequential(
            ResBlock(128, dilation=1), nn.MaxPool1d(2),
            ResBlock(128, dilation=2), nn.MaxPool1d(2),
            ResBlock(128, dilation=4), nn.MaxPool1d(2),
        )
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d_model=128, nhead=4, batch_first=True, dim_feedforward=512), num_layers=2
        )
        self.decoder = nn.Sequential(
            nn.ConvTranspose1d(128, 64, kernel_size=4, stride=4), nn.GELU(),
            nn.ConvTranspose1d(64, 32, kernel_size=2, stride=2), nn.GELU(),
            nn.Conv1d(32, num_tracks, kernel_size=1)
        )

    def forward(self, x):
        x = self.stem(x)
        x = self.res_blocks(x)
        x = x.permute(0, 2, 1) 
        x = self.transformer(x)
        x = x.permute(0, 2, 1)
        return self.decoder(x)

# --- 3. 실행 함수 ---
def train_numpy_model():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"🚀 장치 설정: {device}")

    # Colab에서 다운로드 받은 파일 경로
    SEQ_PATH = "train_seq.npy"
    SIGNAL_PATH = "train_signal.npy"
    
    try:
        # 데이터셋 초기화 (실제 데이터)
        ds = NumpyGenomicDataset(SEQ_PATH, SIGNAL_PATH, num_samples=2000)
        loader = DataLoader(ds, batch_size=4, shuffle=True)
        
        model = AlphaGenomeMedium().to(device)
        optimizer = optim.AdamW(model.parameters(), lr=1e-3)
        criterion = nn.MSELoss()
        
        print("--- 학습 시작 ---")
        model.train()
        for epoch in range(5):
            with tqdm(loader, unit="batch", desc=f"Ep {epoch+1}") as tepoch:
                for seq, target in tepoch:
                    seq, target = seq.to(device), target.to(device)
                    
                    optimizer.zero_grad()
                    output = model(seq)
                    loss = criterion(output, target)
                    loss.backward()
                    optimizer.step()
                    
                    tepoch.set_postfix(loss=loss.item())
        print("🎉 학습 완료! 저장된 모델: best_alphagenome.pth")
        torch.save(model.state_dict(), "best_alphagenome.pth")
        
    except FileNotFoundError:
        print(f"❌ 오류: '{SEQ_PATH}' 파일을 찾을 수 없습니다.")
        print("   1단계(Colab)를 실행하여 .npy 파일을 먼저 생성하고 다운로드해주세요.")

if __name__ == "__main__":
    train_numpy_model()