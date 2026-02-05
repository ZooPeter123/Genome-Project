import torch
import torch.nn as nn
import torch.nn.functional as F

class ResidualConvBlock(nn.Module):
    """
    논문의 Conv Block을 단순화한 형태.
    DNA 서열의 국소적인 모티프(Motif)를 학습합니다.
    """
    def __init__(self, in_channels, out_channels, kernel_size=5):
        super().__init__()
        self.conv1 = nn.Conv1d(in_channels, out_channels, kernel_size, padding=kernel_size//2)
        self.bn1 = nn.BatchNorm1d(out_channels)
        self.conv2 = nn.Conv1d(out_channels, out_channels, kernel_size, padding=kernel_size//2)
        self.bn2 = nn.BatchNorm1d(out_channels)
        self.gelu = nn.GELU()
        
        # 차원이 다를 경우 잔차 연결(Skip Connection)을 위한 1x1 Conv
        self.resize = nn.Conv1d(in_channels, out_channels, 1) if in_channels != out_channels else nn.Identity()

    def forward(self, x):
        residual = self.resize(x)
        x = self.gelu(self.bn1(self.conv1(x)))
        x = self.bn2(self.conv2(x))
        return self.gelu(x + residual)

class MiniAlphaGenome(nn.Module):
    def __init__(self, input_length=16384, num_tracks=5):
        super().__init__()
        self.input_length = input_length
        
        # --- 1. Encoder (Downsampling) ---
        # 해상도를 줄여가며 특징을 추출 (마치 현미경 배율을 낮추며 큰 그림을 보는 과정)
        # Input: [Batch, 4, L] -> L은 서열 길이
        self.stem = nn.Conv1d(4, 64, kernel_size=7, padding=3)
        
        self.enc1 = ResidualConvBlock(64, 64)
        self.pool1 = nn.MaxPool1d(2) # L -> L/2
        
        self.enc2 = ResidualConvBlock(64, 128)
        self.pool2 = nn.MaxPool1d(2) # L/2 -> L/4
        
        self.enc3 = ResidualConvBlock(128, 256)
        self.pool3 = nn.MaxPool1d(2) # L/4 -> L/8
        
        # --- 2. Transformer Tower (Long-range Interaction) ---
        # 논문에서는 128bp 해상도에서 수행하지만, 여기서는 L/8 해상도에서 수행
        # Enhancer와 Promoter 사이의 먼 거리 상호작용을 계산 
        transformer_dim = 256
        encoder_layer = nn.TransformerEncoderLayer(d_model=transformer_dim, nhead=4, batch_first=True)
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=2)
        
        # --- 3. Decoder (Upsampling) ---
        # U-Net 구조: 인코더의 특징을 가져와(Concatenation) 복원
        self.up3 = nn.ConvTranspose1d(256, 128, kernel_size=2, stride=2) 
        self.dec3 = ResidualConvBlock(128 + 128, 128) # +128은 Skip Connection
        
        self.up2 = nn.ConvTranspose1d(128, 64, kernel_size=2, stride=2)
        self.dec2 = ResidualConvBlock(64 + 64, 64)
        
        self.up1 = nn.ConvTranspose1d(64, 64, kernel_size=2, stride=2)
        self.dec1 = ResidualConvBlock(64 + 64, 64)
        
        # --- 4. Prediction Heads ---
        # 특정 트랙(예: RNA-seq, DNase, ATAC 등)의 신호 예측 [cite: 3353]
        self.output_head = nn.Conv1d(64, num_tracks, kernel_size=1)

    def forward(self, x):
        # x shape: [Batch, 4, Length]
        
        # Encoding
        x = self.stem(x)         # [B, 64, L]
        e1 = self.enc1(x)        # [B, 64, L] -> Skip Connection용 저장
        x = self.pool1(e1)       # [B, 64, L/2]
        
        e2 = self.enc2(x)        # [B, 128, L/2] -> Skip Connection용 저장
        x = self.pool2(e2)       # [B, 128, L/4]
        
        e3 = self.enc3(x)        # [B, 256, L/4] -> Skip Connection용 저장
        x = self.pool3(e3)       # [B, 256, L/8]
        
        # Transformer Processing
        # Transformer는 (Batch, Seq, Feature) 형태를 원하므로 차원 변경
        x = x.permute(0, 2, 1)   # [B, L/8, 256]
        x = self.transformer(x)
        x = x.permute(0, 2, 1)   # [B, 256, L/8]
        
        # Decoding with Skip Connections (U-Net style) 
        x = self.up3(x)          # [B, 128, L/4]
        if x.size(2) != e3.size(2): x = F.interpolate(x, size=e3.size(2)) # 크기 보정
        x = torch.cat([x, e3], dim=1) # Skip Connection 결합
        x = self.dec3(x)
        
        x = self.up2(x)          # [B, 64, L/2]
        if x.size(2) != e2.size(2): x = F.interpolate(x, size=e2.size(2))
        x = torch.cat([x, e2], dim=1)
        x = self.dec2(x)
        
        x = self.up1(x)          # [B, 64, L]
        if x.size(2) != e1.size(2): x = F.interpolate(x, size=e1.size(2))
        x = torch.cat([x, e1], dim=1)
        x = self.dec1(x)
        
        # Final Prediction
        out = self.output_head(x) # [B, Num_Tracks, L]
        return out

# --- 실행 및 테스트 코드 ---
if __name__ == "__main__":
    # GPU 사용 가능 여부 확인
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # 1. 모델 인스턴스화
    # 입력 길이: 16kb, 예측 트랙 수: 2개 (예: RNA-seq coverage, ATAC-seq)
    model = MiniAlphaGenome(input_length=16384, num_tracks=2).to(device)
    
    # 2. 가상의 DNA 데이터 생성 (One-hot encoding)
    # Batch Size: 4, Channels: 4 (ACGT), Length: 16384
    dummy_input = torch.randn(4, 4, 16384).to(device)
    
    # 3. 모델 예측 (Forward Pass)
    output = model(dummy_input)
    
    print(f"Input shape: {dummy_input.shape}")   # [4, 4, 16384]
    print(f"Output shape: {output.shape}")       # [4, 2, 16384] -> 염기쌍 별 예측값
    print("구현 성공! AlphaGenome의 축소판이 정상적으로 작동합니다.")