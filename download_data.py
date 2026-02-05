import os
import requests
import gzip
import shutil
from tqdm import tqdm

def download_file(url, save_name):
    """URL에서 파일을 다운로드하며 진행률(Bar)을 표시합니다."""
    if os.path.exists(save_name):
        print(f"✅ 이미 파일이 존재합니다: {save_name}")
        return

    print(f"⬇️ 다운로드 시작: {save_name} ...")
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    
    with open(save_name, 'wb') as file, tqdm(
        desc=save_name,
        total=total_size,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for data in response.iter_content(chunk_size=1024):
            size = file.write(data)
            bar.update(size)
    print(f"🎉 다운로드 완료: {save_name}")

def decompress_gz(gz_path, out_path):
    """gz 압축을 해제합니다."""
    if os.path.exists(out_path):
        return
    
    print(f"📦 압축 해제 중: {gz_path} -> {out_path}")
    with gzip.open(gz_path, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print("✨ 압축 해제 완료!")

def main():
    # 1. 인간 참조 유전체 (UCSC hg38 - 22번 염색체)
    # 전체 hg38은 너무 크므로(3GB+), 테스트용으로 22번(약 10MB)만 받습니다.
    fasta_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
    fasta_gz = "chr22.fa.gz"
    fasta_final = "chr22.fa"
    
    download_file(fasta_url, fasta_gz)
    decompress_gz(fasta_gz, fasta_final)

    # 2. 실제 실험 데이터 (ENCODE GM12878 DNase-seq Signal)
    # ENCODE 공식 포털의 실제 파일 URL입니다. (약 400~500MB)
    # ENCFF901GZH: GRCh38 기반의 Read-depth normalized signal
    bigwig_url = "https://www.encodeproject.org/files/ENCFF901GZH/@@download/ENCFF901GZH.bigWig"
    bigwig_final = "DNase.bigWig"
    
    download_file(bigwig_url, bigwig_final)

    print("\n✅ 모든 데이터 준비 완료! 이제 'real_alphagenome.py'를 실행하세요.")
    print(f"   - FASTA: {os.path.abspath(fasta_final)}")
    print(f"   - BigWig: {os.path.abspath(bigwig_final)}")

if __name__ == "__main__":
    main()