import random

# ==========================================
# 1. 환경 설정: 코돈표 (The Genetic Code)
# ==========================================
# 3개의 염기(Key)가 하나의 아미노산(Value)에 대응되는 딕셔너리 구조입니다.
# *: 종결 코돈 (Stop Codon)

CODON_TABLE = {
    # Phenylalanine
    'TTT': 'F', 'TTC': 'F',
    # Leucine
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    # Isoleucine
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    # Methionine (Start Codon)
    'ATG': 'M',
    # Valine
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    # Serine
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    # Proline
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    # Threonine
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    # Alanine
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    # Tyrosine
    'TAT': 'Y', 'TAC': 'Y',
    # Stop Codons
    'TAA': '*', 'TAG': '*',
    # Histidine
    'CAT': 'H', 'CAC': 'H',
    # Glutamine
    'CAA': 'Q', 'CAG': 'Q',
    # Asparagine
    'AAT': 'N', 'AAC': 'N',
    # Lysine
    'AAA': 'K', 'AAG': 'K',
    # Aspartic Acid
    'GAT': 'D', 'GAC': 'D',
    # Glutamic Acid
    'GAA': 'E', 'GAG': 'E',
    # Cysteine
    'TGT': 'C', 'TGC': 'C',
    # Tryptophan
    'TGG': 'W',
    # Arginine
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    # Glycine
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    # Stop Codon (Opale)
    'TGA': '*'
}

# ==========================================
# 2. 유전체 생성 엔진 (Genome Generator)
# ==========================================

def generate_dna(length=100, gc_content=0.5):
    """
    지정된 길이와 GC 함량 비율에 따라 무작위 DNA 서열을 생성합니다.
    
    Args:
        length (int): 생성할 DNA 서열의 길이 (bp)
        gc_content (float): 전체 서열 중 G와 C가 차지하는 비율 (0.0 ~ 1.0)
                            인간의 평균 GC content는 약 41% (0.41) 입니다.
    """
    
    # GC 함량에 따른 각 염기의 확률 가중치 계산
    # G와 C는 반반씩 나눠가짐 -> gc_content / 2
    # A와 T는 나머지(1 - gc_content)를 반반씩 나눠가짐
    p_gc = gc_content / 2
    p_at = (1 - gc_content) / 2
    
    bases = ['A', 'T', 'G', 'C']
    weights = [p_at, p_at, p_gc, p_gc] # 확률: [A, T, G, C] 순서
    
    # random.choices를 사용하여 가중치에 기반한 무작위 추출
    # k는 추출할 횟수(=DNA 길이)
    dna_list = random.choices(bases, weights=weights, k=length)
    
    return ''.join(dna_list)

# ==========================================
# 3. 실행 및 테스트 (Main Execution)
# ==========================================
if __name__ == "__main__":
    # 가상의 인간 유전체 조각 생성 (길이: 300bp, GC 함량: 41%)
    virtual_genome = generate_dna(length=300, gc_content=0.41)
    
    print(f"--- Virtual Human Genome Project (Phase 1) ---")
    print(f"Total Length: {len(virtual_genome)} bp")
    print(f"Sequence Preview (First 50bp): {virtual_genome[:50]}...")
    
    # 실제 GC Content 검증
    actual_gc = (virtual_genome.count('G') + virtual_genome.count('C')) / len(virtual_genome)
    print(f"Target GC: 41% | Actual GC: {actual_gc*100:.2f}%")

# (앞서 작성한 CODON_TABLE과 generate_dna 함수는 그대로 유지합니다)

# ==========================================
# 4. 전사 (Transcription): DNA -> mRNA
# ==========================================
def transcribe(dna_seq):
    """
    DNA 서열을 mRNA 서열로 변환합니다.
    생물학적 로직:
    1. Coding strand(5'->3')를 기준으로 한다고 가정.
    2. Thymine(T)을 Uracil(U)로 치환.
    """
    # 파이썬의 문자열 치환 함수인 replace를 사용하면 매우 간단합니다.
    return dna_seq.replace('T', 'U')

# ==========================================
# 5. 번역 (Translation): mRNA -> Protein
# ==========================================
def translate(rna_seq):
    """
    mRNA 서열을 3개씩(Triplet) 읽어 단백질(아미노산 서열)로 변환합니다.
    """
    protein_seq = []
    
    # 0부터 서열 길이까지 3칸씩 건너뛰며 반복 (Reading Frame 1 고정)
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i+3]
        
        # 마지막에 염기가 1~2개 남으면 코돈이 안 되므로 버림 (Drop)
        if len(codon) < 3:
            break
            
        # 딕셔너리에서 아미노산 찾기. 없으면 'X' (Unknown) 반환
        amino_acid = CODON_TABLE.get(codon, 'X')
        protein_seq.append(amino_acid)
        
    return "".join(protein_seq)

# ==========================================
# 6. 통합 테스트 (Integration Test)
# ==========================================
if __name__ == "__main__":
    # 1. DNA 생성 (짧게 30bp만 생성해서 눈으로 확인)
    dna = generate_dna(length=30, gc_content=0.5)
    
    # 2. 전사
    mrna = transcribe(dna)
    
    # 3. 번역
    protein = translate(mrna)
    
    print(f"--- Central Dogma Simulation ---")
    print(f"1. DNA (Coding Strand): {dna}")
    print(f"2. mRNA (T -> U)      : {mrna}")
    
    # 시각적으로 코돈과 아미노산을 매칭시켜 보여주기 위한 포맷팅
    # (아미노산 문자를 가운데 정렬하여 코돈 아래에 위치시킴)
    formatted_protein = "  ".join([aa for aa in protein])
    print(f"3. Protein            :  {formatted_protein}")
    
    print(f"\n[Result Analysis]")
    print(f"- DNA Length: {len(dna)} bp")
    print(f"- Protein Length: {len(protein)} aa")

# ==========================================
# 7. ORF 탐색기 (The Gene Hunter)
# ==========================================
def find_orfs_in_mrna(mrna_seq, min_len=5):
    """
    mRNA 서열 내에서 '개시 코돈(AUG)'부터 '종결 코돈'까지의 
    유효한 단백질 서열(ORF)을 모두 찾아냅니다.
    
    Args:
        mrna_seq (str): 입력 mRNA 서열
        min_len (int): 단백질로 인정할 최소 아미노산 길이 (노이즈 필터링)
        
    Returns:
        list: 발견된 단백질 서열들의 리스트
    """
    found_proteins = []
    n = len(mrna_seq)
    
    # 3가지 Reading Frame을 모두 스캔 (0, +1, +2)
    for frame in range(3):
        # frame 위치부터 시작해서 3칸씩 점프하며 스캔
        i = frame
        while i < n - 2:
            codon = mrna_seq[i:i+3]
            
            # 1. Start Codon 감지 (AUG)
            if codon == 'AUG':
                # 단백질 합성 시작
                temp_protein = []
                
                # Start 코돈 이후로 계속 번역 진행 (inner loop)
                for j in range(i, n - 2, 3):
                    current_codon = mrna_seq[j:j+3]
                    
                    # Stop Codon 감지 (UAA, UAG, UGA)
                    if current_codon in ['UAA', 'UAG', 'UGA']:
                        # 유효한 길이인지 확인 (너무 짧으면 노이즈로 간주)
                        if len(temp_protein) >= min_len:
                            found_proteins.append("".join(temp_protein))
                        
                        # 중요: 유전자를 찾았으면, 리보솜은 떨어져 나감.
                        # 메인 루프(i)를 현재 위치(j)로 점프시켜 중복 탐색 방지
                        i = j 
                        break
                    
                    # 일반 아미노산 번역
                    # (기존 CODON_TABLE이 DNA 기준(T)이므로 U를 T로 바꿔서 조회)
                    codon_dna = current_codon.replace('U', 'T')
                    aa = CODON_TABLE.get(codon_dna, '?')
                    temp_protein.append(aa)
                else:
                    # for문이 break 없이 끝남 = Stop codon 없이 서열이 끝남
                    # (완성되지 않은 단백질이므로 버림)
                    pass
                    
            # 다음 코돈으로 이동
            i += 3
            
    return found_proteins

# ==========================================
# 8. 대규모 시뮬레이션 (Execution)
# ==========================================
if __name__ == "__main__":
    # 실험: 확률적으로 유전자를 발견하기 위해 좀 더 긴 DNA 생성 (3,000 bp)
    print(f"\n--- Phase 3: Searching for Genes in 3,000bp ---")
    
    long_dna = generate_dna(length=3000, gc_content=0.45)
    long_mrna = transcribe(long_dna)
    
    # ORF 탐색 (최소 아미노산 10개 이상인 것만 유전자로 인정)
    discovered_genes = find_orfs_in_mrna(long_mrna, min_len=10)
    
    print(f"Total Sequence Length: {len(long_dna)} bp")
    print(f"Discovered ORFs (Genes): {len(discovered_genes)} found\n")
    
    # 발견된 유전자 목록 출력 (상위 5개만)
    for idx, protein in enumerate(discovered_genes[:5]):
        print(f"Gene {idx+1}: (Length {len(protein)}aa)")
        print(f"   Seq: M{protein}...") # 모든 단백질은 M(메티오닌)으로 시작
        print("-" * 30)

import textwrap # 긴 서열을 보기 좋게 줄바꿈하기 위한 도구

# ... (이전 단계의 CODON_TABLE, generate_dna, transcribe, translate, find_orfs_in_mrna 함수들은 그대로 유지) ...

# ==========================================
# 9. 파일 저장 (Data Archiving)
# ==========================================
def save_to_fasta(data_dict, filename, description=""):
    """
    딕셔너리 형태의 서열 데이터를 생물정보학 표준인 FASTA 포맷으로 저장합니다.
    
    Format:
    >Unique_ID Description
    SEQUENCE_DATA...
    """
    with open(filename, "w", encoding="utf-8") as f:
        for header, sequence in data_dict.items():
            # 헤더 작성 (> 기호로 시작)
            f.write(f">{header} {description}\n")
            
            # 서열 작성 (가독성을 위해 80글자마다 줄바꿈)
            wrapped_seq = textwrap.fill(sequence, width=80)
            f.write(wrapped_seq + "\n")
            
    print(f"✅ File saved successfully: {filename}")

# ==========================================
# 10. 최종 통합 시뮬레이션 (Final Execution)
# ==========================================
if __name__ == "__main__":
    print(f"--- Virtual Genome Project: Final Build ---\n")
    
    # 1. Genome 생성 (10,000 bp - 작은 바이러스 수준의 크기)
    print("1. Synthesizing Virtual Genome...", end="")
    genome_dna = generate_dna(length=10000, gc_content=0.45)
    print(" Done.")
    
    # 2. 전사 (Transcription)
    print("2. Transcribing to mRNA...", end="")
    genome_mrna = transcribe(genome_dna)
    print(" Done.")
    
    # 3. 유전자 발굴 (ORF Finding)
    print("3. Scanning for Genes (ORF)...", end="")
    # 최소 30 아미노산(약 90bp) 이상인 것만 유의미한 단백질로 인정
    found_proteins = find_orfs_in_mrna(genome_mrna, min_len=30)
    print(f" {len(found_proteins)} Potential Genes Found.")
    
    # 4. 데이터 구조화 (List -> Dictionary)
    # 저장하기 좋게 이름을 붙여줍니다 (예: Gene_001, Gene_002...)
    protein_db = {}
    for idx, seq in enumerate(found_proteins):
        gene_id = f"Virtual_Gene_{idx+1:03d}"
        protein_db[gene_id] = seq
        
    # 5. 파일로 저장 (Export)
    print("\n4. Exporting Data to FASTA files...")
    
    # 5-1. 전체 유전체 서열 저장
    save_to_fasta(
        {"Chromosome_1": genome_dna}, 
        "virtual_human_genome.fasta", 
        description="| Homo virtualis | Build 1.0"
    )
    
    # 5-2. 발견된 단백질 서열 저장
    save_to_fasta(
        protein_db, 
        "virtual_proteome.fasta", 
        description="| Predicted Protein | ORF Discovery"
    )
    
    print("\n[Simulation Complete]")
    print(f"📂 결과 파일이 생성되었습니다: 'virtual_human_genome.fasta', 'virtual_proteome.fasta'")
    if len(found_proteins) > 0:
        print(f"🎉 축하합니다! {len(found_proteins)}개의 새로운 가상 단백질을 발견했습니다.")