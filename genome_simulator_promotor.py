import random

# ==========================================
# 1. 기초 설정 (Codon Table & Tools)
# ==========================================
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def generate_random_seq(length):
    return ''.join(random.choices(['A', 'T', 'G', 'C'], k=length))

# ==========================================
# 2. 진핵생물 유전자 생성기 (Synthetic Biology)
# ==========================================
# 1. 유전자 생성기 수정 (더 강력한 신호 심기)
def create_eukaryotic_gene():
    promoter = "TATAAA" 
    exon1 = "ATGAAGCTAGTT"  # M-K-L-V (끝이 GTT로 끝남)
    
    # [수정] 인트론 시작을 'GT'가 아니라 더 명확한 'GTAAGT' (Consensus)로 변경
    # 이렇게 하면 Exon의 GTT와 헷갈리지 않습니다.
    intron = "GTAAGT" + generate_random_seq(20) + "AG" 
    
    exon2 = "TGGCCATAA" 
    enhancer = "GGCCGG"
    
    full_dna = generate_random_seq(10) + enhancer + \
               generate_random_seq(15) + promoter + \
               generate_random_seq(5) + \
               exon1 + intron + exon2 + \
               generate_random_seq(20)
               
    return full_dna, (exon1 + exon2)

# ==========================================
# 3. 핵심 엔진: Regulation & Processing
# ==========================================

def step1_transcription_regulation(dna_seq):
    """
    [Promoter 탐색]
    TATA Box('TATAAA')가 있어야만 전사를 시작합니다.
    """
    print(f"\n[Step 1] Finding Promoter (TATA Box)...")
    
    if "TATAAA" in dna_seq:
        start_idx = dna_seq.find("TATAAA") + 6 # TATA box 뒤부터 전사 시작 가정
        print(f"  ✅ Promoter found at index {dna_seq.find('TATAAA')}! Initiating Transcription.")
        
        # DNA -> Pre-mRNA (T -> U)
        pre_mrna = dna_seq[start_idx:].replace('T', 'U')
        return pre_mrna
    else:
        print("  ❌ No Promoter found. Gene is silenced.")
        return None

# 2. 스플라이싱 함수 수정 (강력한 신호 찾기)
def step2_splicing(pre_mrna):
    print(f"\n[Step 2] Splicing (Removing Introns)...")
    
    # [수정] 찾는 신호를 'GU' -> 'GUAAGU' (RNA 기준)로 변경
    # 훨씬 엄격한 기준을 적용하여 Exon 내부의 GU를 무시함
    DONOR_SIGNAL = "GUAAGU" 
    
    intron_start = pre_mrna.find(DONOR_SIGNAL) # GUAAGU 찾기
    
    # 이하는 동일한 논리 (단, start 위치 보정 필요 없음, 자르는 위치는 그대로)
    if intron_start != -1:
        MIN_INTRON_LEN = 20
        search_start_pos = intron_start + MIN_INTRON_LEN
        intron_end = pre_mrna.find("AG", search_start_pos)
        
        if intron_end != -1:
            intron_seq = pre_mrna[intron_start : intron_end+2]
            mature_mrna = pre_mrna[:intron_start] + pre_mrna[intron_end+2:]
            
            print(f"  ✂️  Splicing Success! (Consensus site found)")
            print(f"      - Intron removed: {intron_seq[:10]}... (Length: {len(intron_seq)})")
            return mature_mrna

    print("  ⚠️ No valid splicing signal found.")
    return pre_mrna

def step3_translation(mature_mrna):
    """
    [Translation]
    Mature mRNA -> Polypeptide
    """
    print(f"\n[Step 3] Translation (Ribosome Activity)...")
    
    # Start Codon (AUG) 찾기
    start_idx = mature_mrna.find("AUG")
    if start_idx == -1:
        return ""
        
    protein_seq = []
    # AUG부터 3개씩 읽기
    for i in range(start_idx, len(mature_mrna), 3):
        codon = mature_mrna[i:i+3]
        if len(codon) < 3: break
        
        aa = CODON_TABLE.get(codon.replace('U','T'), '?') # Table은 DNA 기준이라 치환
        
        if aa == '*': # Stop Codon
            break
        protein_seq.append(aa)
        
    polypeptide = "".join(protein_seq)
    print(f"  ✅ Polypeptide synthesized: {polypeptide}")
    return polypeptide

def step4_post_translational_modification(polypeptide):
    """
    [PTM: 인산화 (Phosphorylation)]
    단백질의 특정 아미노산(Serine 'S', Tyrosine 'Y', Threonine 'T')에 
    인산기(Phosphate)를 붙여 활성화시킵니다.
    """
    print(f"\n[Step 4] Post-translational Modification (PTM)...")
    
    modified_protein = []
    modification_count = 0
    
    for aa in polypeptide:
        # Serine(S)이나 Tyrosine(Y)이 오면 인산화된다고 가정
        if aa in ['S', 'Y', 'T']:
            modified_protein.append(f"{aa}[Phos]")
            modification_count += 1
        else:
            modified_protein.append(aa)
            
    final_structure = "-".join(modified_protein)
    
    if modification_count > 0:
        print(f"  ⚡ Phosphorylation applied to {modification_count} residues.")
    else:
        print(f"  ⚪ No target residues for PTM found.")
        
    return final_structure

# ==========================================
# 4. 실행 (Main Pipeline)
# ==========================================
if __name__ == "__main__":
    print("--- Eukaryotic Genome Simulator v2.0 ---")
    
    # 1. 합성 유전자 생성 (DNA)
    my_gene, answer_key = create_eukaryotic_gene()
    print(f"🧬 Genomic DNA generated (Length: {len(my_gene)} bp)")
    print(f"   Seq: {my_gene}")
    
    # 2. 전사 (with Promoter Check)
    pre_mrna = step1_transcription_regulation(my_gene)
    
    if pre_mrna:
        # 3. 스플라이싱 (Intron Removal)
        mature_mrna = step2_splicing(pre_mrna)
        
        # 4. 번역
        polypeptide = step3_translation(mature_mrna)
        
        # 5. 번역 후 변형 (PTM)
        final_protein = step4_post_translational_modification(polypeptide)
        
        print(f"\n🎉 Final Functional Protein: {final_protein}")
        
        # 검증: 우리가 설계한 Exon이 제대로 붙었나 확인
        # (단, 앞뒤 Junk 서열 때문에 완전히 같진 않고 포함 관계 확인)
        print(f"\n[Verification]")
        print(f"Expected Coding Seq (DNA): {answer_key}")
        # mRNA에서 U를 T로 바꿔서 비교
        is_correct = answer_key in mature_mrna.replace('U', 'T')
        print(f"Splicing Logic Correct? : {is_correct}")