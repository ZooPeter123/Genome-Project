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
# 2. 유전자 생성기 (돌연변이 모드 추가)
# ==========================================
def create_eukaryotic_gene(apply_mutation=False):
    """
    Args:
        apply_mutation (bool): True일 경우 병적인 돌연변이를 심습니다.
    """
    promoter = "TATAAA" 
    exon1 = "ATGAAGCTAGTT"  # M-K-L-V
    
    # [핵심 실험 구간]
    if apply_mutation:
        # Pathogenic Mutation: Donor Site 파괴 (GT -> GC)
        # T를 C로 바꿈으로써 스플라이싱 신호가 깨짐
        intron_signal = "GCAAGT" 
        mutation_status = "🔴 Mutant (Splice Site Destroyed: T->C)"
    else:
        # Wild Type: 정상 Donor Site
        intron_signal = "GTAAGT"
        mutation_status = "🟢 Wild Type (Normal)"

    # 인트론 구성 (신호 + 무작위 서열 + AG)
    intron = intron_signal + generate_random_seq(20) + "AG"
    
    exon2 = "TGGCCATAA" # W-P-*
    
    # 전체 DNA 조립
    full_dna = generate_random_seq(10) + "GGCCGG" + \
               generate_random_seq(15) + promoter + \
               generate_random_seq(5) + \
               exon1 + intron + exon2 + \
               generate_random_seq(20)
               
    return full_dna, mutation_status

# ==========================================
# 3. 처리 엔진 (Regulation & Processing)
# ==========================================
def step1_transcription_regulation(dna_seq):
    if "TATAAA" in dna_seq:
        start_idx = dna_seq.find("TATAAA") + 6
        return dna_seq[start_idx:].replace('T', 'U')
    return None

def step2_splicing(pre_mrna):
    # 우리가 정한 강력한 Consensus Signal
    DONOR_SIGNAL = "GUAAGU" 
    
    intron_start = pre_mrna.find(DONOR_SIGNAL)
    
    if intron_start != -1:
        MIN_INTRON_LEN = 20
        search_start_pos = intron_start + MIN_INTRON_LEN
        intron_end = pre_mrna.find("AG", search_start_pos)
        
        if intron_end != -1:
            # 정상적으로 잘라냄
            mature_mrna = pre_mrna[:intron_start] + pre_mrna[intron_end+2:]
            print(f"  ✂️  Splicing Success! (Consensus '{DONOR_SIGNAL}' found)")
            return mature_mrna

    # [병리 상황] 신호를 못 찾으면? 
    # 인트론이 제거되지 않고 그대로 남음 (Intron Retention)
    print(f"  ⚠️  Splicing FAILED. Signal '{DONOR_SIGNAL}' not found.")
    print(f"  🚨  Intron Retention occurred!")
    return pre_mrna

def step3_translation(mature_mrna):
    start_idx = mature_mrna.find("AUG")
    if start_idx == -1: return ""
        
    protein_seq = []
    for i in range(start_idx, len(mature_mrna), 3):
        codon = mature_mrna[i:i+3]
        if len(codon) < 3: break
        
        aa = CODON_TABLE.get(codon.replace('U','T'), '?')
        if aa == '*': break
        protein_seq.append(aa)
        
    return "".join(protein_seq)

# ==========================================
# 4. 비교 실험 실행 (Comparative Study)
# ==========================================
def run_simulation(case_name, mutate):
    print(f"\n{'='*20} {case_name} {'='*20}")
    
    # 1. 유전자 생성
    my_gene, status = create_eukaryotic_gene(apply_mutation=mutate)
    print(f"Condition: {status}")
    
    # 2. 전사
    pre_mrna = step1_transcription_regulation(my_gene)
    
    # 3. 스플라이싱
    mature_mrna = step2_splicing(pre_mrna)
    
    # 4. 번역
    final_protein = step3_translation(mature_mrna)
    
    print(f"Resulting Protein: {final_protein}")
    print(f"Protein Length: {len(final_protein)} aa")
    
    return final_protein

if __name__ == "__main__":
    print("--- In Silico Mutagenesis Experiment ---")
    
    # Case 1: 정상인 (Wild Type)
    protein_wt = run_simulation("CASE 1: Normal Control", mutate=False)
    
    # Case 2: 환자 (Mutant Type)
    protein_mt = run_simulation("CASE 2: Patient (Mutation +1 T>C)", mutate=True)
    
    # 결과 비교 분석
    print(f"\n{'='*20} Final Clinical Report {'='*20}")
    print(f"1. Normal Protein: {protein_wt}")
    print(f"2. Mutant Protein: {protein_mt}")
    
    if protein_wt == protein_mt:
        print("Diagnosis: Benign Variant (No Phenotypic Change)")
    else:
        print("Diagnosis: Pathogenic Variant (Structure Altered)")
        print("Note: Splicing error caused a Frame-shift or Insertion.")