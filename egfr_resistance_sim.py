from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# ==========================================
# 1. 유전자 데이터 셋팅 (Exon 20 & 21) --> T790M & L858R
# ==========================================

# Exon 20 부분 서열 (T790 포함)
# Normal Codon 790: ACG (Threonine)
# Context: ...GGC ACG GTG... (Gly-Thr-Val)
EXON_20_WT_SEQ = "ACTTCACAGCCCTGCATGTGTTCCTCACAGGAGAGAGCCCTGTTGGCTGCAGTGAGCAGATGTTTGGGAGCCAAGCACTGCCAGGTGCCAGCCCTGTGCGTTCTGGTGGCAGAGGGAGTCGTGGCTGCGTTGGCACTGCGTTCCTCCTGCTGGCCCGTGTGCCGCTGCAGGGTGAGGTGCAGCTGGTGGACGCCGCAGTCGGCCTTGAGCTTGCGTTTGTTCTTCACCGTGCGTTCGGCACGGTGTATAAGTAAGCAGCCTCTGTTCTGCTGC"

# Exon 21 부분 서열 (L858 포함 - 이전 코드 사용)
EXON_21_WT_SEQ = "AGCAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGAAGAGAAAGAA"

# ==========================================
# 2. 환자 클래스 정의 (임상 시뮬레이터)
# ==========================================
class CancerPatient:
    def __init__(self, name):
        self.name = name
        self.exon20 = Seq(EXON_20_WT_SEQ)
        self.exon21 = Seq(EXON_21_WT_SEQ)
        self.status = "Healthy"
        
    def check_biomarkers(self):
        """현재 유전자 상태를 분석하여 리포트 (수정됨)"""
        
        # 1. Exon 21: L858R 확인 (Stable Anchor 방식)
        ANCHOR_MOTIF = "TTTGGGCTGGCCAAA"
        idx = self.exon21.find(ANCHOR_MOTIF)
        
        has_l858r = False
        if idx != -1:
            # Anchor 바로 뒤의 3글자(코돈)를 읽어옵니다.
            target_codon = str(self.exon21[idx+15 : idx+18])
            if target_codon == "CGG": # Arginine
                has_l858r = True
        
        # 2. Exon 20: T790M 확인
        # T790M은 변이 자체가 새로운 서열(GGCATGGTG)을 만들기 때문에 기존 방식도 작동하지만,
        # 정확성을 위해 찾기 방식을 유지합니다.
        t790m_motif = "GGCATGGTG" # Mutated motif
        has_t790m = (self.exon20.find(t790m_motif) != -1)
            
        return has_l858r, has_t790m

    # ... (apply_drug_treatment 함수는 아까 수정한 버전 그대로 사용) ...
    def apply_drug_treatment(self, drug_name):
        l858r, t790m = self.check_biomarkers()
        
        print(f"\n💉 Treating {self.name} with [{drug_name}]...")
        print(f"   - Genotype: L858R={'Positive' if l858r else 'Negative'} | T790M={'Positive' if t790m else 'Negative'}")
        
        if not l858r:
            print("   -> Result: No effect (Target mutation L858R not found).")
            return

        # 약물 반응 로직
        if drug_name in ["Gefitinib", "Erlotinib"]: # (First-Generation TKI)
            if t790m:
                print("   -> ⚠️ Result: FAILED (Resistance). Tumor Progression.")
                print("      Reason: T790M mutation blocks the drug binding pocket (Steric Hindrance).")
            else:
                print("   -> ✅ Result: SUCCESS (Partial/Complete Response). Tumor Shrinkage.")
                print("      Reason: Drug effectively binds to L858R mutant kinase.")
                
        elif drug_name in ["Osimertinib"]: # (Third-Generation TKI)
            if t790m:
                print("   -> ✅ Result: SUCCESS. Tumor Shrinkage.")
                print("      Reason: Osimertinib can bind even with T790M mutation.")
            else:
                print("   -> ✅ Result: SUCCESS. Effective against L858R.")

# ==========================================
# 3. 돌연변이 유발 함수
# [수정 1] 돌연변이 유발 함수 (Stable Motif 사용)
# ==========================================
def induce_l858r(seq):
    # 변하지 않는 고정 서열(Anchor)을 사용하여 위치를 찾습니다.
    # L858 코돈(CTG) 바로 앞의 15개 염기: F-G-L-A-K (TTT GGG CTG GCC AAA)
    ANCHOR_MOTIF = "TTTGGGCTGGCCAAA" 
    
    idx = seq.find(ANCHOR_MOTIF)
    if idx == -1: 
        return seq
        
    mutable_seq = list(str(seq))
    
    # Anchor 끝(idx + 15)부터가 코돈 시작입니다.
    # CTG -> CGG (가운데 T를 G로 변경)
    # 위치: idx + 15(C), idx + 16(T), idx + 17(G)
    mutation_pos = idx + 16
    mutable_seq[mutation_pos] = "G" 
    
    return Seq("".join(mutable_seq))

def induce_t790m(seq):
    # T790 (ACG) -> M790 (ATG)
    # Exon 20에서 'GGC ACG GTG' 찾기
    motif = "GGCACGGTG"
    idx = seq.find(motif)
    if idx == -1: return seq
    
    mutable_seq = list(str(seq))
    # ACG의 가운데 C(index+4)를 T로 변경 -> ATG
    # GGC A[C]G GTG -> GGC A[T]G GTG
    mutable_seq[idx + 4] = "T"
    return Seq("".join(mutable_seq))

# ==========================================
# 4. 시뮬레이션 실행 (Time-course Simulation)
# ==========================================
if __name__ == "__main__":
    print("--- Lung Cancer Treatment Simulation (EGFR L858R & T790M) ---")
    
    # 1. 환자 입장
    patient = CancerPatient("Mr. Kim")
    
    # [Stage 1] 암 발병 (L858R 발생)
    print("\n[Stage 1] Diagnosis: Lung Cancer Detected")
    patient.exon21 = induce_l858r(patient.exon21) # 돌연변이 발생
    
    # 1세대 표적항암제 투여
    patient.apply_drug_treatment("Gefitinib")
    
    # ... 시간 경과 (10~14개월 후) ...
    print("\n... 12 Months Later ...")
    
    # [Stage 2] 내성 획득 (Acquired Resistance: T790M 발생)
    print("\n[Stage 2] Disease Progression (Recurrence)")
    print("   🧬 Evolution: Cancer cells acquired T790M 'Gatekeeper' mutation.")
    patient.exon20 = induce_t790m(patient.exon20) # 내성 돌연변이 발생
    
    # 1세대 항암제 다시 투여 (실패 예측)
    patient.apply_drug_treatment("Gefitinib")
    
    # [Stage 3] 3세대 항암제 변경 (Osimertinib)
    print("\n[Stage 3] Treatment Switch")
    patient.apply_drug_treatment("Osimertinib")

    # [Bonus] 파일 저장 for SnapGene
    print("\n[Saving Data]")
    records = [
        SeqRecord(patient.exon20, id="EGFR_Exon20_T790M", description="Resistance Mutation c.2369C>T"),
        SeqRecord(patient.exon21, id="EGFR_Exon21_L858R", description="Sensitizing Mutation c.2573T>G")
    ]
    SeqIO.write(records, "EGFR_Resistance_Profile.fasta", "fasta")
    print("📂 'EGFR_Resistance_Profile.fasta' saved. Check T790M in Exon 20!")