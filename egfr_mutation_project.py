from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# ==========================================
# 1. 실제 데이터 가져오기 (Real-world Data) --> EGFR Exon 21 
# ==========================================
# NCBI Reference Sequence: NM_005228.5 (EGFR)의 Exon 21 핵심 부위 (Coding Sequence)
# 실제 인간 EGFR 유전자의 2560번째 염기부터의 서열입니다.
# Codon 858 (L) 주변 서열: ...AAA GAT CAA CTG AAT... (K-D-Q-L-N)

EGFR_EXON_21_WT_SEQ = "AGCAAC AAG GAA ATC CTC GAT GAA GCC TAC GTG ATG GCC AGC GTG GAC AAC CCC CAC GTG TGC CGC CTG CTG GGC ATC TGC CTC ACC TCC ACC GTG CAG CTC ATC ACG CAG CTC ATG CCC TTC GGC TGC CTC CTG GAC TAT GTC CGG GAA CAC AAA GAC AAT ATT GGC TCC CAG TAC CTG CTC AAC TGG TGT GTG CAG ATC GCA AAG GGC ATG AAC TAC TTG GAG GAC CGT CGC TTG GTG CAC CGC GAC CTG GCA GCC AGG AAC GTA CTG GTG AAA ACA CCG CAG CAT GTC AAG ATC ACA GAT TTT GGG CTG GCC AAA CTG CTG GGT GCG GAA GAG AAA GAA TAC CAT GCA GAA GGA GGC AAA GTG CCT ATC AAG TGG ATG GCA TTG GAA TCA ATT TTA CAC AGA ATC TAT ACC CAC CAG AGT GAT GTC TGG AGC TAC GGG GTG ACC GTT TGG GAG TTG ATG ACC TTT GGA TCC AAG CCA TAT GAC GGA ATC CCT GCC AGC GAG ATC TCC TCC ATC CTG GAG AAA GGA GAA CGC CTC CCT CAG CCA CCC ATA TGT ACC ATC GAT GTC TAC ATG ATC ATG GTC AAG TGC TGG ATG ATA GAC GCA GAT AGT CGC CCA AAG TTC CGT GAG TTG ATC ATC GAA TTC TCC AAA ATG GCC CGA GAC CCC CAG CGC TAC CTT GTC ATT CAG GGG GAT GAA AGA ATG CAT TTG CCA AGT CCT ACA GAC TCC AAC TTC TAC CGT GCC CTG ATG GAT GAA GAA GAC ATG GAC GAC GTG GTG GAT GCC GAC GAG TAC CTC ATC CCA CAG CAG GGC TTC TTC AGC AGC CCC TCC ACG TCA CGG ACT CCC CTC CTG AGC TCT CTG AGT GCA ACC AGC AAC AAT TCC ACC GTG GCT TGC ATT GAT AGA AAT GGG CTG CAA AGC TGT CCC ATC AAG GAA GAC AGC TTC TTG CAG CGA TAC AGC TCA GAC CCC ACA GGC GCC TTG ACT GAG GAC AGC ATA GAC GAC ACC TTC CTC CCA GTG CCT GAA TAC ATA AAC CAG TCC GTT CCC AAA AGG CCC GCT GGC TCT GTG CAG AAT CCT GTC TAT CAC AAT CAG CCT CTG AAC CCC GCG CCC AGC AGA GAC CCA CAC TAC CAG GAC CCC CAC AGC ACT GCA GTG GGC AAC CCC GAG TAT CTC AAC ACT GTC CAG CCC ACC TGT GTC AAC AGC ACA TTC GAC AGC CCT GCC CAC TGG GCC CAG AAA GGC AGC CAC CAA ATT AGC CTG GAC AAC CCT GAC TAC CAG CAG GAC TTC TTT CCC AAG GAA GCC AAG CCA AAT GGC ATC TTT AAG GGC TCC ACA GCT GAA AAT GCA GAA TAC CTA AGG GTC GCG CCA CAA AGC AGT GAA TTT ATT GGA GCA TGA".replace(" ", "")

# BioPython 객체로 변환
wild_type_dna = Seq(EGFR_EXON_21_WT_SEQ)

# ==========================================
# 2. 돌연변이 유발 (Mutagenesis: L858R)
# ==========================================
def induce_l858r_mutation(dna_seq):
    """
    EGFR L858R Mutation:
    Codon 858: CTG (Leucine) -> CGG (Arginine)
    Nucleotide change: c.2573 T > G
    """
    # Exon 21 내에서 L858 코돈(CTG)의 위치를 찾습니다.
    # 위 서열에서 'TTT GGG CTG GCC AAA CTG CTG' 패턴 중 두 번째 CTG가 858번입니다.
    # 정확한 인덱싱을 위해 찾기 (실제 연구에선 좌표가 주어짐)
    
    target_motif = "TTTGGGCTGGCCAAACTG" # F-G-L-A-K-L (853~858)
    motif_index = dna_seq.find(target_motif)
    
    if motif_index == -1:
        print("Error: Target motif not found!")
        return None
        
    # L858의 CTG 위치 계산 (Motif 끝부분)
    l858_codon_index = motif_index + len(target_motif) - 3 # 마지막 CTG
    
    # 돌연변이 전 확인
    original_codon = dna_seq[l858_codon_index : l858_codon_index+3]
    print(f"\n[Target Validation]")
    print(f"  - Loc: Index {l858_codon_index}")
    print(f"  - Codon: {original_codon} (Should be CTG)")
    
    # 돌연변이 생성 (T -> G)
    # CTG -> CGG
    # 문자열은 불변(Immutable)이므로 리스트로 변환 후 수정
    dna_list = list(str(dna_seq))
    dna_list[l858_codon_index + 1] = "G" # 중간의 T를 G로 변경
    mutant_seq = Seq("".join(dna_list))
    
    return mutant_seq

# ==========================================
# 3. 실행 및 비교 (Execution)
# ==========================================
if __name__ == "__main__":
    print("--- Clinical Genomics: EGFR L858R Analysis ---")
    
    # 1. Wild Type (정상) 분석
    print("\n1. Wild Type (Normal EGFR)")
    wt_protein = wild_type_dna.translate()
    print(f"  - DNA Length: {len(wild_type_dna)} bp")
    print(f"  - Protein Segment: ...{wt_protein[130:145]}...") # L858 주변 보여주기
    
    # 2. Mutant Type (폐암) 생성
    mutant_dna = induce_l858r_mutation(wild_type_dna)
    
    if mutant_dna:
        print("\n2. Mutant Type (L858R Positive)")
        mt_protein = mutant_dna.translate()
        
        # 3. 단백질 비교 (L -> R 확인)
        # 위에서 motif index를 기반으로 대략 130~140번째 아미노산 쯤에 위치함
        # 정확한 비교를 위해 차이점 출력
        
        for i, (aa1, aa2) in enumerate(zip(wt_protein, mt_protein)):
            if aa1 != aa2:
                print(f"\n[!!! MUTATION DETECTED !!!]")
                print(f"  - Position: Amino Acid {i+1} (Local Count)")
                print(f"  - Change:   {aa1} (Leucine) --> {aa2} (Arginine)")
                print(f"  - Clinical: This is the 'L858R' driver mutation.")
                
        # 4. 파일 저장 (SnapGene용)
        # FASTA 파일로 저장해서 시각화 준비
        records = [
            SeqRecord(wild_type_dna, id="EGFR_Exon21_WT", description="Homo sapiens EGFR exon 21 partial cds"),
            SeqRecord(mutant_dna, id="EGFR_Exon21_L858R", description="Point mutation c.2573T>G")
        ]
        SeqIO.write(records, "EGFR_L858R_Analysis.fasta", "fasta")
        print(f"\n📂 File saved: 'EGFR_L858R_Analysis.fasta'")
        print("   -> Open this file in SnapGene Viewer to visualize the T>G change!")