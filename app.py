import streamlit as st
import pandas as pd
import time
from Bio.Seq import Seq

# ==========================================
# 1. 기초 데이터 및 로직 설정
# ==========================================
# L858R Target Motif (Exon 21)
WT_SEQ_21 = "TTTGGGCTGGCCAAACTG"  # Normal (L)
MT_SEQ_21 = "TTTGGGCTGGCCAAACGG"  # Mutant (R)

# T790M Target Motif (Exon 20)
WT_SEQ_20 = "GGCACGGTG" # Normal (T)
MT_SEQ_20 = "GGCATGGTG" # Mutant (M)

def get_tumor_response(drug, has_l858r, has_t790m):
    """
    약물과 유전자 상태에 따른 종양 반응률(Response Rate) 반환
    반환값: (반응 유형, 월별 종양 크기 변화율)
    """
    if not has_l858r:
        return "No Effect", 1.05 # 암이 계속 자람 (5% 성장)

    if drug in ["Gefitinib", "Erlotinib"]:
        if has_t790m:
            return "Resistance", 1.1 # 내성으로 인해 급격히 성장 (10% 성장)
        else:
            return "Response", 0.7 # 치료 효과 좋음 (30% 감소)
    
    elif drug == "Osimertinib":
        return "Response", 0.6 # 강력한 치료 효과 (40% 감소)
        
    return "Unknown", 1.0

# ==========================================
# 2. Streamlit UI 구성
# ==========================================
st.set_page_config(page_title="EGFR Simulator", page_icon="🧬", layout="wide")

st.title("🧬 NSCLC: EGFR Mutation & Treatment Simulator")
st.markdown("---")

# [사이드바] 환자 상태 설정
st.sidebar.header("Patient Configuration")
patient_name = st.sidebar.text_input("Patient Name", "Mr. Kim")

st.sidebar.subheader("Genomic Profile")
l858r_status = st.sidebar.checkbox("L858R Mutation (Driver)", value=True)
t790m_status = st.sidebar.checkbox("T790M Mutation (Resistance)", value=False)

st.sidebar.subheader("Treatment Plan")
drug_choice = st.sidebar.selectbox("Select TKI Drug", ["None", "Gefitinib", "Erlotinib", "Osimertinib"])

# [메인 화면] 1. 유전자 서열 시각화
col1, col2 = st.columns(2)

with col1:
    st.subheader("🔍 Exon 21 Status (Driver)")
    if l858r_status:
        st.error(f"Mutant: ...{MT_SEQ_21}...")
        st.caption("Codon 858: CTG (L) -> **CGG (R)**")
    else:
        st.success(f"Wildtype: ...{WT_SEQ_21}...")
        st.caption("Codon 858: CTG (Leucine)")

with col2:
    st.subheader("🛡️ Exon 20 Status (Gatekeeper)")
    if t790m_status:
        st.error(f"Mutant: ...{MT_SEQ_20}...")
        st.caption("Codon 790: ACG (T) -> **ATG (M)**")
    else:
        st.success(f"Wildtype: ...{WT_SEQ_20}...")
        st.caption("Codon 790: ACG (Threonine)")

st.markdown("---")

# [메인 화면] 2. 치료 시뮬레이션 및 차트
st.subheader("📈 Clinical Course Simulation")

if st.button("Run Simulation"):
    # 시뮬레이션 데이터 생성
    months = list(range(0, 13)) # 0~12개월
    tumor_sizes = [100.0] # 초기 크기 100%
    
    response_type, rate = get_tumor_response(drug_choice, l858r_status, t790m_status)
    
    # 12개월치 데이터 계산
    current_size = 100.0
    for _ in range(12):
        if drug_choice == "None":
            current_size *= 1.05 # 치료 안하면 자연 성장
        else:
            current_size *= rate
        
        # 크기 제한 (0보다 작을 순 없음, 너무 커지면 사망 가정)
        current_size = max(0, min(current_size, 200))
        tumor_sizes.append(current_size)

    # 데이터프레임 변환
    chart_data = pd.DataFrame({
        "Month": months,
        "Tumor Burden (%)": tumor_sizes
    })

    # 결과 메시지 출력
    if response_type == "Response":
        st.success(f"✅ **[SUCCESS]** {drug_choice} is effective! Tumor is shrinking.")
    elif response_type == "Resistance":
        st.warning(f"⚠️ **[FAILURE]** Resistance detected. {drug_choice} is blocked by T790M.")
    elif response_type == "No Effect":
        st.error(f"❌ **[FAILURE]** Drug target (L858R) not found.")
    
    # 라인 차트 그리기
    st.line_chart(chart_data, x="Month", y="Tumor Burden (%)")
    
    # 상세 리포트
    with st.expander("See Molecular Pathology Report"):
        st.write(f"Patient: **{patient_name}**")
        st.write(f"Prescription: **{drug_choice}**")
        st.write(f"Molecular Profile: L858R ({l858r_status}), T790M ({t790m_status})")
        st.write("Conclusion: The graphical trend represents the expected tumor volume change based on the interaction between the TKI drug and the kinase domain structure.")

else:
    st.info("Adjust the settings on the sidebar and click 'Run Simulation'.")