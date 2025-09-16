import streamlit as st
from Bio import SeqIO
from difflib import SequenceMatcher
import io

def check_circularity(seq, min_overlap=100, identity_threshold=0.9):
    start = seq[:min_overlap]
    end = seq[-min_overlap:]
    ratio = SequenceMatcher(None, start, end).ratio()
    return ratio >= identity_threshold, ratio

# App Title
st.set_page_config(page_title="Circular DNA Scanner", layout="wide")
st.title("Circular DNA Scanner")

# --- Parameter controls (shown on relevant tabs) ---
min_overlap_default = 100
identity_threshold_default = 0.9

# Tabs: Home, Upload, Result, Visualization, Documentation
tabs = st.tabs(["Home", "Upload", "Result", "Visualization", "Documentation"])

# ---- HOME TAB ----
with tabs[0]:
    st.header("Welcome to Circular DNA Scanner")
    st.write("""
        This tool allows you to scan DNA sequences for circularity. You can upload a FASTA file containing one or more sequences,
        or paste a sequence directly. The scanner will compare the ends of each sequence to check for circularity.
    """)
    st.markdown("""
    - **Upload:** Upload FASTA files or paste sequences.
    - **Result:** See the circularity prediction results for your sequences.
    - **Visualization:** Visualize sequence alignments and similarity.
    - **Documentation:** Learn how the tool works and how to use it.
    """)

# ---- UPLOAD TAB ----
with tabs[1]:
    st.header("Upload or Paste Sequence")
    min_overlap = st.slider("Minimum Overlap (bases)", min_value=20, max_value=500, value=min_overlap_default, step=10)
    identity_threshold = st.slider("Identity Threshold", min_value=0.5, max_value=1.0, value=identity_threshold_default, step=0.01)

    upload_option = st.radio("Choose input method:", ["Upload FASTA file", "Paste Sequence"])
    results = []

    if upload_option == "Upload FASTA file":
        uploaded_file = st.file_uploader("Choose a FASTA file", type=["fa", "fasta", "txt"])
        if st.button("Download example FASTA"):
            example_fasta = ">Example1\nATGCGTACGTTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
            st.download_button("Click to download example FASTA", example_fasta, file_name="example.fasta")
        if uploaded_file is not None:
            fasta_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
            for record in SeqIO.parse(fasta_io, "fasta"):
                seq = str(record.seq).upper()
                is_circular, score = check_circularity(seq, min_overlap, identity_threshold)
                results.append({
                    "Sequence ID": record.id,
                    "Length (bp)": len(seq),
                    "Circular?": is_circular,
                    "Similarity": f"{score:.2f}"
                })
            st.session_state['results'] = results
            st.success("Sequences processed! Go to the Result tab to view results.")

    else:  # Paste Sequence
        pasted_seq = st.text_area("Paste DNA sequence (only bases)", height=150)
        seq_id = st.text_input("Sequence ID (optional)", value="PastedSeq1")
        if st.button("Check Pasted Sequence"):
            if pasted_seq.strip():
                seq = pasted_seq.replace("\n", "").replace(" ", "").upper()
                is_circular, score = check_circularity(seq, min_overlap, identity_threshold)
                result = {
                    "Sequence ID": seq_id,
                    "Length (bp)": len(seq),
                    "Circular?": is_circular,
                    "Similarity": f"{score:.2f}"
                }
                st.session_state['results'] = [result]
                st.success("Sequence processed! Go to the Result tab to view results.")
            else:
                st.warning("Please paste a DNA sequence.")

# ---- RESULT TAB ----
with tabs[2]:
    st.header("Results")
    results = st.session_state.get('results', [])
    if results:
        st.write("## Circularity Prediction Results")
        st.table(results)
    else:
        st.info("No results to show. Please process sequences in the Upload tab.")

# ---- VISUALIZATION TAB ----
with tabs[3]:
    st.header("Visualization")
    results = st.session_state.get('results', [])
    if results:
        for result in results:
            st.subheader(f"Visualization for {result['Sequence ID']}")
            st.write(f"Sequence Length: {result['Length (bp)']} bp")
            st.write(f"Circular? **{result['Circular?']}** (Similarity: {result['Similarity']})")
            # Simple visualization: similarity bar
            similarity = float(result['Similarity'])
            st.progress(similarity)
    else:
        st.info("No results to visualize. Please process sequences in the Upload tab.")

# ---- DOCUMENTATION TAB ----
with tabs[4]:
    st.header("Documentation")
    st.markdown("""
    ## About Circular DNA Scanner
    This tool predicts whether DNA sequences are circular by comparing the first and last N bases for similarity.
    - **Minimum Overlap:** Number of bases compared at sequence ends.
    - **Identity Threshold:** Required similarity ratio to consider sequence circular.

    ### How to Use
    1. Go to the **Upload** tab and upload a FASTA file or paste a DNA sequence.
    2. Adjust parameters as needed.
    3. Process your sequence(s).
    4. View results in the **Result** tab.
    5. See similarity visualization in the **Visualization** tab.

    ### Output
    - **Sequence ID:** Name of the sequence.
    - **Length (bp):** Sequence length in base pairs.
    - **Circular?:** Whether the sequence is predicted to be circular.
    - **Similarity:** Similarity score between sequence ends.

    ### Notes
    - Multi-sequence FASTA files are supported.
    - Adjust parameters for sensitivity.
    - Contact the developer for help or feature requests.
    """)
