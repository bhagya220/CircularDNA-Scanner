import streamlit as st
from Bio import SeqIO
from difflib import SequenceMatcher
import io

def check_circularity(seq, min_overlap=100, identity_threshold=0.9):
    start = seq[:min_overlap]
    end = seq[-min_overlap:]
    ratio = SequenceMatcher(None, start, end).ratio()
    return ratio >= identity_threshold, ratio

st.title("Circular DNA Scanner")
st.write("Scan DNA sequences for circularity by uploading a FASTA file or pasting a sequence.")

# --- Parameter controls ---
min_overlap = st.slider("Minimum Overlap (bases)", min_value=20, max_value=500, value=100, step=10)
identity_threshold = st.slider("Identity Threshold", min_value=0.5, max_value=1.0, value=0.9, step=0.01)

tab1, tab2 = st.tabs(["Upload FASTA file", "Paste Sequence"])

with tab1:
    uploaded_file = st.file_uploader("Choose a FASTA file", type=["fa", "fasta", "txt"])
    if st.button("Download example FASTA"):
        example_fasta = ">Example1\nATGCGTACGTTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATGCGTACGTTAGCTAGC\n"
        st.download_button("Click to download example FASTA", example_fasta, file_name="example.fasta")

    if uploaded_file is not None:
        fasta_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        results = []
        for record in SeqIO.parse(fasta_io, "fasta"):
            seq = str(record.seq).upper()
            is_circular, score = check_circularity(seq, min_overlap, identity_threshold)
            results.append({
                "Sequence ID": record.id,
                "Length (bp)": len(seq),
                "Circular?": is_circular,
                "Similarity": f"{score:.2f}"
            })
        st.write("## Results")
        st.table(results)

with tab2:
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
            st.write("## Result for Pasted Sequence")
            st.table([result])
        else:
            st.warning("Please paste a DNA sequence.")

st.markdown("---")
st.info("Tip: Adjust overlap and threshold to change sensitivity. You can upload multi-sequence FASTA files.")
