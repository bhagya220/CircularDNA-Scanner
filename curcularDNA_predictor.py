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
st.write("Upload a FASTA file to scan for circular DNA.")

uploaded_file = st.file_uploader("Choose a FASTA file", type=["fa", "fasta", "txt"])

if uploaded_file is not None:
    fasta_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    results = []
    for record in SeqIO.parse(fasta_io, "fasta"):
        seq = str(record.seq).upper()
        is_circular, score = check_circularity(seq)
        results.append({
            "Sequence ID": record.id,
            "Length (bp)": len(seq),
            "Circular?": is_circular,
            "Similarity": f"{score:.2f}"
        })
    st.write("## Results")
    st.table(results)
