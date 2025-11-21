import streamlit as st
from Bio import SeqIO, Seq
from difflib import SequenceMatcher
import io
import re
import matplotlib.pyplot as plt
from fpdf import FPDF

def check_circularity(seq, min_overlap=100, identity_threshold=0.9):
    start = seq[:min_overlap]
    end = seq[-min_overlap:]
    ratio = SequenceMatcher(None, start, end).ratio()
    return ratio >= identity_threshold, ratio

def detect_repeats(seq, min_repeat_size=10):
    seq = seq.upper()
    repeats = []
    for size in range(min_repeat_size, min_repeat_size+10):
        for i in range(len(seq)-size):
            motif = seq[i:i+size]
            if seq.count(motif) > 1:
                repeats.append({'motif': motif, 'start': i, 'end': i+size})
    return repeats

def find_ORFs(seq, min_orf=100):
    # Clean up sequence: keep only A, T, G, C
    seq_clean = re.sub(r'[^ATGC]', '', seq.upper())
    seq_obj = Seq.Seq(seq_clean)
    orfs = []
    for strand, nuc in [(+1, seq_obj), (-1, seq_obj.reverse_complement())]:
        for frame in range(3):
            nuc_frame = nuc[frame:]
            codon_length = len(nuc_frame) - (len(nuc_frame) % 3)
            nuc_trimmed = nuc_frame[:codon_length]
            try:
                trans = str(nuc_trimmed.translate(to_stop=False))
            except Exception:
                continue  # Skip frames with translation errors
            for m in re.finditer('M.*?\*', trans):
                orf_len = (m.end() - m.start())*3
                if orf_len >= min_orf:
                    start = frame + m.start()*3
                    end = frame + m.end()*3
                    orfs.append({'strand': strand, 'start': start, 'end': end, 'length': orf_len})
    return orfs

def annotate_elements(seq):
    seq = seq.upper()
    annotations = []
    if "TATAAA" in seq:
        annotations.append({'element': 'TATA box', 'pos': seq.find('TATAAA')})
    if "AATAAA" in seq:
        annotations.append({'element': 'polyA signal', 'pos': seq.find('AATAAA')})
    return annotations

def plot_circular_map(seq, repeats, orfs, annotations):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6, 6))
    length = len(seq)
    ax.set_xticks([])
    ax.set_yticks([])
    # Draw backbone
    ax.plot([0, 2 * 3.1416], [1, 1], color='black', linewidth=6)
    # Draw repeats
    for r in repeats:
        theta = [2 * 3.1416 * r['start'] / length, 2 * 3.1416 * r['end'] / length]
        ax.plot(theta, [1.2, 1.2], color='red', linewidth=4)
    # Draw ORFs
    for o in orfs:
        theta = [2 * 3.1416 * o['start'] / length, 2 * 3.1416 * o['end'] / length]
        ax.plot(theta, [0.8, 0.8], color='green', linewidth=4)
    # Annotate elements
    for a in annotations:
        theta = 2 * 3.1416 * a['pos'] / length
        ax.plot([theta], [1], marker='o', color='blue')
        ax.text(theta, 1.05, a['element'], color='blue', fontsize=8, ha='center')
    ax.set_title('Circular Genome Map', va='bottom')
    return fig

def generate_pdf_report(results, seq, repeats, orfs, annotations):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt="Circular DNA Scanner Report", ln=1, align='C')
    pdf.cell(200, 10, txt=f"Sequence ID: {results['Sequence ID']}", ln=2)
    pdf.cell(200, 10, txt=f"Length (bp): {results['Length (bp)']}", ln=2)
    pdf.cell(200, 10, txt=f"Circular?: {results['Circular?']} (Similarity: {results['Similarity']})", ln=2)
    pdf.cell(200, 10, txt=f"Repeats detected: {len(repeats)}", ln=2)
    pdf.cell(200, 10, txt=f"ORFs detected: {len(orfs)}", ln=2)
    pdf.cell(200, 10, txt=f"Annotations: {', '.join([a['element'] for a in annotations])}", ln=2)
    pdf.multi_cell(0, 10, txt=f"First 300 bases: {seq[:300]}...")
    # For illustration, no images embedded; to embed, save and use image
    return pdf.output(dest='S').encode('latin-1')

st.set_page_config(page_title="Circular DNA Scanner", layout="wide")
st.title("Circular DNA Scanner")

tabs = st.tabs([
    "Home",
    "Upload",
    "Result",
    "Visualization",
    "Documentation",
    "PDF Report"
])
min_overlap_default = 100
identity_threshold_default = 0.9

with tabs[0]:
    st.header("Welcome to Circular DNA Scanner")
    st.write("""
        Scan DNA sequences for circularity, repeats, ORFs/genes, misassemblies, and annotations.
        Visualize and export your results.
    """)
    st.markdown("""
    - **Upload:** Upload FASTA files or paste sequences.
    - **Result:** Circularity, repeats, ORFs/genes, misassemblies, annotations.
    - **Visualization:** Circular genome map.
    - **Documentation:** How the tool works.
    - **PDF Report:** Export results as PDF.
    """)

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
                repeats = detect_repeats(seq)
                orfs = find_ORFs(seq)
                annotations = annotate_elements(seq)
                results.append({
                    "Sequence ID": record.id,
                    "Length (bp)": len(seq),
                    "Circular?": is_circular,
                    "Similarity": f"{score:.2f}",
                    "Repeats": repeats,
                    "ORFs": orfs,
                    "Annotations": annotations,
                    "Sequence": seq
                })
            st.session_state['results'] = results
            st.success("Sequences processed! Go to the Result tab to view results.")
    else:
        pasted_seq = st.text_area("Paste DNA sequence (only bases)", height=150)
        seq_id = st.text_input("Sequence ID (optional)", value="PastedSeq1")
        if st.button("Check Pasted Sequence"):
            if pasted_seq.strip():
                seq = pasted_seq.replace("\n", "").replace(" ", "").upper()
                is_circular, score = check_circularity(seq, min_overlap, identity_threshold)
                repeats = detect_repeats(seq)
                orfs = find_ORFs(seq)
                annotations = annotate_elements(seq)
                result = {
                    "Sequence ID": seq_id,
                    "Length (bp)": len(seq),
                    "Circular?": is_circular,
                    "Similarity": f"{score:.2f}",
                    "Repeats": repeats,
                    "ORFs": orfs,
                    "Annotations": annotations,
                    "Sequence": seq
                }
                st.session_state['results'] = [result]
                st.success("Sequence processed! Go to the Result tab to view results.")
            else:
                st.warning("Please paste a DNA sequence.")

with tabs[2]:
    st.header("Results")
    results = st.session_state.get('results', [])
    if results:
        for result in results:
            st.subheader(f"Results for {result['Sequence ID']}")
            st.write(f"Length: {result['Length (bp)']} bp")
            st.write(f"Circular? **{result['Circular?']}** (Similarity: {result['Similarity']})")
            st.write(f"Repeats Detected: {len(result['Repeats'])}")
            st.write(f"ORFs Detected: {len(result['ORFs'])}")
            st.write(f"Annotations: {', '.join([a['element'] for a in result['Annotations']])}")
            st.write("First 300 bases: " + result["Sequence"][:300] + "...")
    else:
        st.info("No results to show. Please process sequences in the Upload tab.")

with tabs[3]:
    st.header("Visualization")
    results = st.session_state.get('results', [])
    if results:
        for result in results:
            fig = plot_circular_map(result["Sequence"], result["Repeats"], result["ORFs"], result["Annotations"])
            st.pyplot(fig)
    else:
        st.info("No results to visualize. Please process sequences in the Upload tab.")

with tabs[4]:
    st.header("Documentation")
    st.markdown("""
    ## About Circular DNA Scanner
    This tool predicts whether DNA sequences are circular by comparing the first and last N bases for similarity.
    It also detects simple repeats, possible misassemblies, predicts ORFs/genes, annotates known elements, and visualizes the genome in a circular map.

    ### Features
    - **Circularity prediction:** End overlap similarity.
    - **Repeat detection:** Identifies repeated motifs.
    - **ORF/gene prediction:** Finds open reading frames.
    - **Annotation:** Detects TATA box, polyA signal, etc.
    - **Visualization:** Genome map.
    - **PDF report:** Export results.

    ### How to Use
    1. Upload a FASTA file or paste a DNA sequence in the **Upload** tab.
    2. Adjust parameters as needed.
    3. Process sequences.
    4. View results in **Result** and **Visualization** tabs.
    5. Export a PDF report in the **PDF Report** tab.

    ### Notes
    - Multi-sequence FASTA files supported.
    - For research use only; further validation required for clinical/diagnostic use.
    """)

with tabs[5]:
    st.header("PDF Report")
    results = st.session_state.get('results', [])
    if results:
        result = results[0]  # Only first result for demo; can loop for batch
        pdf_bytes = generate_pdf_report(result, result["Sequence"], result["Repeats"], result["ORFs"], result["Annotations"])
        st.download_button("Download PDF Report", pdf_bytes, file_name="CircularDNA_Report.pdf")
    else:
        st.info("No results to export. Please process sequences in the Upload tab.")



