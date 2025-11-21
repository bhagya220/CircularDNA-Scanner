import streamlit as st
from Bio import SeqIO, Seq
from difflib import SequenceMatcher
import io
import re
import matplotlib.pyplot as plt
from fpdf import FPDF

# ATTRACTIVE PNG IMAGE ICONS (change to your own for projects)
DNA_ICON = "https://static.thenounproject.com/png/1201536-200.png"  # PNG double helix
CIRCULAR_ICON = "https://cdn.pixabay.com/photo/2013/07/12/13/53/dna-147571_1280.png"  # PNG
PLASMID_ICON = "https://cdn-icons-png.flaticon.com/512/6167/6167011.png"  # PNG circular DNA

# MODERN, COLORFUL, STYLISH BACKGROUND AND WIDGETS
st.markdown("""
    <style>
    .stApp {
        background: linear-gradient(135deg, #22223b 0%, #4a4e69 40%, #9f86c0 100%) !important;
        min-height: 100vh;
    }
    .main > div {
        background: rgba(245,247,250,0.82) !important;
        border-radius: 20px !important;
        margin-top: 18px !important;
        padding: 12px 20px 12px 20px !important;
        box-shadow: 0 12px 40px 0 rgba(50,50,100,0.13) !important;
    }
    .stMarkdown h1, .stMarkdown h2, .stMarkdown h3, .stMarkdown h4 {
        color: #393e62 !important;
    }
    .stButton>button, .stDownloadButton>button {
        background: linear-gradient(90deg, #7267cb 15%, #5f72bd 75%);
        color: #f7fdff;
        border-radius: 6px;
        border: none;
        font-weight: 700;
        font-size: 1.08em;
        box-shadow: 0px 2px 18px  #aab4ed33;
        transition: 0.18s;
    }
    .stButton>button:hover, .stDownloadButton>button:hover {
        background: linear-gradient(90deg, #9f86c0, #5f72bd);
        color: #fff;
        border: none;
        transform: translateY(-2px) scale(1.03);
    }
    /* Sidebar style */
    .stSidebar {
        background: #393e62 !important;
        color: #fff !important;
    }
    /* Carded sidebar area (fallback for new Streamlit) */
    section[data-testid="stSidebar"] div[class^="css-"] {
        background: linear-gradient(160deg, #3d326b 0%, #7267cb 100%)!important;
        color: #fff !important;
        border-radius: 20px !important;
        margin: 8px;
        padding: 16px 8px;
    }
    /* Tab selector styling */
    .stTabs [data-baseweb="tab"] {
        color: #444c6e !important;
        font-weight: 500;
        font-size:1.11rem;
        background: #f4f5fa99;
        border-radius: 10px 10px 0 0;
        margin-right: 4px;
        margin-bottom: 0px;
        padding: 7px 20px 6px 20px !important;
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(90deg, #7267cbcc 0%, #c7b1e8cc 80%)!important;
        color: #1d1a29 !important;
        border-bottom: 4px solid #5271ff;
        font-weight: 700;
        box-shadow: 0 8px 22px -16px #888aad99;
    }
    img, .element-container img {
        background-color: #f6f5fd;
        border-radius: 14px;
        border: 2px solid #9f86c032;
        box-shadow: 0 2px 8px 1px #6c708fff;
        padding: 7px;
        margin-bottom: 4px;
        max-width: 100%;
        display: block;
        margin-left: auto;
        margin-right: auto;
    }
    </style>
""", unsafe_allow_html=True)

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
                continue
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
    ax.plot([0, 2 * 3.1416], [1, 1], color='black', linewidth=6)
    for r in repeats:
        theta = [2 * 3.1416 * r['start'] / length, 2 * 3.1416 * r['end'] / length]
        ax.plot(theta, [1.2, 1.2], color='red', linewidth=4)
    for o in orfs:
        theta = [2 * 3.1416 * o['start'] / length, 2 * 3.1416 * o['end'] / length]
        ax.plot(theta, [0.8, 0.8], color='green', linewidth=4)
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
    return pdf.output(dest='S').encode('latin-1')

# -- SIDEBAR WITH LOGO AND ABOUT --
with st.sidebar:
    st.image(DNA_ICON, width=180, use_column_width=True)
    st.markdown("### 🧬 Circular DNA Scanner")
    st.info(
        "Detects circularity, repeats, ORFs, and elements in DNA sequences. "
        "Visualize and export results. [GitHub](https://github.com/bhagya220/CircularDNA-Scanner) "
    )
    st.markdown("---")
    st.write("**Author:** bhagya220")
    st.write("Version: 1.0")
    st.image(PLASMID_ICON, width=110, caption="Plasmid example", use_column_width=False)

# -- MAIN PAGE TITLE WITH IMAGE --
st.set_page_config(page_title="Circular DNA Scanner", layout="wide")
col1, col2 = st.columns([7, 1])
with col1:
    st.title("🧬 Circular DNA Scanner")
with col2:
    st.image(CIRCULAR_ICON, width=90, use_column_width=False)

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
    c1, c2 = st.columns([2,1])
    with c1:
        st.header("Welcome to Circular DNA Scanner")
        st.write(
            "Scan DNA sequences for circularity, repeats, ORFs/genes, misassemblies, and annotations. "
            "Visualize and export your results."
        )
        st.markdown("""
        - **Upload:** Upload FASTA files or paste sequences.
        - **Result:** Circularity, repeats, ORFs/genes, misassemblies, annotations.
        - **Visualization:** Circular genome map.
        - **Documentation:** How the tool works.
        - **PDF Report:** Export results as PDF.
        """)
    with c2:
        st.image(DNA_ICON, width=150, caption="DNA double helix", use_column_width=True)
        st.image(PLASMID_ICON, width=90, caption="Circular DNA", use_column_width=True)

with tabs[1]:
    st.header("Upload or Paste Sequence")
    st.image(CIRCULAR_ICON, width=70, use_column_width=False)
    min_overlap = st.slider("Minimum Overlap (bases)", min_value=20, max_value=500, value=min_overlap_default, step=10)
    identity_threshold = st.slider("Identity Threshold", min_value=0.5, max_value=1.0, value=identity_threshold_default, step=0.01)
    upload_option = st.radio("Choose input method:", ["Upload FASTA file", "Paste Sequence"])
    results = []
    if upload_option == "Upload FASTA file":
        uploaded_file = st.file_uploader("Choose a FASTA file", type=["fa", "fasta", "txt"])
        if st.button("Download example FASTA"):
            example_fasta = (
                ">Example1\nATGCGTACGTTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA" +
                "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA" +
                "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
            )
            st.download_button("Click to download example FASTA", example_fasta, file_name="example.fasta")
        if uploaded_file is not None:
            fasta_io = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
            with st.spinner("Processing uploaded sequences..."):
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
                st.balloons()
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
            c1, c2 = st.columns([4,1])
            with c1:
                st.subheader(f"Results for {result['Sequence ID']}")
                st.write(f"**Length:** {result['Length (bp)']} bp")
                st.write(f"**Circular?** {'🟢 Yes' if result['Circular?'] else '🔴 No'} (Similarity: {result['Similarity']})")
                st.write(f"**Repeats Detected:** {len(result['Repeats'])}")
                st.write(f"**ORFs Detected:** {len(result['ORFs'])}")
                st.write(f"**Annotations:** {', '.join([a['element'] for a in result['Annotations']]) or 'None'}")
                st.code(result["Sequence"][:300] + "...", language="text")
            with c2:
                st.image(PLASMID_ICON, width=80, use_column_width=True)
    else:
        st.info("No results to show. Please process sequences in the Upload tab.")

with tabs[3]:
    st.header("Visualization")
    st.image(PLASMID_ICON, width=70, use_column_width=False)
    results = st.session_state.get('results', [])
    if results:
        for result in results:
            st.subheader(f"Circular Map for {result['Sequence ID']}")
            fig = plot_circular_map(result["Sequence"], result["Repeats"], result["ORFs"], result["Annotations"])
            st.pyplot(fig)
    else:
        st.info("No results to visualize. Please process sequences in the Upload tab.")

with tabs[4]:
    st.header("Documentation")
    st.image(DNA_ICON, width=90, use_column_width=False)
    st.markdown("""
    ## About Circular DNA Scanner
    This tool predicts whether DNA sequences are circular by comparing the first and last N bases for similarity.
    It also detects simple repeats, possible misassemblies, predicts ORFs/genes, annotates known elements, and visualizes the genome in a circular map.

    ### Features
    - **🧬 Circularity prediction:** End overlap similarity.
    - **🔄 Repeat detection:** Identifies repeated motifs.
    - **🟢 ORF/gene prediction:** Finds open reading frames.
    - **🔷 Annotation:** Detects TATA box, polyA signal, etc.
    - **🖼️ Visualization:** Genome map.
    - **📄 PDF report:** Export results.

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
        st.image("https://cdn-icons-png.flaticon.com/512/337/337946.png", width=48, caption="Download your PDF!", use_column_width=False)
    else:
        st.info("No results to export. Please process sequences in the Upload tab.")
