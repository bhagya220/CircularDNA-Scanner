print ("This tool will Scan cDNA in genomes")


from Bio import SeqIO
from difflib import SequenceMatcher
import argparse

def check_circularity(seq, min_overlap=100, identity_threshold=0.9):
    """
    Simple circularity check:
    Compare first N bases with last N bases.
    If similarity >= threshold → considered circular.
    """
    start = seq[:min_overlap]
    end = seq[-min_overlap:]
    ratio = SequenceMatcher(None, start, end).ratio()
    return ratio >= identity_threshold, ratio

def main(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        is_circular, score = check_circularity(seq)
        print(f"Sequence: {record.id}")
        print(f"  Length: {len(seq)} bp")
        print(f"  Circular?: {is_circular} (similarity={score:.2f})")
        print("-" * 40)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simple Circular DNA Predictor")
    parser.add_argument("fasta", help="Input FASTA file")
    args = parser.parse_args()
    main(args.fasta)
