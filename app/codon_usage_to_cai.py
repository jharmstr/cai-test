import re
import json
from collections import defaultdict

# --- Paste your codon-usage text (RNA alphabet with U) ---
USAGE_TEXT = """
UUU 26.1(170666)  UCU 23.5(153557)  UAU 18.8(122728)  UGU  8.1( 52903)
UUC 18.4(120510)  UCC 14.2( 92923)  UAC 14.8( 96596)  UGC  4.8( 31095)
UUA 26.2(170884)  UCA 18.7(122028)  UAA  1.1(  6913)  UGA  0.7(  4447)
UUG 27.2(177573)  UCG  8.6( 55951)  UAG  0.5(  3312)  UGG 10.4( 67789)

CUU 12.3( 80076)  CCU 13.5( 88263)  CAU 13.6( 89007)  CGU  6.4( 41791)
CUC  5.4( 35545)  CCC  6.8( 44309)  CAC  7.8( 50785)  CGC  2.6( 16993)
CUA 13.4( 87619)  CCA 18.3(119641)  CAA 27.3(178251)  CGA  3.0( 19562)
CUG 10.5( 68494)  CCG  5.3( 34597)  CAG 12.1( 79121)  CGG  1.7( 11351)

AUU 30.1(196893)  ACU 20.3(132522)  AAU 35.7(233124)  AGU 14.2( 92466)
AUC 17.2(112176)  ACC 12.7( 83207)  AAC 24.8(162199)  AGC  9.8( 63726)
AUA 17.8(116254)  ACA 17.8(116084)  AAA 41.9(273618)  AGA 21.3(139081)
AUG 20.9(136805)  ACG  8.0( 52045)  AAG 30.8(201361)  AGG  9.2( 60289)

GUU 22.1(144243)  GCU 21.2(138358)  GAU 37.6(245641)  GGU 23.9(156109)
GUC 11.8( 76947)  GCC 12.6( 82357)  GAC 20.2(132048)  GGC  9.8( 63903)
GUA 11.8( 76927)  GCA 16.2(105910)  GAA 45.6(297944)  GGA 10.9( 71216)
GUG 10.8( 70337)  GCG  6.2( 40358)  GAG 19.2(125717)  GGG  6.0( 39359)
""".strip()

# --- Standard genetic code (RNA codons -> amino acids). Stops marked '*'.
RNA_CODON_TO_AA = {
    # U-start
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L","UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*","UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    # C-start
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L","CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q","CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    # A-start
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M","ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K","AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    # G-start
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V","GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E","GGU":"G","GGC":"G","GGA":"G","GGG":"G"
}

def parse_usage_block(text: str):
    """
    Parse lines like 'UUU 19.7(   101)' into { 'UUU': {'count':101, 'per_thou':19.7}, ... }.
    Picks up both tokens even if spacing is irregular.
    """
    codon_re = re.compile(r"\b([ACGU]{3})\b")
    # number outside parens (per-thousand), then parentheses with an integer
    value_re = re.compile(r"([0-9.]+)\s*\(\s*([0-9]+)\s*\)")
    tokens = text.split()
    # Easier: scan the whole text with joint regex "CODON ... value(value)"
    items = {}
    for m in re.finditer(r"([ACGU]{3})\s*([0-9.]+)\s*\(\s*([0-9]+)\s*\)", text):
        codon, per_thou, count = m.group(1), float(m.group(2)), int(m.group(3))
        items[codon] = {"count": count, "per_thou": per_thou}
    # Some tables might omit parentheses for a few codons; try to catch lone CODON + number
    if not items:
        # Fallback: naive pairwise scan
        i = 0
        while i < len(tokens)-1:
            if codon_re.fullmatch(tokens[i]) and re.fullmatch(r"[0-9.]+", tokens[i+1]):
                codon = tokens[i]; per_thou = float(tokens[i+1]); count = None
                items[codon] = {"count": count, "per_thou": per_thou}
                i += 2
            else:
                i += 1
    return items

def compute_RSCU_and_CAI_weights(items: dict):
    """
    items: dict codon -> {'count': int|None, 'per_thou': float|None}
    Returns: (RSCU dict (DNA), CAI weights dict (DNA))
    """
    # choose signal = counts if present else per_thou
    signal = {}
    for codon_rna, vals in items.items():
        if vals.get("count") is not None:
            signal[codon_rna] = float(vals["count"])
        else:
            signal[codon_rna] = float(vals.get("per_thou", 0.0))

    # group by amino acid (exclude stops)
    aa_to_codons = defaultdict(list)
    for codon_rna, aa in RNA_CODON_TO_AA.items():
        if aa != "*":
            aa_to_codons[aa].append(codon_rna)

    # compute RSCU
    RSCU_rna = {}
    for aa, codons in aa_to_codons.items():
        total = sum(signal.get(c, 0.0) for c in codons)
        s = len(codons)
        denom = (total / s) if s > 0 else 1.0
        for c in codons:
            RSCU_rna[c] = (signal.get(c, 0.0) / denom) if denom > 0 else 0.0

    # CAI weights: RSCU / max_RSCU within AA
    CAI_rna = {}
    for aa, codons in aa_to_codons.items():
        max_r = max(RSCU_rna.get(c, 0.0) for c in codons) or 1.0
        for c in codons:
            CAI_rna[c] = (RSCU_rna.get(c, 0.0) / max_r)

    # convert RNA->DNA (U->T), keep stops excluded
    def rna_to_dna(c): return c.replace("U", "T")
    RSCU_dna = {rna_to_dna(c): round(v, 6) for c, v in RSCU_rna.items()}
    CAI_dna  = {rna_to_dna(c): round(v, 6) for c, v in CAI_rna.items()}

    return RSCU_dna, CAI_dna

def main():
    items = parse_usage_block(USAGE_TEXT)
    if not items:
        raise SystemExit("Failed to parse any codons. Check the formatting of USAGE_TEXT.")
    RSCU, CAI = compute_RSCU_and_CAI_weights(items)

    # Pretty print for direct paste into your app
    print("# --- CAI weights (E. coli) from the provided usage table (DNA alphabet) ---")
    # As Python dict literal
    as_pairs = ", ".join([f"'{k}': {v}" for k, v in sorted(CAI.items())])
    print(f"ECOLI_CAI_FROM_TABLE = {{ {as_pairs} }}\n")

    # Also dump JSON (if you prefer)
    print("# JSON (same content):")
    print(json.dumps(CAI, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()