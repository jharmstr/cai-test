from __future__ import annotations

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo
from Bio import SeqIO
from collections import defaultdict, Counter
import math, pandas as pd, io, matplotlib.pyplot as plt, seaborn as sns, re
from typing import Dict, List, Tuple, Optional, Callable

# ---------------- Codon table / utilities ----------------
CODON_TABLE: Dict[str, List[str]] = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG']
}
SENSE_CODONS = sorted({c for v in CODON_TABLE.values() for c in v})
AA_FOR = {c: aa for aa, codons in CODON_TABLE.items() for c in codons}

def normalize_seq(s: str) -> str:
    s = (s or "").upper().replace(" ", "").replace("\n", "").replace("U", "T")
    s = re.sub(r"[^ACGT]", "", s)
    n = len(s) - (len(s) % 3)
    return s[:n]

# ---------------- Built-in CAI/tAI weights ----------------
ECOLI_CAI = {
    "AAA": 1.0,
  "AAC": 1.0,
  "AAG": 0.364706,
  "AAT": 0.896,
  "ACA": 0.282051,
  "ACC": 1.0,
  "ACG": 0.504274,
  "ACT": 0.350427,
  "AGA": 0.052632,
  "AGC": 1.0,
  "AGG": 0.06015,
  "AGT": 0.435294,
  "ATA": 0.121795,
  "ATC": 0.596154,
  "ATG": 1.0,
  "ATT": 1.0,
  "CAA": 0.43662,
  "CAC": 0.82716,
  "CAG": 1.0,
  "CAT": 1.0,
  "CCA": 0.248175,
  "CCC": 0.240876,
  "CCG": 1.0,
  "CCT": 0.313869,
  "CGA": 0.165414,
  "CGC": 1.0,
  "CGG": 0.157895,
  "CGT": 0.81203,
  "CTA": 0.1125,
  "CTC": 0.225,
  "CTG": 1.0,
  "CTT": 0.254167,
  "GAA": 1.0,
  "GAC": 0.541237,
  "GAG": 0.419643,
  "GAT": 1.0,
  "GCA": 0.548223,
  "GCC": 0.822335,
  "GCG": 1.0,
  "GCT": 0.279188,
  "GGA": 0.274854,
  "GGC": 1.0,
  "GGG": 0.25731,
  "GGT": 0.637427,
  "GTA": 0.437037,
  "GTC": 0.444444,
  "GTG": 1.0,
  "GTT": 0.637037,
  "TAC": 0.872093,
  "TAT": 1.0,
  "TCA": 0.470588,
  "TCC": 0.329412,
  "TCG": 0.482353,
  "TCT": 0.341176,
  "TGC": 1.0,
  "TGG": 1.0,
  "TGT": 0.731707,
  "TTA": 0.325,
  "TTC": 0.762376,
  "TTG": 0.254167,
  "TTT": 1.0
}
ECOLI_TAI = {
    'TTT': 0.188, 'TTC': 0.222, 'TTA': 0.111, 'TTG': 0.111, 'CTT': 0.188, 'CTC': 0.111,
    'CTA': 0.111, 'CTG': 0.111, 'ATT': 0.188, 'ATC': 0.333, 'ATA': 0.097, 'ATG': 1.00,
    'GTT': 0.188, 'GTC': 0.222, 'GTA': 0.556, 'GTG': 0.188, 'TCT': 0.188, 'TCC': 0.222,
    'TCA': 0.111, 'TCG': 0.111, 'AGT': 0.188, 'AGC': 0.111, 'CCT': 0.188, 'CCC': 0.111,
    'CCA': 0.111, 'CCG': 0.111, 'ACT': 0.188, 'ACC': 0.222, 'ACA': 0.222, 'ACG': 0.111,
    'GCT': 0.188, 'GCC': 0.222, 'GCA': 0.333, 'GCG': 0.188, 'TAT': 0.188, 'TAC': 0.333,
    'CAT': 0.188, 'CAC': 0.111, 'CAA': 0.222, 'CAG': 0.222, 'AAT': 0.188, 'AAC': 0.444,
    'AAA': 0.667, 'AAG': 0.188, 'GAT': 0.188, 'GAC': 0.333, 'GAA': 0.444, 'GAG': 0.188,
    'TGT': 0.188, 'TGC': 0.111, 'TGG': 0.111, 'CGT': 0.333, 'CGC': 0.360, 'CGA': 0.063,
    'CGG': 0.111, 'AGA': 0.222, 'AGG': 0.111, 'GGT': 0.188, 'GGC': 0.444, 'GGA': 0.222,
    'GGG': 0.111
}

# S. cerevisiae (default weights)
SCER_CAI = {
    "AAA": 1.0,
  "AAC": 0.695763,
  "AAG": 0.73592,
  "AAT": 1.0,
  "ACA": 0.87596,
  "ACC": 0.627873,
  "ACG": 0.392727,
  "ACT": 1.0,
  "AGA": 1.0,
  "AGC": 0.414999,
  "AGG": 0.433481,
  "AGT": 0.602161,
  "ATA": 0.590443,
  "ATC": 0.569731,
  "ATG": 1.0,
  "ATT": 1.0,
  "CAA": 1.0,
  "CAC": 0.570573,
  "CAG": 0.443874,
  "CAT": 1.0,
  "CCA": 1.0,
  "CCC": 0.37035,
  "CCG": 0.289173,
  "CCT": 0.737732,
  "CGA": 0.140652,
  "CGC": 0.122181,
  "CGG": 0.081614,
  "CGT": 0.30048,
  "CTA": 0.493425,
  "CTC": 0.200171,
  "CTG": 0.385723,
  "CTT": 0.450947,
  "GAA": 1.0,
  "GAC": 0.537565,
  "GAG": 0.421948,
  "GAT": 1.0,
  "GCA": 0.765478,
  "GCC": 0.595246,
  "GCG": 0.291693,
  "GCT": 1.0,
  "GGA": 0.456194,
  "GGC": 0.409349,
  "GGG": 0.252125,
  "GGT": 1.0,
  "GTA": 0.533315,
  "GTC": 0.533454,
  "GTG": 0.487629,
  "GTT": 1.0,
  "TAC": 0.787074,
  "TAT": 1.0,
  "TCA": 0.794676,
  "TCC": 0.605137,
  "TCG": 0.364366,
  "TCT": 1.0,
  "TGC": 0.587774,
  "TGG": 1.0,
  "TGT": 1.0,
  "TTA": 0.962331,
  "TTC": 0.706116,
  "TTG": 1.0,
  "TTT": 1.0
}
SCER_TAI = {
    'TTT':0.069,'TTC':0.518,'TTA':0.362,'TTG':0.793,'CTT':0.007,'CTC':0.052,'CTA':0.155,'CTG':0.118,
    'ATT':0.673,'ATC':0.219,'ATA':0.104,'ATG':0.518,'GTT':0.725,'GTC':0.219,'GTA':0.104,'GTG':0.182,
    'TCT':0.569,'TCC':0.219,'TCA':0.155,'TCG':0.170,'AGT':0.014,'AGC':0.014,'CCT':0.104,'CCC':0.219,
    'CCA':0.518,'CCG':0.393,'ACT':0.569,'ACC':0.219,'ACA':0.207,'ACG':0.209,'GCT':0.569,'GCC':0.219,
    'GCA':0.259,'GCG':0.197,'TAT':0.055,'TAC':0.414,'CAT':0.048,'CAC':0.362,'CAA':0.466,'CAG':0.406,
    'AAT':0.069,'AAC':0.518,'AAA':0.362,'AAG':1.00,'GAT':0.104,'GAC':0.776,'GAA':0.725,'GAG':0.654,
    'TGT':0.028,'TGC':0.207,'TGG':0.311,'CGT':0.311,'CGC':0.219,'CGA':0.219,'CGG':0.052,'AGA':0.569,
    'AGG':0.485,'GGT':0.111,'GGC':0.828,'GGA':0.155,'GGG':0.222
}

#COPT ratio values

ECO_COPT_RATIO = {
    'GCT': 1.717, 'GCC': 0.639, 'GCA': 0.986, 'GCG': 0.900,
    'CGT': 2.934, 'CGC': 1.189, 'CGA': 0.276, 'CGG': 0.210, 'AGA': 0.144, 'AGG': 0.153,
    'AAT': 0.271, 'AAC': 3.693,
    'GAT': 0.455, 'GAC': 2.198,
    'TGT': 0.622, 'TGC': 1.608,
    'CAA': 0.507, 'CAG': 1.972,
    'GAA': 1.331, 'GAG': 0.751,
    'GGT': 1.840, 'GGC': 1.474, 'GGA': 0.234, 'GGG': 0.327,
    'CAT': 0.381, 'CAC': 2.622,
    'ATT': 0.539, 'ATC': 3.406, 'ATA': 0.130,
    'CTT': 0.542, 'CTC': 0.744, 'CTA': 0.285, 'CTG': 3.759, 'TTA': 0.301, 'TTG': 0.527,
    'AAA': 1.136, 'AAG': 0.881,
    'TTT': 0.285, 'TTC': 3.510,
    'CCT': 0.526, 'CCC': 0.278, 'CCA': 0.737, 'CCG': 2.883,
    'TCT': 2.360, 'TCC': 2.142, 'TCA': 0.418, 'TCG': 0.581, 'AGT': 0.343, 'AGC': 1.086,
    'ACT': 1.963, 'ACC': 1.886, 'ACA': 0.293, 'ACG': 0.413,
    'TAT': 0.348, 'TAC': 2.876,
    'GTT': 1.668, 'GTC': 0.579, 'GTA': 1.221, 'GTG': 0.707,
    'ATG': 1.000,  # Met
    'TGG': 1.000,  # Trp
}
SCER_COPT_RATIO = {
    'GCT': 3.919, 'GCC': 1.318, 'GCA': 0.150, 'GCG': 0.122,
    'CGT': 1.554, 'CGC': 0.157, 'CGA': 0.042, 'CGG': 0.040, 'AGA': 4.273, 'AGG': 0.122,
    'AAT': 0.193, 'AAC': 5.188,
    'GAT': 0.485, 'GAC': 2.061,
    'TGT': 3.671, 'TGC': 0.272,
    'CAA': 6.547, 'CAG': 0.153,
    'GAA': 4.168, 'GAG': 0.240,
    'GGT': 10.496, 'GGC': 0.304, 'GGA': 0.120, 'GGG': 0.111,
    'CAT': 0.281, 'CAC': 3.561,
    'ATT': 1.261, 'ATC': 2.669, 'ATA': 0.108,
    'CTT': 0.542, 'CTC': 0.744, 'CTA': 0.285, 'CTG': 3.759, 'TTA': 0.301, 'TTG': 0.527,
    'AAA': 0.208, 'AAG': 4.811,
    'TTT': 0.251, 'TTC': 3.985,
    'CCT': 0.423, 'CCC': 0.144, 'CCA': 6.153, 'CCG': 0.126,
    'TCT': 2.819, 'TCC': 2.565, 'TCA': 0.331, 'TCG': 0.191, 'AGT': 0.364, 'AGC': 0.413,
    'ACT': 2.056, 'ACC': 2.640, 'ACA': 1.018, 'ACG': 0.126,
    'TAT': 0.230, 'TAC': 4.340,
    'GTT': 1.950, 'GTC': 2.696, 'GTA': 0.135, 'GTG': 0.221,
    'ATG': 1.000,  # Met
    'TGG': 1.000,  # Trp
}

def copt_ratio_table_for_species(species_key: str) -> Dict[str, float]:
    if species_key == "scer":
        return {c: SCER_COPT_RATIO.get(c, 1.0) for c in SENSE_CODONS}
    return {c: ECO_COPT_RATIO.get(c, 1.0) for c in SENSE_CODONS}

# ---------------- CAI / tAI / ENC calculators ----------------
def calculate_index(seq: str, weights: Dict[str, float]) -> float:
    seq = normalize_seq(seq)
    log_sum, count = 0.0, 0
    for i in range(0, len(seq) - 2, 3):
        w = weights.get(seq[i:i+3], 0.0)
        if w > 0.0:
            log_sum += math.log(w); count += 1
    return math.exp(log_sum / count) if count else 0.0

def calculate_cai(seq: str, cai_w: Dict[str, float]) -> float: return calculate_index(seq, cai_w)
def calculate_tai(seq: str, tai_w: Dict[str, float]) -> float: return calculate_index(seq, tai_w)

def calculate_enc(seq: str) -> float:
    seq = normalize_seq(seq)
    aa_codon_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = AA_FOR.get(codon)
        if aa:
            aa_codon_counts[aa][codon] += 1
    Fk_list: List[Tuple[int, float]] = []
    for codons in aa_codon_counts.values():
        k = len(codons)
        if k <= 1: continue
        n = sum(codons.values())
        fk = sum((c / n) ** 2 for c in codons.values())
        if n > 1:
            Fk = (n * fk - 1) / (n - 1)
            if Fk > 0: Fk_list.append((k, Fk))
    if not Fk_list: return 61.0
    try: return float(2 + sum(k for k, _ in Fk_list) / sum(1 / Fk for _, Fk in Fk_list))
    except ZeroDivisionError: return 61.0

# ---------------- Copt (ratio/log2) + Binary % ----------------
def copt_ratio_for_gene(seq: str, codon_ratio: Dict[str, float]) -> Tuple[float, float]:
    seq = normalize_seq(seq)
    logs = []
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if c in AA_FOR:
            r = max(codon_ratio.get(c, 1.0), 1e-12)
            logs.append(math.log(r, 2))
    if not logs:
        return 1.0, 0.0
    mean_log2 = sum(logs) / len(logs)
    gm_ratio = float(2 ** mean_log2)
    return gm_ratio, float(mean_log2)

BINARY_THRESHOLD = 1.0  # codon considered "optimal" if ratio >= 1

def copt_percent_for_gene(seq: str, codon_ratio: Dict[str, float], threshold: float = BINARY_THRESHOLD) -> float:
    seq = normalize_seq(seq)
    total, optimal = 0, 0
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if c in AA_FOR:
            total += 1
            if codon_ratio.get(c, 1.0) >= threshold:
                optimal += 1
    return (optimal / total * 100.0) if total > 0 else 0.0

# ---------------- Sequence helpers ----------------
def parse_sequences_from_text(text: str) -> List[Tuple[str, str]]:
    text = (text or "").strip()
    if not text: return []
    if text.startswith(">"):
        entries, sid, buf = [], None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if sid is not None and buf: entries.append((sid, "".join(buf)))
                sid, buf = line[1:].strip() or f"seq{len(entries)+1}", []
            else:
                buf.append(line.strip())
        if sid is not None: entries.append((sid, "".join(buf)))
        return [(sid, normalize_seq(seq)) for sid, seq in entries if normalize_seq(seq)]
    else:
        seqs = [normalize_seq(l) for l in text.splitlines() if l.strip()]
        return [(f"seq{i+1}", s) for i, s in enumerate(seqs) if s]

# ---------------- UI ----------------
help_modal = ui.modal(
    ui.h3("Help & Notes"),
    ui.p("This app computes CAI, tAI, ENC, and Copt metrics for coding sequences. See readme at https://github.com/jharmstr/codon_metrics_app for more details"),
    ui.h4("CAI"),
    ui.tags.ul(
        ui.tags.li("Codon Adaptation Index from Sharp et al."),
        ui.tags.li("Weights from public data"),
    ),
    ui.h4("tAI"),
    ui.tags.ul(
        ui.tags.li("tRNA Adaptation Index from dos Reis et al."),
        ui.tags.li("Weights from public data"),
    ),
    ui.h4("ENC"),
    ui.tags.ul(
        ui.tags.li("Effective Number of Codons from Wright"),
    ),
    ui.h4("Copt"),
    ui.tags.ul(
        ui.tags.li("Per-codon values are High/Low usage ratios from Zhou et al."),
        ui.tags.li("Copt_ratio = geometric mean of per-codon ratios; Copt_log2 = average log2 ratio."),
        ui.tags.li(f"Copt (%) = % of codons with ratio ≥ {BINARY_THRESHOLD:.2f} (binary optimal/non-optimal)."),
        ui.tags.li("Single-codon AAs (ATG, TGG) are neutral (ratio = 1.0)."),
    ),
    easy_close=True,
    footer=ui.input_action_button("help_close", "Close"),
    size="l",
)

app_ui = ui.page_fillable(
    ui.navset_bar(
        ui.nav_panel(
            "App",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_action_button("help_open", "Help"),
                    ui.hr(),
                    ui.h4("Input mode"),
                    ui.input_radio_buttons(
                        "input_mode", None,
                        choices={"upload": "Upload FASTA", "paste": "Paste sequences", "table": "Upload table (CSV)"},
                        selected="upload"
                    ),
                    ui.panel_conditional("input.input_mode == 'upload'",
                        ui.input_file("fasta", "Upload FASTA", accept=[".fa",".fasta",".fna"])
                    ),
                    ui.panel_conditional("input.input_mode == 'paste'",
                        ui.input_text_area("seq_text", "Paste FASTA or one sequence per line", rows=8)
                    ),
                    ui.panel_conditional("input.input_mode == 'table'",
                        ui.TagList(
                            ui.input_file("table_file", "Upload CSV table", accept=[".csv"]),
                            ui.output_ui("table_col_picker"),
                            ui.input_select("fixed_species", "Set species for ALL rows (optional; overrides column)",
                                            choices={"": "(none)", "ecoli": "E. coli", "scer": "S. cerevisiae"},
                                            selected="")
                        )
                    ),
                    ui.hr(),
                    ui.h4("Species(for expression optimization)"),
                    ui.input_radio_buttons("species", None,
                        choices={"ecoli": "E. coli (built-in)", "scer": "S. cerevisiae (built-in)"},
                        selected="ecoli"),
                    ui.input_switch("combined", "Combine plots", value=False),
                    ui.download_button("download_csv", "Download metrics (CSV)"),
                    ui.download_button("download_plot", "Download plot (PNG)"),
                    ui.panel_conditional("input.input_mode == 'table'",
                        ui.download_button("download_augmented", "Download augmented table (CSV)")
                    ),
                    width=380
                ),
                ui.layout_columns(
                    ui.card(ui.card_header("Codon Metrics Table"), ui.output_data_frame("metrics_df")),
                    ui.card(ui.card_header("Copy-paste table (TSV)"), ui.output_ui("metrics_tsv_ui")),
                ),
                ui.card(ui.card_header("CAI / tAI / ENC / Copt"),
                        ui.output_plot("metrics_plot", height="460px"))
            ),
        ),
        ui.nav_spacer(),
        ui.nav_control(ui.a("CAI • tAI • ENC • Copt", href="#")),
        title="Codon Metrics (CAI, tAI, ENC, Copt)"
    )
)

# ---------------- Common helpers for metrics rows ----------------
def active_sets_for_species(species_key: str) -> Tuple[Dict[str,float], Dict[str,float], Dict[str,float]]:
    if species_key == "scer":
        return SCER_CAI, SCER_TAI, copt_ratio_table_for_species("scer")
    return ECOLI_CAI, ECOLI_TAI, copt_ratio_table_for_species("ecoli")

def compute_all_metrics(seq: str, species_key: str) -> Dict[str, float]:
    cai_w, tai_w, copt_ratio = active_sets_for_species(species_key)
    ratio_gm, log2_mean = copt_ratio_for_gene(seq, copt_ratio)
    copt_percent = copt_percent_for_gene(seq, copt_ratio, BINARY_THRESHOLD)
    return {
        "CAI": round(calculate_cai(seq, cai_w), 4),
        "tAI": round(calculate_tai(seq, tai_w), 4),
        "ENC": round(calculate_enc(seq), 2),
        "Copt_ratio": round(ratio_gm, 4),
        "Copt_log2": round(log2_mean, 4),
        "Copt (%)": round(copt_percent, 2),
    }

# ---------------- Server ----------------
def server(input: Inputs, output: Outputs, session: Session):
    # Help modal controls
    @reactive.effect
    @reactive.event(input.help_open)
    def _open_help(): ui.modal_show(help_modal)
    @reactive.effect
    @reactive.event(input.help_close)
    def _close_help(): ui.modal_remove()

    # ===== FASTA / Paste modes =====
    @reactive.calc
    def metrics_df_calc_fp() -> pd.DataFrame:
        sp = input.species()
        try:
            seqs: List[Tuple[str, str]] = []
            if input.input_mode() == "upload":
                files = input.fasta()
                if files:
                    for rec in SeqIO.parse(files[0]["datapath"], "fasta"):
                        seqs.append((rec.id, str(rec.seq)))
            elif input.input_mode() == "paste":
                seqs.extend(parse_sequences_from_text(input.seq_text() or ""))

            if not seqs:
                return pd.DataFrame(columns=[
                    "ID","CAI","tAI","ENC","Copt_ratio","Copt_log2","Copt (%)"
                ])

            rows = []
            for sid, s in seqs:
                m = compute_all_metrics(s, sp)
                rows.append({"ID": sid, **m})
            return pd.DataFrame(rows)[
                ["ID","CAI","tAI","ENC","Copt_ratio","Copt_log2","Copt (%)"]
            ]
        except Exception as e:
            return pd.DataFrame({"Error":[str(e)]})

    # ===== TABLE mode =====
    @reactive.calc
    def table_df_raw() -> Optional[pd.DataFrame]:
        if input.input_mode() != "table":
            return None
        f = input.table_file()
        if not f:
            return None
        try:
            return pd.read_csv(f[0]["datapath"])
        except Exception as e:
            return pd.DataFrame({"Error":[f"Failed to read CSV: {e}"]})

    @output
    @render.ui
    def table_col_picker():
        df = table_df_raw()
        if df is None or df.empty or "Error" in df.columns:
            return ui.help_text("Upload a CSV to choose columns.")
        cols = list(df.columns)
        # cheap guesses
        seq_guess = next((c for c in cols if c.lower() in ("cds","sequence","seq","cds_seq","CDS")), cols[0])
        id_guess  = next((c for c in cols if c.lower() in ("id","tx_id","gene","name","locus","accession")), cols[0])
        species_guess = next((c for c in cols if c.lower() in ("species","organism","host")), "")
        return ui.TagList(
            ui.input_select("id_col", "ID column", choices=cols, selected=id_guess),
            ui.input_select("seq_col", "Sequence column", choices=cols, selected=seq_guess),
            ui.input_select("species_col", "Species column (optional: 'ecoli' or 'scer')",
                            choices=["(none)"] + cols, selected=species_guess or "(none)")
        )

    @reactive.calc
    def metrics_df_calc_table() -> pd.DataFrame:
        df = table_df_raw()
        if df is None or df.empty:
            return pd.DataFrame(columns=[
                "ID","CAI","tAI","ENC","Copt_ratio","Copt_log2","Copt (%)"
            ])
        if "Error" in df.columns:
            return df

        id_col = input.id_col() if "id_col" in input else None
        seq_col = input.seq_col() if "seq_col" in input else None
        species_col_sel = input.species_col() if "species_col" in input else "(none)"
        fixed_species = input.fixed_species() if "fixed_species" in input else ""

        if not id_col or id_col not in df.columns:
            id_col = df.columns[0]
        if not seq_col or seq_col not in df.columns:
            # No usable sequence column
            return pd.DataFrame({"Error":[
                "Select a valid sequence column in the sidebar (Table mode)."
            ]})

        rows = []
        for _, r in df.iterrows():
            sid = str(r[id_col])
            seq = normalize_seq(str(r[seq_col]))
            if not seq:
                rows.append({"ID": sid, "CAI": float("nan"), "tAI": float("nan"),
                             "ENC": float("nan"), "Copt_ratio": float("nan"),
                             "Copt_log2": float("nan"), "Copt (%)": float("nan")})
                continue

            # Determine species for this row
            if fixed_species:
                sp = fixed_species
            else:
                if species_col_sel and species_col_sel != "(none)" and species_col_sel in df.columns:
                    val = str(r[species_col_sel]).strip().lower()
                    sp = "scer" if val.startswith("scer") or val in {"s.cer","s cerevisiae","s. cerevisiae","yeast"} else \
                         "ecoli" if val.startswith("e") or "coli" in val else "ecoli"
                else:
                    sp = "ecoli"

            m = compute_all_metrics(seq, sp)
            rows.append({"ID": sid, **m})
        return pd.DataFrame(rows)[
            ["ID","CAI","tAI","ENC","Copt_ratio","Copt_log2","Copt (%)"]
        ]

    # ===== Unified outputs (show metrics table depending on mode) =====
    def metrics_df_current() -> pd.DataFrame:
        mode = input.input_mode()
        if mode == "table":
            return metrics_df_calc_table()
        return metrics_df_calc_fp()

    # Interactive table
    @render.data_frame
    def metrics_df():
        return render.DataTable(metrics_df_current())

    # Copy-paste TSV
    @render.ui
    def metrics_tsv_ui():
        df = metrics_df_current()
        cols = ["ID","CAI","tAI","ENC","Copt_ratio","Copt_log2","Copt (%)"]
        if df.empty:
            tsv = "\t".join(cols)
        else:
            # if an Error column is present, show it plainly
            if "Error" in df.columns:
                tsv = df.to_csv(sep="\t", index=False, lineterminator="\n")
            else:
                tsv = df[cols].to_csv(sep="\t", index=False, lineterminator="\n")
        return ui.TagList(
            ui.input_action_button("copy_btn", "Copy to clipboard"),
            ui.tags.textarea(
                tsv, id="tsv_area", readonly=True,
                style=("width:100%; height:260px; font-family: ui-monospace, SFMono-Regular, "
                       "Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;")
            ),
            ui.tags.script(
                """
                (function(){
                  const btn = document.getElementById("copy_btn");
                  const ta  = document.getElementById("tsv_area");
                  if (!btn || !ta) return;
                  btn.onclick = async function(){
                    ta.focus(); ta.select(); ta.setSelectionRange(0, ta.value.length);
                    try { await navigator.clipboard.writeText(ta.value); }
                    catch (e) { document.execCommand("copy"); }
                  };
                })();
                """
            )
        )

    # Plot
    def make_boxplot_figure(df: pd.DataFrame, combined: bool):
        if df.empty or ("Error" in df.columns):
            fig = plt.figure(figsize=(4, 2))
            plt.text(0.5, 0.5, "Provide sequences", ha="center", va="center")
            plt.axis("off"); return fig
        metrics = ["CAI","tAI","ENC","Copt_ratio","Copt_log2","Copt (%)"]
        if combined:
            dfm = df.melt(id_vars=["ID"], value_vars=metrics, var_name="Metric", value_name="Value")
            plt.figure(figsize=(10, 4)); sns.boxplot(x="Metric", y="Value", data=dfm)
            plt.ylabel("Value"); plt.title("Codon Usage Metrics"); plt.tight_layout()
            return plt.gcf()
        fig, axes = plt.subplots(1, len(metrics), figsize=(4*len(metrics), 4))
        for ax, col in zip(axes, metrics):
            sns.boxplot(y=df[col], ax=ax, width=0.3)
            ax.set_title(col); ax.set_ylabel(col); ax.set_xticks([])
        plt.tight_layout(); return fig

    @render.plot(alt="Boxplots of CAI, tAI, ENC, and Copt")
    def metrics_plot():
        return make_boxplot_figure(metrics_df_current(), combined=bool(input.combined()))

    # Downloads (metrics only)
    @render.download(filename="codon_metrics.csv")
    def download_csv():
        df = metrics_df_current()
        with io.StringIO() as s:
            df.to_csv(s, index=False); yield s.getvalue()

    @render.download(filename="codon_metrics_boxplot.png")
    def download_plot():
        fig = make_boxplot_figure(metrics_df_current(), combined=bool(input.combined()))
        with io.BytesIO() as buf:
            fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
            plt.close(fig); yield buf.getvalue()

    # Download augmented (TABLE mode only): merges original CSV + computed columns by ID
    @render.download(filename=lambda: augmented_name())
    def download_augmented():
        # guard: only in table mode
        if input.input_mode() != "table":
            yield b""; return
        raw = table_df_raw()
        metrics = metrics_df_calc_table()
        if raw is None or raw.empty or metrics.empty or "Error" in metrics.columns:
            with io.StringIO() as s:
                s.write("Error: no augmented output available.\n")
                yield s.getvalue().encode("utf-8")
            return
        # Merge on ID
        # If ID column name != 'ID', we still export both (original ID col + computed 'ID').
        # Prefer left-merge to preserve row order and all original columns.
        id_col = input.id_col() if "id_col" in input else None
        if id_col and id_col in raw.columns:
            merged = raw.merge(metrics, left_on=id_col, right_on="ID", how="left")
        else:
            # no id column selected; append metrics in order (best effort)
            merged = pd.concat([raw.reset_index(drop=True), metrics.reset_index(drop=True)], axis=1)
        with io.StringIO() as s:
            merged.to_csv(s, index=False)
            yield s.getvalue().encode("utf-8")

    def augmented_name() -> str:
        f = input.table_file()
        if f:
            name = f[0]["name"]
            base = name.rsplit(".", 1)[0]
            return f"{base}_with_metrics.csv"
        return "table_with_metrics.csv"

app = App(app_ui, server)