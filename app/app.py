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
    'TTT': 0.58, 'TTC': 1.00, 'TTA': 0.02, 'TTG': 0.07, 'CTT': 0.13, 'CTC': 0.20,
    'CTA': 0.02, 'CTG': 1.00, 'ATT': 0.49, 'ATC': 1.00, 'ATA': 0.03, 'ATG': 1.00,
    'GTT': 0.35, 'GTC': 0.47, 'GTA': 0.07, 'GTG': 1.00, 'TCT': 0.22, 'TCC': 0.39,
    'TCA': 0.15, 'TCG': 0.06, 'AGT': 0.15, 'AGC': 1.00, 'CCT': 0.42, 'CCC': 0.29,
    'CCA': 0.28, 'CCG': 1.00, 'ACT': 0.23, 'ACC': 1.00, 'ACA': 0.36, 'ACG': 0.47,
    'GCT': 0.37, 'GCC': 1.00, 'GCA': 0.28, 'GCG': 0.76, 'TAT': 0.43, 'TAC': 1.00,
    'CAT': 0.43, 'CAC': 1.00, 'CAA': 0.27, 'CAG': 1.00, 'AAT': 0.47, 'AAC': 1.00,
    'AAA': 0.44, 'AAG': 1.00, 'GAT': 0.63, 'GAC': 1.00, 'GAA': 0.68, 'GAG': 1.00,
    'TGT': 0.44, 'TGC': 1.00, 'TGG': 1.00, 'CGT': 0.36, 'CGC': 1.00, 'CGA': 0.07,
    'CGG': 0.11, 'AGA': 0.02, 'AGG': 0.02, 'GGT': 0.41, 'GGC': 1.00, 'GGA': 0.25,
    'GGG': 0.50
}
ECOLI_TAI = {
    'TTT': 0.43, 'TTC': 1.00, 'TTA': 0.17, 'TTG': 0.32, 'CTT': 0.22, 'CTC': 0.38,
    'CTA': 0.07, 'CTG': 1.00, 'ATT': 0.38, 'ATC': 0.69, 'ATA': 0.10, 'ATG': 1.00,
    'GTT': 0.31, 'GTC': 0.44, 'GTA': 0.09, 'GTG': 1.00, 'TCT': 0.32, 'TCC': 0.51,
    'TCA': 0.27, 'TCG': 0.23, 'AGT': 0.25, 'AGC': 0.55, 'CCT': 0.31, 'CCC': 0.29,
    'CCA': 0.25, 'CCG': 1.00, 'ACT': 0.28, 'ACC': 1.00, 'ACA': 0.38, 'ACG': 0.47,
    'GCT': 0.37, 'GCC': 1.00, 'GCA': 0.29, 'GCG': 0.69, 'TAT': 0.37, 'TAC': 1.00,
    'CAT': 0.41, 'CAC': 1.00, 'CAA': 0.36, 'CAG': 1.00, 'AAT': 0.48, 'AAC': 1.00,
    'AAA': 0.38, 'AAG': 1.00, 'GAT': 0.54, 'GAC': 1.00, 'GAA': 0.59, 'GAG': 1.00,
    'TGT': 0.45, 'TGC': 1.00, 'TGG': 1.00, 'CGT': 0.27, 'CGC': 1.00, 'CGA': 0.09,
    'CGG': 0.13, 'AGA': 0.05, 'AGG': 0.05, 'GGT': 0.39, 'GGC': 1.00, 'GGA': 0.21,
    'GGG': 0.47
}

# S. cerevisiae (embedded weights you asked for)
SCER_CAI = {
    'TTT':0.31,'TTC':1.00,'TTA':0.06,'TTG':0.10,'CTT':0.18,'CTC':0.28,'CTA':0.05,'CTG':1.00,
    'ATT':0.45,'ATC':1.00,'ATA':0.10,'ATG':1.00,'GTT':0.36,'GTC':0.63,'GTA':0.11,'GTG':1.00,
    'TCT':0.35,'TCC':0.73,'TCA':0.28,'TCG':0.09,'AGT':0.22,'AGC':1.00,'CCT':0.36,'CCC':0.32,
    'CCA':0.31,'CCG':1.00,'ACT':0.44,'ACC':1.00,'ACA':0.52,'ACG':0.22,'GCT':0.42,'GCC':1.00,
    'GCA':0.51,'GCG':0.19,'TAT':0.42,'TAC':1.00,'CAT':0.46,'CAC':1.00,'CAA':0.24,'CAG':1.00,
    'AAT':0.52,'AAC':1.00,'AAA':0.53,'AAG':1.00,'GAT':0.60,'GAC':1.00,'GAA':0.67,'GAG':1.00,
    'TGT':0.59,'TGC':1.00,'TGG':1.00,'CGT':0.23,'CGC':1.00,'CGA':0.08,'CGG':0.14,'AGA':0.07,
    'AGG':0.08,'GGT':0.53,'GGC':1.00,'GGA':0.34,'GGG':0.58
}
SCER_TAI = {
    'TTT':0.36,'TTC':1.00,'TTA':0.12,'TTG':0.23,'CTT':0.24,'CTC':0.43,'CTA':0.10,'CTG':1.00,
    'ATT':0.41,'ATC':0.78,'ATA':0.14,'ATG':1.00,'GTT':0.35,'GTC':0.58,'GTA':0.13,'GTG':1.00,
    'TCT':0.40,'TCC':0.69,'TCA':0.30,'TCG':0.16,'AGT':0.27,'AGC':1.00,'CCT':0.35,'CCC':0.31,
    'CCA':0.29,'CCG':1.00,'ACT':0.38,'ACC':1.00,'ACA':0.47,'ACG':0.22,'GCT':0.39,'GCC':1.00,
    'GCA':0.46,'GCG':0.23,'TAT':0.41,'TAC':1.00,'CAT':0.46,'CAC':1.00,'CAA':0.30,'CAG':1.00,
    'AAT':0.49,'AAC':1.00,'AAA':0.44,'AAG':1.00,'GAT':0.56,'GAC':1.00,'GAA':0.61,'GAG':1.00,
    'TGT':0.51,'TGC':1.00,'TGG':1.00,'CGT':0.25,'CGC':1.00,'CGA':0.12,'CGG':0.17,'AGA':0.10,
    'AGG':0.12,'GGT':0.49,'GGC':1.00,'GGA':0.32,'GGG':0.56
}

# ---------------- Zhou Table S5 per-codon High/Low ratios (first number per codon) ----------------
ECO_ZHOUS5_RATIO = {
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
SCER_ZHOUS5_RATIO = {
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
        return {c: SCER_ZHOUS5_RATIO.get(c, 1.0) for c in SENSE_CODONS}
    return {c: ECO_ZHOUS5_RATIO.get(c, 1.0) for c in SENSE_CODONS}

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

# ---------------- Zhou Copt (ratio/log2) + Binary % ----------------
def zhou_copt_ratio_for_gene(seq: str, codon_ratio: Dict[str, float]) -> Tuple[float, float]:
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
    ui.p("This app computes CAI, tAI, ENC, and Copt metrics for coding sequences."),
    ui.h4("Copt (Zhou) used here"),
    ui.tags.ul(
        ui.tags.li("Per-codon values are High/Low usage ratios from Zhou et al., Table S5 (first number per codon)."),
        ui.tags.li("Copt_ratio (Zhou) = geometric mean of per-codon ratios; Copt_log2 (Zhou) = average log2 ratio."),
        ui.tags.li(f"Copt (%) = % of codons with ratio ≥ {BINARY_THRESHOLD:.2f} (binary optimal/non-optimal)."),
        ui.tags.li("Single-codon AAs (ATG, TGG) are neutral (ratio = 1.0).")
    ),
    easy_close=True,
    footer=ui.input_action_button("help_close", "Close"),
    size="l"
)

app_ui = ui.page_fillable(
    ui.navset_bar(
        ui.nav_panel(
            "App",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_action_button("help_open", "Help / What’s Copt?"),
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
                    ui.h4("Species for CAI/tAI (FASTA/Paste modes)"),
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
        ui.nav_control(ui.a("⭐  CAI • tAI • ENC • Copt", href="#")),
        title="Codon Metrics (CAI, tAI, ENC, Copt — Zhou)"
    )
)

# ---------------- Common helpers for metrics rows ----------------
def active_sets_for_species(species_key: str) -> Tuple[Dict[str,float], Dict[str,float], Dict[str,float]]:
    if species_key == "scer":
        return SCER_CAI, SCER_TAI, copt_ratio_table_for_species("scer")
    return ECOLI_CAI, ECOLI_TAI, copt_ratio_table_for_species("ecoli")

def compute_all_metrics(seq: str, species_key: str) -> Dict[str, float]:
    cai_w, tai_w, copt_ratio = active_sets_for_species(species_key)
    ratio_gm, log2_mean = zhou_copt_ratio_for_gene(seq, copt_ratio)
    copt_percent = copt_percent_for_gene(seq, copt_ratio, BINARY_THRESHOLD)
    return {
        "CAI": round(calculate_cai(seq, cai_w), 4),
        "tAI": round(calculate_tai(seq, tai_w), 4),
        "ENC": round(calculate_enc(seq), 2),
        "Copt_ratio (Zhou)": round(ratio_gm, 4),
        "Copt_log2 (Zhou)": round(log2_mean, 4),
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
                    "ID","CAI","tAI","ENC","Copt_ratio (Zhou)","Copt_log2 (Zhou)","Copt (%)"
                ])

            rows = []
            for sid, s in seqs:
                m = compute_all_metrics(s, sp)
                rows.append({"ID": sid, **m})
            return pd.DataFrame(rows)[
                ["ID","CAI","tAI","ENC","Copt_ratio (Zhou)","Copt_log2 (Zhou)","Copt (%)"]
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
                "ID","CAI","tAI","ENC","Copt_ratio (Zhou)","Copt_log2 (Zhou)","Copt (%)"
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
                             "ENC": float("nan"), "Copt_ratio (Zhou)": float("nan"),
                             "Copt_log2 (Zhou)": float("nan"), "Copt (%)": float("nan")})
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
            ["ID","CAI","tAI","ENC","Copt_ratio (Zhou)","Copt_log2 (Zhou)","Copt (%)"]
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
        cols = ["ID","CAI","tAI","ENC","Copt_ratio (Zhou)","Copt_log2 (Zhou)","Copt (%)"]
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
        metrics = ["CAI","tAI","ENC","Copt_ratio (Zhou)","Copt_log2 (Zhou)","Copt (%)"]
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

    @render.plot(alt="Boxplots of CAI, tAI, ENC, and Copt (Zhou)")
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