# app.py
# Run: shiny run --reload app.py
# Requires: pip install shiny biopython pandas matplotlib seaborn


from __future__ import annotations

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo
from Bio import SeqIO
from collections import defaultdict, Counter
import math
import pandas as pd
import io
import matplotlib.pyplot as plt
import seaborn as sns
import re
from typing import Dict, List, Tuple

# ---------------- Built-in weights ----------------
ECOLI_CAI = {
    'TTT': 0.58, 'TTC': 1.00, 'TTA': 0.02, 'TTG': 0.07,
    'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.02, 'CTG': 1.00,
    'ATT': 0.49, 'ATC': 1.00, 'ATA': 0.03, 'ATG': 1.00,
    'GTT': 0.35, 'GTC': 0.47, 'GTA': 0.07, 'GTG': 1.00,
    'TCT': 0.22, 'TCC': 0.39, 'TCA': 0.15, 'TCG': 0.06,
    'AGT': 0.15, 'AGC': 1.00,
    'CCT': 0.42, 'CCC': 0.29, 'CCA': 0.28, 'CCG': 1.00,
    'ACT': 0.23, 'ACC': 1.00, 'ACA': 0.36, 'ACG': 0.47,
    'GCT': 0.37, 'GCC': 1.00, 'GCA': 0.28, 'GCG': 0.76,
    'TAT': 0.43, 'TAC': 1.00, 'CAT': 0.43, 'CAC': 1.00,
    'CAA': 0.27, 'CAG': 1.00, 'AAT': 0.47, 'AAC': 1.00,
    'AAA': 0.44, 'AAG': 1.00, 'GAT': 0.63, 'GAC': 1.00,
    'GAA': 0.68, 'GAG': 1.00, 'TGT': 0.44, 'TGC': 1.00,
    'TGG': 1.00,
    'CGT': 0.36, 'CGC': 1.00, 'CGA': 0.07, 'CGG': 0.11,
    'AGA': 0.02, 'AGG': 0.02,
    'GGT': 0.41, 'GGC': 1.00, 'GGA': 0.25, 'GGG': 0.50
}
ECOLI_TAI = {
    'TTT': 0.43, 'TTC': 1.00, 'TTA': 0.17, 'TTG': 0.32,
    'CTT': 0.22, 'CTC': 0.38, 'CTA': 0.07, 'CTG': 1.00,
    'ATT': 0.38, 'ATC': 0.69, 'ATA': 0.10, 'ATG': 1.00,
    'GTT': 0.31, 'GTC': 0.44, 'GTA': 0.09, 'GTG': 1.00,
    'TCT': 0.32, 'TCC': 0.51, 'TCA': 0.27, 'TCG': 0.23,
    'AGT': 0.25, 'AGC': 0.55,
    'CCT': 0.31, 'CCC': 0.29, 'CCA': 0.25, 'CCG': 1.00,
    'ACT': 0.28, 'ACC': 1.00, 'ACA': 0.38, 'ACG': 0.47,
    'GCT': 0.37, 'GCC': 1.00, 'GCA': 0.29, 'GCG': 0.69,
    'TAT': 0.37, 'TAC': 1.00, 'CAT': 0.41, 'CAC': 1.00,
    'CAA': 0.36, 'CAG': 1.00, 'AAT': 0.48, 'AAC': 1.00,
    'AAA': 0.38, 'AAG': 1.00, 'GAT': 0.54, 'GAC': 1.00,
    'GAA': 0.59, 'GAG': 1.00, 'TGT': 0.45, 'TGC': 1.00,
    'TGG': 1.00,
    'CGT': 0.27, 'CGC': 1.00, 'CGA': 0.09, 'CGG': 0.13,
    'AGA': 0.05, 'AGG': 0.05,
    'GGT': 0.39, 'GGC': 1.00, 'GGA': 0.21, 'GGG': 0.47
}

# Codon table for ENC
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
SENSE_CODONS: List[str] = sorted({c for codons in CODON_TABLE.values() for c in codons})

# ---------------- Calculators ----------------
def normalize_seq(seq: str) -> str:
    s = seq.upper().replace(" ", "").replace("\n", "").replace("U", "T")
    return re.sub(r"[^ACGT]", "", s)

def calculate_index(seq: str, weights: Dict[str, float]) -> float:
    seq = normalize_seq(seq)
    log_sum, count = 0.0, 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        w = weights.get(codon, 0.0)
        if w > 0.0:
            log_sum += math.log(w)
            count += 1
    return math.exp(log_sum / count) if count else 0.0

def calculate_cai(seq: str, cai_weights: Dict[str, float]) -> float:
    return calculate_index(seq, cai_weights)

def calculate_tai(seq: str, tai_weights: Dict[str, float]) -> float:
    return calculate_index(seq, tai_weights)

def calculate_enc(seq: str) -> float:
    seq = normalize_seq(seq)
    aa_codon_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        for aa, codons in CODON_TABLE.items():
            if codon in codons:
                aa_codon_counts[aa][codon] += 1

    Fk_list: List[Tuple[int, float]] = []
    for codons in aa_codon_counts.values():
        k = len(codons)
        if k <= 1:
            continue
        n = sum(codons.values())
        fk = sum((c / n) ** 2 for c in codons.values())
        if n > 1:
            Fk = (n * fk - 1) / (n - 1)
            if Fk > 0:
                Fk_list.append((k, Fk))

    if not Fk_list:
        return 61.0
    try:
        return float(2 + sum(k for k, _ in Fk_list) / sum(1 / Fk for _, Fk in Fk_list))
    except ZeroDivisionError:
        return 61.0

def compute_codon_metrics_from_fasta(fasta_path: str, cai_w: Dict[str, float], tai_w: Dict[str, float]) -> pd.DataFrame:
    results = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        results.append({
            "ID": record.id,
            "CAI": round(calculate_cai(seq, cai_w), 4),
            "tAI": round(calculate_tai(seq, tai_w), 4),
            "ENC": round(calculate_enc(seq), 2),
        })
    return pd.DataFrame(results)

def parse_sequences_from_text(text: str) -> List[Tuple[str, str]]:
    """Accepts FASTA text or newline-delimited raw sequences. Returns [(ID, seq)]."""
    text = text.strip()
    if not text:
        return []
    if text.startswith(">"):
        entries = []
        sid, buf = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if sid is not None and buf:
                    entries.append((sid, "".join(buf)))
                sid = line[1:].strip() or f"seq{len(entries)+1}"
                buf = []
            else:
                buf.append(line.strip())
        if sid is not None:
            entries.append((sid, "".join(buf)))
        return [(sid, normalize_seq(seq)) for sid, seq in entries if normalize_seq(seq)]
    else:
        seqs = [normalize_seq(l) for l in text.splitlines() if l.strip()]
        return [(f"seq{i+1}", s) for i, s in enumerate(seqs) if s]

def compute_codon_metrics_from_text(text: str, cai_w: Dict[str, float], tai_w: Dict[str, float]) -> pd.DataFrame:
    seqs = parse_sequences_from_text(text)
    results = []
    for sid, seq in seqs:
        results.append({
            "ID": sid,
            "CAI": round(calculate_cai(seq, cai_w), 4),
            "tAI": round(calculate_tai(seq, tai_w), 4),
            "ENC": round(calculate_enc(seq), 2),
        })
    return pd.DataFrame(results)

def read_weights_csv(fileinfo: List[FileInfo] | None) -> Dict[str, float] | None:
    if not fileinfo:
        return None
    path = fileinfo[0]["datapath"]
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if "codon" not in cols or "weight" not in cols:
        raise ValueError("Weights CSV must have columns: codon, weight")
    df = df.rename(columns={cols["codon"]: "codon", cols["weight"]: "weight"})
    df["codon"] = df["codon"].str.upper().str.replace("U", "T")
    sub = df[df["codon"].isin(SENSE_CODONS)][["codon", "weight"]].dropna()
    return dict(zip(sub["codon"], sub["weight"]))

def pick_weights(species: str, cai_up: Dict[str, float] | None, tai_up: Dict[str, float] | None) -> Tuple[Dict[str, float], Dict[str, float]]:
    if species == "ecoli":
        return ECOLI_CAI, ECOLI_TAI
    if cai_up is None or tai_up is None:
        eq = {c: 1.0 for c in SENSE_CODONS}
        return eq, eq
    return cai_up, tai_up

def make_boxplot_figure(df: pd.DataFrame, combined: bool = False):
    if df.empty or ("Error" in df.columns):
        fig = plt.figure(figsize=(4, 2))
        plt.text(0.5, 0.5, "Provide valid sequences to see results", ha="center", va="center")
        plt.axis("off")
        return fig
    if combined:
        df_m = df.melt(id_vars=["ID"], value_vars=["CAI", "tAI", "ENC"],
                       var_name="Metric", value_name="Value")
        plt.figure(figsize=(6, 4))
        sns.boxplot(x="Metric", y="Value", data=df_m)
        plt.ylabel("Value")
        plt.title("Codon Usage Metrics")
        plt.tight_layout()
        return plt.gcf()
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    sns.boxplot(y=df["CAI"], ax=axes[0], width=0.3); axes[0].set_title("CAI"); axes[0].set_ylabel("CAI"); axes[0].set_xticks([])
    sns.boxplot(y=df["tAI"], ax=axes[1], width=0.3); axes[1].set_title("tAI"); axes[1].set_ylabel("tAI"); axes[1].set_xticks([])
    sns.boxplot(y=df["ENC"], ax=axes[2], width=0.3); axes[2].set_title("ENC"); axes[2].set_ylabel("ENC"); axes[2].set_xticks([])
    plt.tight_layout()
    return fig

# ---------------- UI ----------------
app_ui = ui.page_fillable(
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_radio_buttons(
                "input_mode", "Input method",
                choices={"upload": "Upload FASTA", "paste": "Paste sequences"},
                selected="upload"
            ),
            ui.panel_conditional(
                "input.input_mode == 'upload'",
                ui.input_file("fasta", "Upload FASTA", accept=[".fa", ".fasta", ".fna"])
            ),
            ui.panel_conditional(
                "input.input_mode == 'paste'",
                ui.input_text_area(
                    "seq_text", "Paste sequences",
                    placeholder="FASTA format preferred (e.g., >id1\nATGGCC...)\nOr one DNA/RNA sequence per line.",
                    rows=10
                )
            ),
            ui.hr(),
            ui.input_radio_buttons(
                "species", "Weights",
                choices={"ecoli": "E. coli K-12 (built-in)", "custom": "Custom (upload CSVs)"},
                selected="ecoli"
            ),
            ui.panel_conditional(
                "input.species == 'custom'",
                ui.help_text("CSV must have columns: codon, weight; codons like TTT (DNA)."),
                ui.input_file("cai_csv", "Upload CAI weights (CSV)", accept=[".csv"]),
                ui.input_file("tai_csv", "Upload tAI weights (CSV)", accept=[".csv"]),
                ui.help_text("If either file is missing, equal weights (1.0) will be used.")
            ),
            ui.input_switch("combined", "Combine plots into one", value=False),
            ui.hr(),
            ui.download_button("download_csv", "Download metrics (CSV)"),
            ui.download_button("download_plot", "Download plot (PNG)"),
            width=320
        ),
        ui.layout_columns(
            ui.card(
                ui.card_header("Codon Metrics Table"),
                ui.output_data_frame("metrics_df")
            ),
            ui.card(
                ui.card_header("Copy-paste table (TSV)"),
                ui.output_ui("metrics_tsv_ui")
            ),
        ),
        ui.card(
            ui.card_header("CAI / tAI / ENC Boxplots"),
            ui.output_plot("metrics_plot", height="420px")
        ),
    ),
    title="Codon Metrics (CAI, tAI, ENC)"
)

# ---------------- Server ----------------
def server(input: Inputs, output: Outputs, session: Session):

    # Load custom weight CSVs (reactive)
    @reactive.calc
    def custom_cai():
        try:
            return read_weights_csv(input.cai_csv())
        except Exception as e:
            return {"__error__": str(e)}

    @reactive.calc
    def custom_tai():
        try:
            return read_weights_csv(input.tai_csv())
        except Exception as e:
            return {"__error__": str(e)}

    @reactive.calc
    def active_weights() -> Tuple[Dict[str, float], Dict[str, float]]:
        species = input.species()
        cai_w = custom_cai()
        tai_w = custom_tai()
        if isinstance(cai_w, dict) and "__error__" in cai_w:
            cai_w = None
        if isinstance(tai_w, dict) and "__error__" in tai_w:
            tai_w = None
        return pick_weights(species, cai_w, tai_w)

    # Central metrics DF (reactive)
    @reactive.calc
    def metrics_df_calc() -> pd.DataFrame:
        cai_w, tai_w = active_weights()
        mode = input.input_mode()
        try:
            if mode == "upload":
                files = input.fasta()
                if not files:
                    return pd.DataFrame(columns=["ID", "CAI", "tAI", "ENC"])
                path = files[0]["datapath"]
                df = compute_codon_metrics_from_fasta(path, cai_w, tai_w)
            else:
                text = input.seq_text() or ""
                df = compute_codon_metrics_from_text(text, cai_w, tai_w)
            if not df.empty:
                df = df[["ID", "CAI", "tAI", "ENC"]]
            return df
        except Exception as e:
            # Show the error in-table instead of spinning forever
            return pd.DataFrame({"Error": [str(e)]})

    # Interactive table (id matches function name)
    @render.data_frame
    def metrics_df():
        df = metrics_df_calc()
        return render.DataTable(df)

    # Copy-paste TSV UI (readonly textarea + copy button)
    @render.ui
    def metrics_tsv_ui():
        df = metrics_df_calc()
        if df.empty:
            tsv = "ID\tCAI\ttAI\tENC"
        elif "Error" in df.columns:
            tsv = df.to_csv(sep="\t", index=False, lineterminator="\n")
        else:
            tsv = df.to_csv(sep="\t", index=False, lineterminator="\n")
        return ui.TagList(
            ui.input_action_button("copy_btn", "Copy to clipboard"),
            ui.tags.textarea(
                tsv,
                id="tsv_area",
                readonly=True,
                style=(
                    "width:100%; height:240px; "
                    "font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "
                    "'Liberation Mono', 'Courier New', monospace;"
                ),
            ),
            ui.tags.script(
                """
                (function(){
                  const btn = document.getElementById("copy_btn");
                  const ta  = document.getElementById("tsv_area");
                  if (!btn || !ta) return;
                  btn.onclick = async function(){
                    ta.focus();
                    ta.select();
                    ta.setSelectionRange(0, ta.value.length);
                    try { await navigator.clipboard.writeText(ta.value); }
                    catch (e) { document.execCommand("copy"); }
                  };
                })();
                """
            )
        )

    # Plot
    @render.plot(alt="Boxplots of CAI, tAI, and ENC")
    def metrics_plot():
        df = metrics_df_calc()
        return make_boxplot_figure(df, combined=bool(input.combined()))

    # Downloads
    @render.download(filename="codon_metrics.csv")
    def download_csv():
        df = metrics_df_calc()
        with io.StringIO() as s:
            df.to_csv(s, index=False)
            yield s.getvalue()

    @render.download(filename="codon_metrics_boxplot.png")
    def download_plot():
        df = metrics_df_calc()
        fig = make_boxplot_figure(df, combined=bool(input.combined()))
        with io.BytesIO() as buf:
            fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
            plt.close(fig)
            yield buf.getvalue()

app = App(app_ui, server)