import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import hypergeom

# ── Config ────────────────────────────────────────────────────────────────────

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]
NET_DIR         = "network_enrichment_results"

COMPARISONS = {
    "diabetic_empty_42-nondiabetic_empty_42": "Empty defect",
    "diabetic_PCL_42-nondiabetic_PCL_42":     "PCL scaffold",
}

labels = [
    "Tissue-level\nDEPs",
    "Tissue-level DEPs +\nconnector proteins",
    "First-level\nDEPs",
]

COLORS = {
    "Empty defect": "#4A90C4",
    "PCL scaffold":  "#E07B54",
}
LINE_COLOR = "#1A1A1A"
BAR_WIDTH  = 0.32
GROUP_GAP  = 0.08

# Load data that is shared across all threshold combinations
ai_de     = pd.read_csv("de_analysis_results_AI_RobNorm_scaled/de_results_raw.csv")
as_de     = pd.read_csv("de_analysis_results_AS_RobNorm_scaled/de_results_raw.csv")
bone_caps = set(pd.read_csv("input/bone_caps_meta_analysis.csv")["Gene"].dropna().str.upper())

# Collect results for supplement export
_results = []

# ── Loop over all threshold combinations ──────────────────────────────────────

for fc in FC_THRESHOLDS:
    for stab in STAB_THRESHOLDS:
        THRESH_DIR = f"valid_DEPs/FC{fc}_Stab{stab}"

        print(f"\n{'='*60}")
        print(f"FC={fc}  Stab={stab}  ({THRESH_DIR})")
        print(f"{'='*60}")

        # ── Compute protein counts and p-values from source files ─────────────

        # Reads:
        #   - summary_statistics.csv  → tissue-level / first-level counts and overlaps
        #   - DE results              → total tested proteome size per comparison
        #   - exception_proteins.csv  → connector proteins from both-levels network enrichment
        #   - bone_caps_meta_analysis → reference gene list for hypergeometric test

        summary = pd.read_csv(os.path.join(THRESH_DIR, "summary_statistics.csv")).set_index("Comparison")

        PROTEIN_N     = {}
        protein_pvals = {}

        for comp_key, comp_label in COMPARISONS.items():
            row = summary.loc[comp_key]

            # Total proteome size = all proteins tested in either fraction for this comparison
            ai_prots = set(ai_de[ai_de["Comparison"] == comp_key]["Protein.IDs"].dropna())
            as_prots = set(as_de[as_de["Comparison"] == comp_key]["Protein.IDs"].dropna())
            M = len(ai_prots | as_prots)       # population size
            n = int(row["Bone_Caps_Tested"])   # successes in population

            # Set 1: tissue-level DEPs (first-level + second-level combined)
            s1 = int(row["Validated_Both_Fractions"] + row["Validated_Exclusive"] + row["Validated_Second_Level"])
            o1 = int(row["Overlap_First_Level"] + row["Overlap_Second_Level"])
            # hypergeom.sf(k, M, n, N) with n the number of successes in M and N the number of draws
            p1 = hypergeom.sf(o1 - 1, M, n, s1)

            # Set 2: tissue-level DEPs + connector proteins (from both-levels network)
            exc_path  = os.path.join(NET_DIR, comp_key, "both_levels", "exception_proteins.csv")
            exc_prots = set(pd.read_csv(exc_path)["Gene"].dropna().str.upper())
            s2 = s1 + len(exc_prots)
            o2 = o1 + len(exc_prots & bone_caps)
            p2 = hypergeom.sf(o2 - 1, M, n, s2)

            # Set 3: first-level DEPs only
            s3 = int(row["Validated_Both_Fractions"] + row["Validated_Exclusive"])
            o3 = int(row["Overlap_First_Level"])
            p3 = hypergeom.sf(o3 - 1, M, n, s3)

            PROTEIN_N[comp_label]     = [s1, s2, s3]
            protein_pvals[comp_label] = [p1, p2, p3]

            expected1 = round((s1 * n) / M, 2) if M > 0 else 0
            expected2 = round((s2 * n) / M, 2) if M > 0 else 0
            expected3 = round((s3 * n) / M, 2) if M > 0 else 0

            for set_label, s, o, expected, p in [
                ("Tissue-level DEPs",                    s1, o1, expected1, p1),
                ("Tissue-level DEPs + connector proteins", s2, o2, expected2, p2),
                ("First-level DEPs",                     s3, o3, expected3, p3),
            ]:
                _results.append({
                    "FC threshold":             fc,
                    "Stability threshold":      stab,
                    "Comparison":               comp_label,
                    "Protein set":              set_label,
                    "Proteome size (M)":        M,
                    "Bone CAPs in proteome (n)": n,
                    "Protein set size (N)":     s,
                    "Bone CAPs overlap (k)":    o,
                    "Expected overlap":         expected,
                    "p-value":                  p,
                })

            print(f"\n{comp_label}:")
            print(f"  Tissue-level DEPs                  n={s1:>4}  overlap={o1:>3}  p={p1:.3e}")
            print(f"  Tissue-level DEPs + connectors     n={s2:>4}  overlap={o2:>3}  p={p2:.3e}")
            print(f"  First-level DEPs only              n={s3:>4}  overlap={o3:>3}  p={p3:.3e}")

        # ── Figure ────────────────────────────────────────────────────────────

        x          = np.arange(len(labels))
        comparisons = list(PROTEIN_N.keys())
        n_comp     = len(comparisons)

        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.5, 4.5))
        fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.28, wspace=0.32)

        # ── Panel (a): protein set sizes ──────────────────────────────────────

        for ci, comp in enumerate(comparisons):
            offset    = (ci - (n_comp - 1) / 2) * (BAR_WIDTH + GROUP_GAP / 2)
            positions = x + offset
            ax_a.bar(positions, PROTEIN_N[comp], width=BAR_WIDTH,
                     color=COLORS[comp], alpha=0.90,
                     edgecolor="white", linewidth=0.6, zorder=2)
            for pos, n in zip(positions, PROTEIN_N[comp]):
                ax_a.text(pos, n + 3, f"n={n}", ha="center", va="bottom",
                          fontsize=7.5, color=COLORS[comp], fontweight="bold")

        ymax_a = max(v for vals in PROTEIN_N.values() for v in vals) * 1.28
        ax_a.set_ylabel("Number of proteins", fontsize=10, labelpad=6)
        ax_a.set_ylim(0, ymax_a)
        ax_a.set_xlim(-0.6, len(labels) - 0.4)
        ax_a.tick_params(axis="y", labelsize=9)
        ax_a.yaxis.grid(True, linestyle=":", linewidth=0.7, alpha=0.5, zorder=0)
        ax_a.set_axisbelow(True)
        ax_a.spines[["top", "right"]].set_visible(False)
        ax_a.set_xticks(x)
        ax_a.set_xticklabels(labels, fontsize=9)
        ax_a.text(-0.01, 1.02, "(a)", transform=ax_a.transAxes,
                  fontsize=11, fontweight="bold", va="bottom", ha="left")

        legend_handles = [
            mpatches.Patch(color=COLORS["Empty defect"], label="Empty defect"),
            mpatches.Patch(color=COLORS["PCL scaffold"],  label="PCL scaffold"),
        ]

        # ── Panel (b): overlap significance — lollipop chart ─────────────────

        all_neg_log = [-np.log10(p) for pvals in protein_pvals.values() for p in pvals]
        ymax_b      = max(all_neg_log) * 1.22

        for ci, comp in enumerate(comparisons):
            neg_log   = [-np.log10(p) for p in protein_pvals[comp]]
            offset    = (ci - (n_comp - 1) / 2) * (BAR_WIDTH + GROUP_GAP / 2)
            positions = x + offset
            for pos, val in zip(positions, neg_log):
                val_clipped = min(val, ymax_b)
                ax_b.scatter(pos, val_clipped, color=COLORS[comp],
                             s=100, zorder=3, edgecolors="white", linewidths=0.6)
                if val > ymax_b:
                    ax_b.annotate(
                        "", xy=(pos, ymax_b - 0.02), xytext=(pos, ymax_b - 0.4),
                        arrowprops=dict(arrowstyle="-|>", color=COLORS[comp], lw=2.0),
                        zorder=6,
                    )

        sig_line = -np.log10(0.05)
        ax_b.axhline(sig_line, color=LINE_COLOR, linewidth=1.8, linestyle="--", zorder=5)
        ax_b.text(-0.55, sig_line + ymax_b * 0.02, "p = 0.05",
                  ha="left", va="bottom", fontsize=8, color=LINE_COLOR)

        ax_b.set_ylabel("−log₁₀(p-value)\n(overlap with bone-healing reference proteins)", fontsize=10, labelpad=6)
        ax_b.set_ylim(0, ymax_b)
        ax_b.set_xlim(-0.6, len(labels) - 0.4)
        ax_b.tick_params(axis="y", labelsize=9)
        ax_b.yaxis.grid(True, linestyle=":", linewidth=0.7, alpha=0.5, zorder=0)
        ax_b.set_axisbelow(True)
        ax_b.spines[["top", "right"]].set_visible(False)
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(labels, fontsize=9)
        ax_b.text(-0.01, 1.02, "(b)", transform=ax_b.transAxes,
                  fontsize=11, fontweight="bold", va="bottom", ha="left")

        # ── Shared legend centered below both panels ──────────────────────────

        fig.legend(handles=legend_handles, fontsize=9, frameon=False,
                   loc="lower center", bbox_to_anchor=(0.53, 0.01), ncol=2,
                   handlelength=1.2, handletextpad=0.5, columnspacing=1.2,
                   title="Comparison",
                   title_fontproperties={"weight": "bold", "size": 10})

        # ── Save ──────────────────────────────────────────────────────────────

        out_path = os.path.join(THRESH_DIR, "validation_protein_overlap.png")
        plt.savefig(out_path, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"Saved {out_path}")

# ── Supplement export ─────────────────────────────────────────────────────────

PROTEIN_SET_ORDER = [
    "Tissue-level DEPs",
    "Tissue-level DEPs + connector proteins",
    "First-level DEPs",
]

results_df = pd.DataFrame(_results)

SUPP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "supplements")
os.makedirs(SUPP_DIR, exist_ok=True)
OUT_FILE = os.path.join(SUPP_DIR, "protein_bone_overlap_supplement.xlsx")
with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
    for comp_label in COMPARISONS.values():
        subset = results_df[results_df["Comparison"] == comp_label].drop(columns="Comparison")
        pivot  = subset.pivot(
            index=["FC threshold", "Stability threshold"],
            columns="Protein set",
            values=["Protein set size (N)", "Bone CAPs overlap (k)", "Expected overlap", "p-value"],
        )
        pivot = pivot.reindex(columns=pd.MultiIndex.from_product([
            ["Protein set size (N)", "Bone CAPs overlap (k)", "Expected overlap", "p-value"],
            PROTEIN_SET_ORDER,
        ]))
        pivot.columns = [f"{metric} | {pset}" for metric, pset in pivot.columns]
        pivot = pivot.reset_index()

        pivot.to_excel(writer, sheet_name=comp_label, index=False)

print(f"Saved → {OUT_FILE}")
print(f"  Sheets: {list(COMPARISONS.values())}")