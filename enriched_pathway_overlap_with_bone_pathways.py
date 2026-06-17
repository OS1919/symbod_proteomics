import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import hypergeom

# ── Config ────────────────────────────────────────────────────────────────────

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]
BASE_DIR  = "valid_DEPs"
BONE_FILE = "input/bone_enrichments_meta_analysis.csv"

COMPARISONS = {
    "diabetic_empty_42-nondiabetic_empty_42": "Empty defect",
    "diabetic_PCL_42-nondiabetic_PCL_42":     "PCL scaffold",
}

SET_DIRS = ["set1_tissue_level", "set2_tissue_plus_network", "set3_first_level"]

labels = [
    "Pathways enriched\nfrom\ntissue-level DEPs",
    "Pathways enriched\nfrom\ntissue-level DEPs\n+ connector proteins",
    "Pathways enriched\nfrom\nfirst-level DEPs",
]

COLORS = {
    "Empty defect": "#4A90C4",
    "PCL scaffold":  "#E07B54",
}
LINE_COLOR = "#1A1A1A"
BAR_WIDTH  = 0.32
GROUP_GAP  = 0.08

# Load bone reference pathways once (shared across all threshold combinations)
bone_pathways = set(pd.read_csv(BONE_FILE)["term_id"].dropna())

print(f"Bone meta-analysis: {len(bone_pathways)} pathways\n")

# ── Loop over all threshold combinations ──────────────────────────────────────

for fc in FC_THRESHOLDS:
    for stab in STAB_THRESHOLDS:
        thresh_dir = os.path.join(BASE_DIR, f"FC{fc}_Stab{stab}", "enrichment")

        if not os.path.isdir(thresh_dir):
            print(f"Skipping FC{fc}_Stab{stab} — enrichment directory not found")
            continue

        print(f"\n{'='*60}")
        print(f"FC={fc}  Stab={stab}  ({thresh_dir})")
        print(f"{'='*60}")

        # ── Collect counts and p-values per comparison × set ─────────────────

        PATHWAY_N     = {}
        pathway_pvals = {}

        for comp_key, comp_label in COMPARISONS.items():
            n_sig_list = []
            p_list     = []

            for set_dir in SET_DIRS:
                fname = f"{comp_key}_{set_dir}_enrichment_complete.csv"
                fpath = os.path.join(thresh_dir, set_dir, fname)

                if not os.path.exists(fpath):
                    print(f"  Missing: {fpath}")
                    n_sig_list.append(0)
                    p_list.append(1.0)
                    continue

                df     = pd.read_csv(fpath)
                sig_df = df[df["p.adjust"] < 0.05]

                complete_pathways = set(df["ID"].dropna())
                sig_pathways      = set(sig_df["ID"].dropna())
                overlap_complete  = bone_pathways & complete_pathways
                overlap_sig       = bone_pathways & sig_pathways

                M = len(complete_pathways)
                n = len(overlap_complete)
                N = len(sig_pathways)
                k = len(overlap_sig)

                expected = (N * n) / M if M > 0 else 0
                p        = hypergeom.sf(k - 1, M, n, N)

                n_sig_list.append(N)
                p_list.append(p)

                print(f"\n  {comp_label} / {set_dir}:")
                print(f"    Terms tested={M}  bone in tested={n}  significant={N}  bone in significant={k}  p={p:.2e}")

            PATHWAY_N[comp_label]     = n_sig_list
            pathway_pvals[comp_label] = p_list

        # ── Figure ────────────────────────────────────────────────────────────

        x          = np.arange(len(labels))
        comparisons = list(PATHWAY_N.keys())
        n_comp     = len(comparisons)

        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.5, 4.5))
        fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.28, wspace=0.32)

        # ── Panel (a): significant pathway counts ─────────────────────────────

        for ci, comp in enumerate(comparisons):
            offset    = (ci - (n_comp - 1) / 2) * (BAR_WIDTH + GROUP_GAP / 2)
            positions = x + offset
            ax_a.bar(positions, PATHWAY_N[comp], width=BAR_WIDTH,
                     color=COLORS[comp], alpha=0.90,
                     edgecolor="white", linewidth=0.6, zorder=2)
            for pos, n in zip(positions, PATHWAY_N[comp]):
                ax_a.text(pos, n + 3, f"n={n}", ha="center", va="bottom",
                          fontsize=7.5, color=COLORS[comp], fontweight="bold")

        ymax_a = max(v for vals in PATHWAY_N.values() for v in vals) * 1.28
        ymax_a = max(ymax_a, 10)  # minimum sensible axis height
        ax_a.set_ylabel("Number of significant pathways", fontsize=10, labelpad=6)
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

        all_neg_log = [-np.log10(p) for pvals in pathway_pvals.values() for p in pvals]
        ymax_b      = max(max(all_neg_log) * 1.22, -np.log10(0.05) * 1.5)

        for ci, comp in enumerate(comparisons):
            neg_log   = [-np.log10(p) for p in pathway_pvals[comp]]
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

        ax_b.set_ylabel("−log₁₀(p-value)\n(overlap with bone-healing reference pathways)", fontsize=10, labelpad=6)
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

        out_path = os.path.join(BASE_DIR, f"FC{fc}_Stab{stab}", "enriched_pathway_overlap.png")
        plt.savefig(out_path, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"\nSaved {out_path}")

print(f"\n{'='*60}\nDONE\n{'='*60}")