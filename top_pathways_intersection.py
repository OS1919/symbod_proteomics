import pandas as pd
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
FC_THRESH = 5.6
STAB_THRESH = 1.2
BASE_DIR = "valid_DEPs"
BONE_FILE = "input/bone_enrichments_meta_analysis.csv"

# Comparisons to analyze
COMPARISONS = [
    "diabetic_empty_42-nondiabetic_empty_42",
    "diabetic_PCL_42-nondiabetic_PCL_42"
]
# ============================================================================

thresh_dir = os.path.join(BASE_DIR, f"FC{FC_THRESH}_Stab{STAB_THRESH}", "enrichment")

# Load bone pathways
bone_df = pd.read_csv(BONE_FILE)
bone_pathway_ids = set(bone_df['term_id'].dropna())

sets = ['set1_first_level', 'set2_first_plus_network', 'set3_first_plus_second']

# Collect all results
all_results = []

for comparison in COMPARISONS:
    print(f"\n{'=' * 70}")
    print(f"COMPARISON: {comparison}")
    print(f"{'=' * 70}")

    comp_top10 = {}

    for set_name in sets:
        # Load enrichment results
        filepath = os.path.join(thresh_dir, set_name,
                                f"{comparison}_{set_name}_enrichment_complete.csv")

        if not os.path.exists(filepath):
            print(f"WARNING: File not found: {filepath}")
            continue

        df = pd.read_csv(filepath)

        # Get top 10 significant pathways
        sig_df = df[df['p.adjust'] < 0.05].copy()
        sig_df = sig_df.sort_values('p.adjust').reset_index(drop=True)
        top10 = sig_df.head(10)

        comp_top10[set_name] = top10

        print(f"\n{set_name}: {len(sig_df)} significant, Top 10:")
        for i, row in top10.iterrows():
            bone_marker = "★" if row['ID'] in bone_pathway_ids else " "
            print(f"  {i + 1:2d}. {bone_marker} {row['Description'][:50]:50s} p={row['p.adjust']:.2e}")

    # Find pathways in all 3 sets (robust)
    if len(comp_top10) == 3:
        ids_set1 = set(comp_top10['set1_first_level']['ID'])
        ids_set2 = set(comp_top10['set2_first_plus_network']['ID'])
        ids_set3 = set(comp_top10['set3_first_plus_second']['ID'])

        shared_all = ids_set1 & ids_set2 & ids_set3

        print(f"\n{'=' * 70}")
        print(f"PATHWAYS IN TOP 10 OF ALL 3 SETS: {len(shared_all)}")
        print(f"{'=' * 70}")

        if shared_all:
            for pathway_id in shared_all:
                # Get description from set1
                row = comp_top10['set1_first_level'][
                    comp_top10['set1_first_level']['ID'] == pathway_id].iloc[0]
                bone_marker = "★ BONE" if pathway_id in bone_pathway_ids else ""

                print(f"\n{row['Description']} {bone_marker}")
                print(f"  ID: {pathway_id}, Category: {row['Category']}")

                # Show p-values across sets
                for set_name in sets:
                    p_val = comp_top10[set_name][
                        comp_top10[set_name]['ID'] == pathway_id]['p.adjust'].values[0]
                    print(f"    {set_name:30s}: p={p_val:.2e}")

                # Store for output
                all_results.append({
                    'Comparison': comparison,
                    'Pathway_ID': pathway_id,
                    'Description': row['Description'],
                    'Category': row['Category'],
                    'In_Bone_Meta': pathway_id in bone_pathway_ids,
                    'Set1_pvalue': comp_top10['set1_first_level'][
                        comp_top10['set1_first_level']['ID'] == pathway_id]['p.adjust'].values[0],
                    'Set2_pvalue': comp_top10['set2_first_plus_network'][
                        comp_top10['set2_first_plus_network']['ID'] == pathway_id]['p.adjust'].values[0],
                    'Set3_pvalue': comp_top10['set3_first_plus_second'][
                        comp_top10['set3_first_plus_second']['ID'] == pathway_id]['p.adjust'].values[0]
                })

# Save results
output_df = pd.DataFrame(all_results)
output_file = os.path.join(thresh_dir, "top10_shared_pathways_summary.csv")
output_df.to_csv(output_file, index=False)

print(f"\n\n{'=' * 70}")
print(f"SUMMARY SAVED TO: {output_file}")
print(f"{'=' * 70}")
print(f"\nTotal robust pathways (in all 3 sets): {len(output_df)}")
print(f"Validated by bone meta-analysis: {output_df['In_Bone_Meta'].sum()}")

# ============================================================================
# COMBINED HEATMAP FIGURE
# ============================================================================
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import numpy as np

output_df = pd.DataFrame(all_results)

col_labels = ["First-level DEPs", "First-level DEPs +\nconnector\nproteins", "First-level DEPs +\nsecond-level DEPs + \n(AR ≥ 5.6 & ΔAR ≤ 1.2)"]
comp_labels = {"diabetic_empty_42-nondiabetic_empty_42": "Empty defect",
               "diabetic_PCL_42-nondiabetic_PCL_42": "PCL scaffold"}

# Get all unique terms across both comparisons, preserving order by Set1 p-value
all_terms = []
for comparison in COMPARISONS:
    comp_df = output_df[output_df['Comparison'] == comparison].sort_values('Set1_pvalue')
    for _, row in comp_df.iterrows():
        if row['Pathway_ID'] not in [t['id'] for t in all_terms]:
            all_terms.append({
                'id': row['Pathway_ID'],
                'description': row['Description'],
                'bone': row['In_Bone_Meta'],
            })

# Build one matrix per comparison (rows = all_terms, cols = 3 sets)
matrices = {}
for comparison in COMPARISONS:
    comp_df = output_df[output_df['Comparison'] == comparison]
    mat = np.full((len(all_terms), 3), np.nan)
    for i, term in enumerate(all_terms):
        row = comp_df[comp_df['Pathway_ID'] == term['id']]
        if not row.empty:
            mat[i, 0] = -np.log10(row['Set1_pvalue'].values[0])
            mat[i, 1] = -np.log10(row['Set2_pvalue'].values[0])
            mat[i, 2] = -np.log10(row['Set3_pvalue'].values[0])
    matrices[comparison] = mat

# Shared color scale across both panels
vmax = max(np.nanmax(m) for m in matrices.values()) * 1.05

# Y-axis labels with bone marker
y_labels = [
    ("★ " if t['bone'] else "   ") + t['description']
    for t in all_terms
]

n_terms = len(all_terms)

row_labels = [
    "First-level DEPs",
    "First-level DEPs +\nconnector proteins",
    "First-level DEPs +\nsecond-level DEPs\n(AR ≥ 5.6 & ΔAR ≤ 1.2)",
]
x_labels = [("★ " if t['bone'] else "") + t['description'] for t in all_terms]

fig, axes = plt.subplots(2, 1, figsize=(max(12, n_terms * 2.0), 14),
                         sharex=True)
fig.subplots_adjust(hspace=0.35, bottom=0.30, left=0.22, right=0.93, top=0.95)

for ax, comparison in zip(axes, COMPARISONS):
    mat = matrices[comparison].T  # (3 sets, n_terms pathways)

    bg = np.zeros_like(mat)
    ax.imshow(bg, aspect='auto', cmap='Greys', vmin=0, vmax=1, alpha=0.15)

    masked = np.ma.masked_invalid(mat)
    im = ax.imshow(masked, aspect='auto', cmap='Blues', vmin=0, vmax=vmax)

    ax.set_yticks(range(3))
    ax.set_yticklabels(row_labels, fontsize=12)
    ax.set_title(comp_labels[comparison], fontsize=14, fontweight='bold', pad=8)

    # Cell annotations
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            if np.isnan(val):
                ax.text(j, i, "x", ha='center', va='center',
                        fontsize=11, color='#aaaaaa')
            else:
                ax.text(j, i, f"{val:.1f}", ha='center', va='center',
                        fontsize=12,
                        color='white' if val > vmax * 0.6 else '#333333')

# X-labels only on bottom panel
axes[1].set_xticks(range(n_terms))
axes[1].set_xticklabels(x_labels, fontsize=12, rotation=45, ha='right')

# Shared colorbar
sm = mcm.ScalarMappable(cmap='Blues', norm=mcolors.Normalize(vmin=0, vmax=vmax))
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02)
cbar.set_label("−log₁₀(adjusted p-value)", fontsize=12)
cbar.ax.tick_params(labelsize=11)

axes[1].annotate(
    "★ Validated by bone healing meta-analysis (Schmidt et al.)  |  x = term not in top 10 for this comparison",
    xy=(0, 0), xycoords='axes fraction',
    xytext=(0, -0.55), textcoords='axes fraction',
    fontsize=12, style='italic', annotation_clip=False, va='top',
)

out_path = os.path.join(thresh_dir, "heatmap_combined.png")
plt.savefig(out_path, dpi=180, bbox_inches='tight')
plt.close()
print(f"Saved: {out_path}")