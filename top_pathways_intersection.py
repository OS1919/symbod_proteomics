import pandas as pd
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
FC_THRESH = 5.6
STAB_THRESH = 5
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