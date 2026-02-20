import pandas as pd
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
os.chdir('/home/ole/symbod_proteomics')
FC_THRESH = 5.6
STAB_THRESH = 5
BASE_DIR = 'valid_DEPs'

COMPARISONS = [
    'diabetic_empty_42-nondiabetic_empty_42',
    'diabetic_PCL_42-nondiabetic_PCL_42'
]

# ============================================================================
# ANALYZE TOP 10 DRUG OVERLAP ACROSS SETS
# ============================================================================

results_dir = os.path.join(BASE_DIR, f'FC{FC_THRESH}_Stab{STAB_THRESH}', 'drug_repurposing')

print(f"\n{'=' * 70}")
print("TOP 10 DRUG OVERLAP ANALYSIS")
print(f"{'=' * 70}")

# Store all results for table
all_comparison_results = []

for comparison in COMPARISONS:
    print(f"\n{comparison}")
    print("-" * 70)

    # Store data for each set
    set_data = {}

    for set_name in ['set1_first_level', 'set2_first_plus_network', 'set3_first_plus_second']:
        # Read the CSV file for this set
        csv_file = os.path.join(results_dir, f'{comparison}_{set_name}_drugs.csv')

        if not os.path.exists(csv_file):
            print(f"WARNING: File not found: {csv_file}")
            continue

        df = pd.read_csv(csv_file)
        set_data[set_name] = df

        # Get top 10 by score
        top10 = df.nlargest(10, 'Score')

        print(f"\n{set_name}: {len(top10)} drugs")
        for idx, row in top10.iterrows():
            print(f"  {row['Drug']:30s} (score={row['Score']:.4f}, targets={row['Num_Targets']})")

    # Get top 10 drugs for overlap analysis
    top10_by_set = {}
    for set_name, df in set_data.items():
        top10 = df.nlargest(10, 'Score')
        top10_by_set[set_name] = set(top10['Drug'].tolist())

    # Find intersection across all 3 sets
    if len(top10_by_set) == 3:
        set1 = top10_by_set['set1_first_level']
        set2 = top10_by_set['set2_first_plus_network']
        set3 = top10_by_set['set3_first_plus_second']

        # All 3 sets
        shared_all = set1 & set2 & set3

        print(f"\n{'=' * 70}")
        print("OVERLAP SUMMARY:")
        print(f"{'=' * 70}")
        print(f"Drugs in top 10 of ALL 3 sets: {len(shared_all)}")
        if shared_all:
            for drug in sorted(shared_all):
                print(f"  - {drug}")

        # Create results table for drugs in all 3 sets
        if shared_all:
            print(f"\n{'=' * 70}")
            print("DETAILED SCORES FOR SHARED DRUGS:")
            print(f"{'=' * 70}")

            for drug in sorted(shared_all):
                # Get scores from each set
                score_set1 = set_data['set1_first_level'][
                    set_data['set1_first_level']['Drug'] == drug]['Score'].values[0]
                score_set2 = set_data['set2_first_plus_network'][
                    set_data['set2_first_plus_network']['Drug'] == drug]['Score'].values[0]
                score_set3 = set_data['set3_first_plus_second'][
                    set_data['set3_first_plus_second']['Drug'] == drug]['Score'].values[0]

                avg_score = (score_set1 + score_set2 + score_set3) / 3

                # Get targets from set1 (should be same across sets)
                targets = set_data['set1_first_level'][
                    set_data['set1_first_level']['Drug'] == drug]['Targets'].values[0]
                num_targets = set_data['set1_first_level'][
                    set_data['set1_first_level']['Drug'] == drug]['Num_Targets'].values[0]
                status = set_data['set1_first_level'][
                    set_data['set1_first_level']['Drug'] == drug]['Status'].values[0]

                print(f"\n{drug}:")
                print(f"  Set1 score: {score_set1:.4f}")
                print(f"  Set2 score: {score_set2:.4f}")
                print(f"  Set3 score: {score_set3:.4f}")
                print(f"  Average:    {avg_score:.4f}")
                print(f"  Targets:    {num_targets} ({targets})")
                print(f"  Status:     {status}")

                # Store for table
                all_comparison_results.append({
                    'Comparison': comparison,
                    'Drug': drug,
                    'Set1_Score': score_set1,
                    'Set2_Score': score_set2,
                    'Set3_Score': score_set3,
                    'Average_Score': avg_score,
                    'Num_Targets': num_targets,
                    'Targets': targets,
                    'Status': status
                })

# Create final results table
if all_comparison_results:
    results_df = pd.DataFrame(all_comparison_results)
    results_df = results_df.sort_values(['Comparison', 'Average_Score'], ascending=[True, False])

    # Save to CSV
    output_file = os.path.join(results_dir, 'top10_shared_drugs_summary.csv')
    results_df.to_csv(output_file, index=False)

    print(f"\n{'=' * 70}")
    print("RESULTS TABLE SAVED")
    print(f"{'=' * 70}")
    print(f"File: {output_file}")
    print(f"\nDrugs appearing in top 10 of all 3 sets:")
    print(results_df[['Comparison', 'Drug', 'Average_Score', 'Num_Targets', 'Status']].to_string(index=False))
else:
    print(f"\n{'=' * 70}")
    print("No drugs found in top 10 of all 3 sets for any comparison")
    print(f"{'=' * 70}")

print(f"\n{'=' * 70}")
print("✓ Top 10 overlap analysis complete!")