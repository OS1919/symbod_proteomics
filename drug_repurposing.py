import drugstone
import pandas as pd
import json
import time
import os

drugstone.set_api('https://api.stable.drugst.one/')
drugstone.print_license()
drugstone.accept_license()

# ============================================================================
# CONFIGURATION
# ============================================================================
FC_THRESH = 5.6
STAB_THRESH = 1.2
BASE_DIR = 'valid_DEPs'

# Comparisons to analyze
COMPARISONS = [
    'diabetic_empty_42-nondiabetic_empty_42',
    'diabetic_PCL_42-nondiabetic_PCL_42'
]

# Drug search parameters
PARAMETERS = {
    "target": "drug",
    "algorithm": "trustrank",
    "ppi_dataset": "nedrex",
    "damping_factor": 0.85,
    "includeIndirectDrugs": False,
    "includeNonApprovedDrugs": True,
    "filterPaths": True, # include only shortest connections in the result
    "resultSize": 100
}


# ============================================================================

def get_genes(comparison, protein_set, base_dir, fc_thresh, stab_thresh):
    """Get genes for a specific protein set"""

    # Load validated DEPs
    validated_file = os.path.join(base_dir, f'FC{fc_thresh}_Stab{stab_thresh}', 'validated_DEPs.csv')
    df = pd.read_csv(validated_file)
    df = df[df['Comparison'] == comparison]

    if protein_set == 'set1_first_level':
        # First-level only
        df = df[df['Validation'].str.contains('Exclusive to |(same direction)', regex=True)]
        genes = [gene.split(';')[0] for gene in df['Orthologs'].dropna().tolist()]

    elif protein_set == 'set2_first_plus_network':
        # First-level + exception proteins
        df = df[df['Validation'].str.contains('Exclusive to |(same direction)', regex=True)]
        genes = [gene.split(';')[0] for gene in df['Orthologs'].dropna().tolist()]

        # Add exception proteins
        exception_file = f'network_enrichment_results/{comparison}/exception_proteins.csv'
        if os.path.exists(exception_file):
            exception_df = pd.read_csv(exception_file)
            exception_genes = exception_df['Gene'].dropna().tolist()
            genes.extend(exception_genes)
        else:
            print(f"WARNING: Exception file not found: {exception_file}")

    elif protein_set == 'set3_first_plus_second':
        # First-level + second-level (no filtering)
        genes = [gene.split(';')[0] for gene in df['Orthologs'].dropna().tolist()]

    else:
        raise ValueError(f"Unknown protein set: {protein_set}")

    genes = list(set(genes))  # Remove duplicates
    return genes


def run_drug_search(genes, comparison, set_name, output_dir):
    """Run drug repurposing search"""

    print(f"\n{'=' * 70}")
    print(f"Running drug search: {comparison} - {set_name}")
    print(f"{'=' * 70}")
    print(f"Number of seed genes: {len(genes)}")

    # Run DrugStone
    task = drugstone.new_task(genes, PARAMETERS)
    result = task.get_result()

    # Save files with unique names
    json_file = os.path.join(output_dir, f'{comparison}_{set_name}_drugs.json')

    # Save directly
    result.download_json()

    # Wait and rename
    time.sleep(2)
    os.rename('result.json', json_file)

    print(f"Saved: {json_file}")

    # Parse results
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Extract drug information
    drug_data = []
    for drug_name, drug_info in data['drugs'].items():
        drug_record = {
            'Comparison': comparison,
            'Set': set_name,
            'Drug': drug_info['label'],
            'Score': drug_info['score'],
            'DrugBank_ID': drug_info['drugId'],
            'Status': drug_info['status'],
            'Num_Targets': len(drug_info['hasEdgesTo']),
            'Targets': ', '.join(drug_info['hasEdgesTo']),
            'Score_Raw': drug_info['score_raw']
        }
        drug_data.append(drug_record)

    drug_df = pd.DataFrame(drug_data)
    drug_df = drug_df.sort_values(['Score', 'Num_Targets'], ascending=[False, False])

    # Save CSV
    csv_file = os.path.join(output_dir, f'{comparison}_{set_name}_drugs.csv')
    drug_df.to_csv(csv_file, index=False)
    print(f"Saved: {csv_file}")

    # Print summary
    print(f"\nSummary:")
    print(f"  Total drugs: {len(drug_df)}")
    print(f"  Multi-target drugs (3+): {len(drug_df[drug_df['Num_Targets'] >= 3])}")
    print(
        f"  Top drug: {drug_df.iloc[0]['Drug']} (score={drug_df.iloc[0]['Score']:.4f}, targets={drug_df.iloc[0]['Num_Targets']})")

    return drug_df


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

output_dir = os.path.join(BASE_DIR, f'FC{FC_THRESH}_Stab{STAB_THRESH}', 'drug_repurposing')
os.makedirs(output_dir, exist_ok=True)

all_results = []

for comparison in COMPARISONS:
    for set_name in ['set1_first_level', 'set2_first_plus_network', 'set3_first_plus_second']:
        try:
            # Get genes
            genes = get_genes(comparison, set_name, BASE_DIR, FC_THRESH, STAB_THRESH)

            if len(genes) == 0:
                print(f"WARNING: No genes found for {comparison} - {set_name}")
                continue

            # Run search
            drug_df = run_drug_search(genes, comparison, set_name, output_dir)
            all_results.append(drug_df)

        except Exception as e:
            print(f"ERROR in {comparison} - {set_name}: {e}")
            continue

# Combine all results
if all_results:
    combined_df = pd.concat(all_results, ignore_index=True)
    combined_file = os.path.join(output_dir, 'all_drugs_combined.csv')
    combined_df.to_csv(combined_file, index=False)
    print(f"\n{'=' * 70}")
    print(f"All results saved to: {combined_file}")
    print(f"{'=' * 70}")

    # Summary across all analyses
    print(f"\n=== OVERALL SUMMARY ===")
    print(f"Total analyses: {len(all_results)}")
    print(f"Total unique drugs: {combined_df['Drug'].nunique()}")
    print(f"\nDrugs by comparison:")
    print(combined_df.groupby('Comparison')['Drug'].nunique())
    print(f"\nDrugs by set:")
    print(combined_df.groupby('Set')['Drug'].nunique())

print("\n✓ Drug repurposing analysis complete!")