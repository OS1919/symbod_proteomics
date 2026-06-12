import drugstone
import pandas as pd
import json
import time
import os

drugstone.set_api('https://api.stable.drugst.one/')
drugstone.print_license()
drugstone.accept_license()

# Load the data
df = pd.read_csv('valid_DEPs/FC5.59_Stab1.2/validated_DEPs.csv')

print(f"Total proteins in file: {len(df)}")
print(f"\nValidation categories:")
print(df['Validation'].value_counts())

# Filter by multiple validation types
validation_filters = [
     "Significant in both fractions (same direction)",
     "Exclusive to AI fraction",
     "Exclusive to AS fraction"
]
comparison_filter = "diabetic_PCL_42-nondiabetic_PCL_42"
#comparison_filter = "diabetic_empty_42-nondiabetic_empty_42"
filtered_df = df[
    (df['Validation'].isin(validation_filters)) &
    (df['Comparison'] == comparison_filter)
]

# Extract gene names from filtered dataframe
my_genes = [gene.split(';')[0] for gene in filtered_df['Orthologs'].dropna().tolist()]
my_genes = list(set(my_genes))  # Remove duplicates

print(f"Number of filtered proteins: {len(filtered_df)}")
print(f"Number of seed proteins: {len(my_genes)}")

output_dir = os.path.join('network_enrichment_results', comparison_filter)
os.makedirs(output_dir, exist_ok=True)

enrichment_parameters = {
    "target": "drug-target", # specifies which type of nodes (proteins, drugs, or both) the algorithm should include in its output or even input
    "algorithm": "multisteiner",  # Multi-Steiner Tree algorithm (as in meta-analysis paper)
    "ppiDataset": "NeDRex",
    "num_trees": 5,  # Number of Steiner trees (meta-analysis paper uses 5)
    "tolerance": 5,  # Tolerance parameter (meta-analysis paper uses 5)
    "hubPenalty": 0.5,  # Hub penalty (meta-analysis paper uses 0.5)
    "resultSize": 500  # Allow more proteins to be added
}

print("\n=== Performing Network Enrichment ===")

# Create task for network enrichment
enrichment_task = drugstone.new_task(my_genes, enrichment_parameters)

# Get enriched network result
enrichment_result = enrichment_task.get_result()

# Save enriched network
os.chdir(output_dir)
enrichment_result.download_json()
enrichment_result.download_graph()
os.chdir('../..')  # Changed from '..' to '../..'

print("Waiting for enrichment files to be written...")
time.sleep(3)

# Read the enriched network
json_path = os.path.join(output_dir, 'result.json')
with open(json_path, 'r') as f:
    enriched_data = json.load(f)

# Analyze the enriched network
seed_proteins_lower = set(gene.lower() for gene in my_genes)
exception_proteins = []

print("\n=== Network Enrichment Results ===")

# Extract all proteins (nodes) from the network
if 'genes' in enriched_data:
    all_nodes = enriched_data['genes']
    print(f"Total proteins in enriched network: {len(all_nodes)}")

    # Identify exception proteins (added by algorithm)
    for node_id, node_info in all_nodes.items():
        protein_name = node_info.get('label', node_id)
        if protein_name.lower() not in seed_proteins_lower:
            exception_proteins.append({
                'Gene': protein_name,
                'Node_ID': node_id,
                'Type': 'Exception'
            })

    print(f"Seed proteins (input): {len(my_genes)}")
    print(f"Exception proteins (added): {len(exception_proteins)}")

    # Save exception proteins
    exception_df = pd.DataFrame(exception_proteins)
    exception_csv_path = os.path.join(output_dir, 'exception_proteins.csv')
    exception_df.to_csv(exception_csv_path, index=False)

    # Create a comprehensive network summary
    network_summary = {
        'Total_Nodes': len(all_nodes),
        'Seed_Proteins': len(my_genes),
        'Exception_Proteins': len(exception_proteins),
        'Total_Edges': len(enriched_data.get('edges', {}))
    }

    # Save network summary to file
    summary_path = os.path.join(output_dir, 'network_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("=== Network Summary ===\n")
        for key, value in network_summary.items():
            f.write(f"{key.replace('_', ' ')}: {value}\n")
    print("✓ Network summary saved to 'network_summary.txt'")
