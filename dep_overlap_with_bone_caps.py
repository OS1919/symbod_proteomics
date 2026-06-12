import pandas as pd
from scipy.stats import hypergeom
import os

# Configuration
base_path = 'valid_DEPs/FC5.59_Stab1.2'
comparison = "diabetic_PCL_42-nondiabetic_PCL_42"
#comparison = "diabetic_empty_42-nondiabetic_empty_42"


# Load both fraction files
ai_data = pd.read_csv('de_analysis_results_AI_RobNorm_scaled/de_results_raw.csv')
as_data = pd.read_csv('de_analysis_results_AS_RobNorm_scaled/de_results_raw.csv')

# Filter for specific comparison
ai_comp = ai_data[ai_data['Comparison'] == comparison]
as_comp = as_data[as_data['Comparison'] == comparison]

# Get unique proteins from both fractions
ai_proteins = set(ai_comp['Protein.IDs'].dropna())
as_proteins = set(as_comp['Protein.IDs'].dropna())

# Combine to get all tested proteins
all_tested_proteins = ai_proteins.union(as_proteins)

print(f"AI fraction: {len(ai_proteins)} proteins")
print(f"AS fraction: {len(as_proteins)} proteins")
print(f"Total unique tested: {len(all_tested_proteins)} proteins")

# This is your actual proteome size
total_proteome_size = len(all_tested_proteins)

# Load data files
summary_stats = pd.read_csv(os.path.join(base_path, 'summary_statistics.csv'))
exception = pd.read_csv('network_enrichment_results/' + comparison + '/exception_proteins.csv')
bone_caps = pd.read_csv('input/bone_caps_meta_analysis.csv')

# Prepare exception proteins set
exception_prots = set(exception['Gene'].dropna().str.upper())
bone_caps_prots = set(bone_caps['Gene'].dropna().str.upper())

overlap_with_exceptions = len(exception_prots.intersection(bone_caps_prots))
print('Overlap with exception proteins:', overlap_with_exceptions)

print(f"Total proteome size: {total_proteome_size}")
print(f"Known bone CAPs: {len(bone_caps_prots)}")
print(f"Network exception proteins: {len(exception_prots)}\n")

# Filter for specific comparison
row = summary_stats[summary_stats['Comparison'] == comparison]

if row.empty:
    print(f"ERROR: Comparison '{comparison}' not found!")
else:
    row = row.iloc[0]

    bone_caps_tested = int(row['Bone_Caps_Tested'])
    print('Bone CAPs tested', bone_caps_tested)

    # SET 1: First-level DEPs only
    set1_count = int(row['Validated_Both_Fractions'] + row['Validated_Exclusive'])
    print(f"1st-level valid DEPs: {set1_count}")
    set1_overlap = int(row['Overlap_First_Level'])
    print(f"1st-level valid DEPs bone CAP overlap: {set1_overlap}")
    # hypergeom.sf(k, M, n, N) with n the number of successes in M and N the number of draws
    set1_pval = hypergeom.sf(set1_overlap - 1, total_proteome_size, bone_caps_tested, set1_count)

    print(f"Set 1 (First-level only):")
    print(f"  Total: {set1_count}")
    print(f"  Overlap: {set1_overlap} ({set1_overlap/set1_count*100:.2f}%)")
    print(f"  P-value: {set1_pval:.2e}")

    # SET 2: First-level + Network exception proteins
    set2_count = set1_count + len(exception_prots)
    set2_overlap = set1_overlap  + overlap_with_exceptions
    set2_pval = hypergeom.sf(set2_overlap - 1, total_proteome_size, bone_caps_tested, set2_count)

    print(f"\nSet 2 (First-level + Network):")
    print(f"  Total: {set2_count} ({len(exception_prots)} from network)")
    print(f"  Overlap: {set2_overlap} ({set2_overlap / set2_count * 100:.2f}%)")
    print(f"  P-value: {set2_pval:.2e}")

    # SET 3: First-level + Second-level DEPs
    set3_count = set1_count + int(row['Validated_Dominant'])
    set3_overlap = set1_overlap + int(row['Overlap_Second_Level'])
    set3_pval = hypergeom.sf(set3_overlap - 1, total_proteome_size, bone_caps_tested, set3_count)

    print(f"\nSet 3 (First-level + Second-level):")
    print(f"  Total: {set3_count} (+{int(row['Validated_Dominant'])} from second-level)")
    print(f"  Overlap: {set3_overlap} ({set3_overlap / set3_count * 100:.2f}%)")
    print(f"  P-value: {set3_pval:.2e}")