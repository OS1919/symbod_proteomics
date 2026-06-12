import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from sklearn.metrics import auc

# ============================================================================
# CONFIGURATION - CHANGE THESE
# ============================================================================
FC_THRESH = 5.59
STAB_THRESH = 1.2
BASE_DIR = "valid_DEPs"
BONE_FILE = "input/bone_enrichments_meta_analysis.csv"
# ============================================================================

thresh_dir = os.path.join(BASE_DIR, f"FC{FC_THRESH}_Stab{STAB_THRESH}", "enrichment")

# Read bone meta-analysis pathways
bone_df = pd.read_csv(BONE_FILE)
bone_pathways = set(bone_df['term_id'].dropna())

print(f"\n{'=' * 60}")
print(f"Bone meta-analysis: {len(bone_pathways)} pathways")
print(f"{'=' * 60}\n")

# Loop through set directories (set1, set2, set3)
for set_dir in ['set1_first_level', 'set2_first_plus_network', 'set3_first_plus_second']:
    set_path = os.path.join(thresh_dir, set_dir)

    if not os.path.isdir(set_path):
        continue

    print(f"\n{'=' * 60}")
    print(f"{set_dir.upper()}")
    print(f"{'=' * 60}")

    # Find all enrichment files in this set directory
    for filename in os.listdir(set_path):
        if not filename.endswith('_enrichment_complete.csv'):
            continue

        # Extract comparison name from filename
        # Format: comparison_setname_enrichment_complete.csv
        comparison = filename.replace(f'_{set_dir}_enrichment_complete.csv', '')

        print(f"\n  Comparison: {comparison}")
        print("-" * 40)

        complete_file = os.path.join(set_path, filename)
        complete_df = pd.read_csv(complete_file)

        # Get significant pathways (p.adjust < 0.05)
        sig_df = complete_df[complete_df['p.adjust'] < 0.05].copy()

        sig_pathways = set(sig_df['ID'].dropna())
        overlap_sig = bone_pathways & sig_pathways

        # All tested pathways
        complete_pathways = set(complete_df['ID'].dropna())
        overlap_complete = bone_pathways & complete_pathways

        # ================================================================
        # HYPERGEOMETRIC TEST PARAMETERS
        # ================================================================
        M = len(complete_pathways)
        n = len(overlap_complete)
        N = len(sig_pathways)
        k = len(overlap_sig)

        expected_overlap = (N * n) / M if M > 0 else 0
        # Question: "What's the probability of getting k or MORE bone pathways by chance?"
        # This is P(X ≥ k) where X ~ Hypergeometric(M, n, N)
        #
        # sf(k-1) gives P(X > k-1) = P(X ≥ k)
        # This is a ONE-TAILED test for ENRICHMENT (over-representation)
        # Under the null hypothesis (no enrichment), we'd expect:
        # E[X] = N × (n/M) red marbles in our draw
        p_value = hypergeom.sf(k - 1, M, n, N)
        fold_enrichment = (k / expected_overlap) if expected_overlap > 0 else float('inf')

        # ================================================================
        # PRINT RESULTS
        # ================================================================
        print(f"  ────────────────────────────────────────")
        print(f"  HYPERGEOMETRIC MODEL:")
        print(f"Terms tested:        {M}")
        print(f"Bone-healing terms in tested terms:      {n}")
        print(f"Significant terms:   {N}")
        print(f"Bone-healing terms in significant terms: {k}")
        print(f"p-value:          {p_value:.2e}")

print(f"\n{'=' * 60}")
print("DONE")
print(f"{'=' * 60}\n")