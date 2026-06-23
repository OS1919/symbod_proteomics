"""
Exports the validated DEPs for every threshold combination and comparison into
a single Excel file, one sheet per (threshold combination × comparison).

Output: validated_deps_supplement.xlsx  — 18 sheets (9 thresholds × 2 comparisons)
Sheet naming: "FC{fc} Stab{stab} {comparison label}"
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]

COMPARISONS = {
    "diabetic_empty_42-nondiabetic_empty_42": "Empty defect",
    "diabetic_PCL_42-nondiabetic_PCL_42":     "PCL scaffold",
}

OUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tissue_level_deps_supplement.xlsx")

README = pd.DataFrame([
    ("Sheet name",           "Content"),
    ("FC{fc} Stab{stab} {comparison}", "e.g. FC5.59 Stab1.2 Empty defect. One sheet per threshold × comparison (18 total)."),
    ("",                     ""),
    ("Threshold",            "Description"),
    ("{fc}",                 "FC threshold: minimum comparison-level abundance ratio between fractions required to qualify a protein as a second-level DEP."),
    ("{stab}",               "Stability threshold: minimum ratio of group-level geometric mean FCs required for a second-level DEP to be considered stable across groups."),
    ("",                     ""),
    ("Column name",          "Description"),
    ("Protein.IDs",          "UniProt accessions."),
    ("Gene.Names",           "Gene symbols."),
    ("Orthologs",            "Human ortholog gene symbols."),
    ("Direction",            "Up or Down."),
    ("DEP_in_AI",            "Boolean whether protein is DEP in AI."),
    ("DEP_in_AS",            "Boolean whether protein is DEP in AS."),
    ("FC_between_fractions", "Fold change between fractions."),
    ("Dominant_Fraction",    "Dominant fraction (AI or AS) according to comparison-level abundance ratio."),
    ("Validation",           "Validation reason for the protein being a tissue-level DEP."),
    ("GeoMeanFC_Group1",     "Geometric mean FC across diabetic samples of the comparison."),
    ("GeoMeanFC_Group2",     "Geometric mean FC across non-diabetic samples of the comparison."),
    ("GeoMeanFC_Ratio",      "max(GeoMeanFC_Group1, GeoMeanFC_Group2) / min(GeoMeanFC_Group1, GeoMeanFC_Group2)"),
], columns=["_", "__"])

with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
    README.to_excel(writer, sheet_name="README", index=False, header=False, startrow=2)
    writer.sheets["README"]["A1"] = "Validated Differentially Expressed Proteins (DEPs)"
    for fc in FC_THRESHOLDS:
        for stab in STAB_THRESHOLDS:
            path = os.path.join(ROOT, "valid_DEPs", f"FC{fc}_Stab{stab}", "validated_DEPs.csv")
            df   = pd.read_csv(path)

            for comp_key, comp_label in COMPARISONS.items():
                subset     = df[df["Comparison"] == comp_key].drop(columns=["Comparison", "Comparison_Label"])
                sheet_name = f"FC{fc} Stab{stab} {comp_label}"
                subset.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  {sheet_name}: {len(subset)} validated DEPs")

print(f"\nSaved → {OUT_FILE}")