"""
Exports the RobNorm normalized and RobNorm scaled-normalized abundance data and
metadata for AI and AS fractions into a single Excel file, for use as a supplement.

Output: normalized_data_supplement.xlsx
  - Sheet "AI normalized data"        : RobNorm-normalized log2 abundances for AI fraction
  - Sheet "AS normalized data"        : RobNorm-normalized log2 abundances for AS fraction
  - Sheet "AI scaled normalized data" : RobNorm + scaled log2 abundances for AI fraction
  - Sheet "AS scaled normalized data" : RobNorm + scaled log2 abundances for AS fraction
  - Sheet "AI metadata"               : sample metadata for AI fraction
  - Sheet "AS metadata"               : sample metadata for AS fraction
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

AI_DATA        = os.path.join(ROOT, "normalization_results_AI/RobNorm_normalized_data.csv")
AS_DATA        = os.path.join(ROOT, "normalization_results_AS/RobNorm_normalized_data.csv")
AI_DATA_SCALED = os.path.join(ROOT, "normalization_results_AI/scaled/RobNorm_scaled_normalized_data.csv")
AS_DATA_SCALED = os.path.join(ROOT, "normalization_results_AS/scaled/RobNorm_scaled_normalized_data.csv")
AI_META        = os.path.join(ROOT, "normalization_results_AI/meta_data.csv")
AS_META        = os.path.join(ROOT, "normalization_results_AS/meta_data.csv")

OUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "normalized_data_supplement.xlsx")

ai_data        = pd.read_csv(AI_DATA)
as_data        = pd.read_csv(AS_DATA)
ai_data_scaled = pd.read_csv(AI_DATA_SCALED)
as_data_scaled = pd.read_csv(AS_DATA_SCALED)
ai_meta        = pd.read_csv(AI_META)
as_meta        = pd.read_csv(AS_META)

README = pd.DataFrame([
    ("Sheet name",               "Content"),
    ("AI normalized data",       "RobNorm-normalized log2 abundances, AI fraction."),
    ("AS normalized data",       "RobNorm-normalized log2 abundances, AS fraction."),
    ("AI scaled normalized data","RobNorm-normalized and scaled log2 abundances, AI fraction."),
    ("AS scaled normalized data","RobNorm-normalized and scaled log2 abundances, AS fraction."),
    ("AI metadata",              "Sample metadata, AI fraction."),
    ("AS metadata",              "Sample metadata, AS fraction."),
    ("",                         ""),
    ("Column name",              "Description"),
    ("Protein.Group",            "Protein group identifier."),
    ("Protein.IDs",              "UniProt accessions."),
    ("Protein.Names",            "Protein names."),
    ("Genes",                    "Gene symbols."),
    ("First.Protein.Description","Leading protein description."),
    ("FilteredProteinIDs",       "Filtered protein IDs."),
    ("RemappedGeneNames",        "Remapped gene names."),
    ("Gene.Names",               "Gene symbols."),
    ("Orthologs",                "Human ortholog gene symbols."),
    ("<sample columns>",         "RobNorm-normalized log2 abundances (scaled in scaled sheets). Missing = not quantified."),
    ("",                         ""),
    ("Metadata column",          "Description"),
    ("sample_name",              "Sample identifier."),
    ("animal_no",                "Animal number."),
    ("fraction",                 "Proteomics fraction."),
    ("genotype",                 "Diabetic or non-diabetic."),
    ("scaffold",                 "Scaffold condition."),
    ("time",                     "Time point (days)."),
    ("intervention",             "Group label (genotype + scaffold + time)."),
], columns=["_", "__"])

with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
    README.to_excel(writer, sheet_name="README", index=False, header=False, startrow=2)
    writer.sheets["README"]["A1"] = "RobNorm-Normalized Protein Abundance Data"
    ai_data.to_excel(writer, sheet_name="AI normalized data",       index=False)
    as_data.to_excel(writer, sheet_name="AS normalized data",       index=False)
    ai_data_scaled.to_excel(writer, sheet_name="AI scaled normalized data", index=False)
    as_data_scaled.to_excel(writer, sheet_name="AS scaled normalized data", index=False)
    ai_meta.to_excel(writer, sheet_name="AI metadata",              index=False)
    as_meta.to_excel(writer, sheet_name="AS metadata",              index=False)

print(f"Saved → {OUT_FILE}")
print(f"  AI normalized data        : {len(ai_data)} proteins × {len(ai_data.columns)} columns")
print(f"  AS normalized data        : {len(as_data)} proteins × {len(as_data.columns)} columns")
print(f"  AI scaled normalized data : {len(ai_data_scaled)} proteins × {len(ai_data_scaled.columns)} columns")
print(f"  AS scaled normalized data : {len(as_data_scaled)} proteins × {len(as_data_scaled.columns)} columns")
print(f"  AI metadata               : {len(ai_meta)} samples")
print(f"  AS metadata               : {len(as_meta)} samples")