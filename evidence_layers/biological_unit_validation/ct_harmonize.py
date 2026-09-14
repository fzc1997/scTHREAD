#!/usr/bin/env python3
"""
Cross-study cell-type label harmonization for the two marrow studies (GSE307660, GSE276974).

The two studies annotate the same lineages with different vocabularies. For CROSS-STUDY effect
sizes and shared-cell-type claims we coarsen both to the common vocabulary (7 shared lineages).
The within-donor permutation TEST is robust to vocab either way (each donor is one study, internally
consistent), but harmonized labels make the effect size and the "shared across studies" narrative correct.

Shared coarse vocabulary: B cell, T cell, NK, Monocyte, Dendritic cell, Erythroid, Plasma cell.
Study-private progenitors (HSPC, Granulo/Myeloid progenitor) map to 'Progenitor' (present only in
GSE276974; kept for within-study tests, dropped from strict cross-study shared set).
"""
HARMONIZE = {
    # GSE307660 vocab
    "B cell": "B cell", "T cell": "T cell", "NK cell": "NK",
    "Monocyte/Myeloid": "Monocyte", "Dendritic cell": "Dendritic cell",
    "Erythroid": "Erythroid", "Malignant plasma cell": "Plasma cell",
    # GSE276974 vocab
    "CD4 T": "T cell", "CD8 T": "T cell", "NK": "NK", "Monocyte": "Monocyte",
    "Plasma cell": "Plasma cell", "HSPC": "Progenitor",
    "Granulo/Myeloid progenitor": "Progenitor",
    # pass-through
    "Unassigned": "Unassigned",
}
SHARED = ["B cell", "T cell", "NK", "Monocyte", "Dendritic cell", "Erythroid", "Plasma cell"]


def harmonize(series):
    """Map a pandas Series of raw cell_type labels to the harmonized vocabulary."""
    return series.map(lambda x: HARMONIZE.get(x, x))
