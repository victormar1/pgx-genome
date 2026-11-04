# -*- coding: utf-8 -*-
"""Table de verite GeT-RM, depuis le classeur consolide du CDC.

Tous les genes nommes du classeur, y compris CYP2B6, DPYD, NAT2, UGT1A1, VKORC1
et les deux colonnes d'alleles HLA.
Source : CDC GeT-RM, Consolidated_PGx-HLA_table, version du 22 janvier 2025.

Usage :
    construire_verite_getrm.py --classeur Consolidated_PGx-HLA_table_1-22-25-V4.xlsx \
        --sortie verites/verite_getrm.tsv
Demande openpyxl.
"""
import argparse, csv, re, collections
from openpyxl import load_workbook

_p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
_p.add_argument("--classeur", required=True)
_p.add_argument("--sortie", required=True)
_a = _p.parse_args()
SRC, OUT = _a.classeur, _a.sortie

# nom court -> nom PharmCAT
def nom(h):
    h = re.sub(r"\s+", " ", h).strip()
    for g in ("CYP1A1","CYP1A2","CYP2A6","CYP2B6","CYP2C8","CYP2C9","CYP2C19","CYP2D6","CYP2E1",
              "CYP3A4","CYP3A5","CYP4F2","DPYD","GSTM1","GSTP1","GSTT1","NAT1","NAT2","NUDT15",
              "SLC15A2","SLC22A2","SLCO1B1","SLCO2B1","TPMT","UGT1A1","UGT2B7","UGT2B15","UGT2B17",
              "HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5",
              "HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1"):
        if h.startswith(g): return g
    if h.startswith("CYP2C Cluster"): return "CYP2C_cluster"
    if h.startswith("VKORC1 NM_024006.5"): return "VKORC1"
    if h.startswith("GGCX"): return None
    return None

wb = load_workbook(SRC, read_only=True, data_only=True)
ws = wb["PGx HLA Genotypes"]
rows = list(ws.iter_rows(values_only=True))
wb.close()
hdr = [str(x).strip() if x is not None else "" for x in rows[1]]
skip = re.compile(r"Reference|URL|accession|ftp|Coriell|Characteriz", re.I)
col_ech = next(i for i, h in enumerate(hdr) if h.startswith("Coriell"))

cols = collections.defaultdict(list)
for i, h in enumerate(hdr):
    if not h or skip.search(h): continue
    g = nom(h)
    if g: cols[g].append(i)

genes = sorted(cols)
lignes = []
for r in rows[2:]:
    s = r[col_ech]
    if not s: continue
    s = str(s).strip().split()[0]
    if not re.match(r"^(HG|NA)\d+$", s): continue
    d = {"sample": s}
    for g in genes:
        vals = [str(r[i]).strip() for i in cols[g] if i < len(r) and r[i] not in (None, "")]
        vals = [v for v in vals if v.lower() not in ("n/a", "na", "nd", "none", "-")]
        if not vals: d[g] = ""
        elif g.startswith("HLA") and len(vals) == 2: d[g] = f"{vals[0]}/{vals[1]}"
        else: d[g] = vals[0]
    lignes.append(d)

with open(OUT, "w", newline="\n", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["sample"] + genes, delimiter="\t")
    w.writeheader(); w.writerows(lignes)

print(f"{len(lignes)} echantillons GeT-RM, {len(genes)} genes")
for g in genes:
    n = sum(1 for d in lignes if d[g])
    if n: print(f"{g:16s} {n:7d}")
print("ecrit :", OUT)
