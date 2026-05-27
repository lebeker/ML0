"""
Convert PubChem XML (PC-Compounds) to CSV — streaming via iterparse.
Processes one PC-Compound at a time, memory usage stays flat regardless of file size.

Usage:
    python to_csv.py [input.xml] [output.csv]

Defaults:
    input  = Compound_000000001_000500000.xml
    output = molecules.csv
"""

import csv
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

INPUT  = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("Compound_000000001_000500000.xml")
OUTPUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("molecules.csv")

NS = "http://www.ncbi.nlm.nih.gov"

# Column order in the output CSV.
# (label, name) → csv column name.  name=None means the property has no sub-name.
PROP_MAP = {
    ("Compound",             "Canonicalized")       : "canonicalized",
    ("Compound Complexity",  None)                  : "complexity",
    ("Count",                "Hydrogen Bond Acceptor"): "hb_acceptors",
    ("Count",                "Hydrogen Bond Donor")  : "hb_donors",
    ("Count",                "Rotatable Bond")       : "rotatable_bonds",
    ("Fingerprint",          "SubStructure Keys")    : "fingerprint",
    ("IUPAC Name",           "Preferred")            : "iupac_name_preferred",
    ("IUPAC Name",           "Allowed")              : "iupac_name_allowed",
    ("IUPAC Name",           "CAS-like Style")       : "iupac_name_cas_style",
    ("IUPAC Name",           "Markup")               : "iupac_name_markup",
    ("IUPAC Name",           "Systematic")           : "iupac_name_systematic",
    ("IUPAC Name",           "Traditional")          : "iupac_name_traditional",
    ("InChI",                "Standard")             : "inchi",
    ("InChIKey",             "Standard")             : "inchikey",
    ("Log P",                "XLogP3-AA")            : "xlogp3",
    ("Mass",                 "Exact")                : "exact_mass",
    ("Molecular Formula",    None)                   : "molecular_formula",
    ("Molecular Weight",     None)                   : "molecular_weight",
    ("SMILES",               "Absolute")             : "smiles_isomeric",
    ("SMILES",               "Connectivity")         : "smiles_canonical",
    ("Topological",          "Polar Surface Area")   : "tpsa",
    ("Weight",               "MonoIsotopic")         : "monoisotopic_weight",
}

FIELDNAMES = ["cid"] + list(PROP_MAP.values())


def tag(local):
    return f"{{{NS}}}{local}"


def _text(elem, local):
    child = elem.find(tag(local))
    return child.text.strip() if child is not None and child.text else None


def extract_compound(elem):
    row = {col: "" for col in FIELDNAMES}

    # CID
    cid_el = elem.find(f".//{tag('PC-CompoundType_id_cid')}")
    row["cid"] = cid_el.text.strip() if cid_el is not None and cid_el.text else ""

    # Properties
    for info in elem.iter(tag("PC-InfoData")):
        urn = info.find(tag("PC-InfoData_urn"))
        if urn is None:
            continue

        urn_inner = urn.find(tag("PC-Urn"))
        if urn_inner is None:
            continue
        label_el = urn_inner.find(tag("PC-Urn_label"))
        name_el  = urn_inner.find(tag("PC-Urn_name"))
        label = label_el.text.strip() if label_el is not None and label_el.text else None
        name  = name_el.text.strip()  if name_el  is not None and name_el.text  else None

        col = PROP_MAP.get((label, name)) or PROP_MAP.get((label, None))
        if col is None:
            continue

        val_block = info.find(tag("PC-InfoData_value"))
        if val_block is None:
            continue

        # Pick whichever value child is present
        for val_tag in ("PC-InfoData_value_sval", "PC-InfoData_value_fval",
                        "PC-InfoData_value_ival", "PC-InfoData_value_binary"):
            child = val_block.find(tag(val_tag))
            if child is not None and child.text:
                row[col] = child.text.strip()
                break

    return row


LIMIT = None  # set to an int to stop early for debugging


def convert(input_path: Path, output_path: Path):
    print(f"Input : {input_path}  ({input_path.stat().st_size / 1e9:.2f} GB)")
    print(f"Output: {output_path}")

    written = 0
    root_ref = None

    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()

        context = ET.iterparse(str(input_path), events=("start", "end"))

        for event, elem in context:
            local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag

            if event == "start" and local == "PC-Compounds":
                root_ref = elem

            elif event == "end" and local == "PC-Compound":
                row = extract_compound(elem)
                writer.writerow(row)
                written += 1
                if written <= 5:
                    print(f"  [{written}] cid={row['cid']}  iupac={row.get('iupac_name_preferred')}  formula={row.get('molecular_formula')}")

                if root_ref is not None:
                    root_ref.remove(elem)
                else:
                    elem.clear()

                if LIMIT and written >= LIMIT:
                    print(f"Stopped at {LIMIT} records (debug mode).")
                    break

                if written % 10_000 == 0:
                    print(f"  {written:,} compounds written…")

    print(f"Done. {written:,} compounds → {output_path}")


if __name__ == "__main__":
    convert(INPUT, OUTPUT)
