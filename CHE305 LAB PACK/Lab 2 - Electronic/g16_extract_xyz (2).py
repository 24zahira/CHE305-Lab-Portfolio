
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract the final (optimized) geometry from a Gaussian .log/.out file to XYZ.

Usage:
    python g16_extract_xyz.py input.log -o output.xyz

Notes:
- Prefers the last 'Standard orientation' block after 'Optimization completed.'
- Falls back to the last available orientation block if no completion marker is found.
- Handles 'Standard orientation', 'Input orientation', and 'Z-Matrix orientation' blocks.
"""

import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple

# Periodic table mapping: atomic number -> symbol
PT = {
    1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F", 10:"Ne",
    11:"Na",12:"Mg",13:"Al",14:"Si",15:"P",16:"S",17:"Cl",18:"Ar",19:"K",20:"Ca",
    21:"Sc",22:"Ti",23:"V",24:"Cr",25:"Mn",26:"Fe",27:"Co",28:"Ni",29:"Cu",30:"Zn",
    31:"Ga",32:"Ge",33:"As",34:"Se",35:"Br",36:"Kr",37:"Rb",38:"Sr",39:"Y",40:"Zr",
    41:"Nb",42:"Mo",43:"Tc",44:"Ru",45:"Rh",46:"Pd",47:"Ag",48:"Cd",49:"In",50:"Sn",
    51:"Sb",52:"Te",53:"I",54:"Xe",55:"Cs",56:"Ba",57:"La",58:"Ce",59:"Pr",60:"Nd",
    61:"Pm",62:"Sm",63:"Eu",64:"Gd",65:"Tb",66:"Dy",67:"Ho",68:"Er",69:"Tm",70:"Yb",
    71:"Lu",72:"Hf",73:"Ta",74:"W",75:"Re",76:"Os",77:"Ir",78:"Pt",79:"Au",80:"Hg",
    81:"Tl",82:"Pb",83:"Bi",84:"Po",85:"At",86:"Rn",87:"Fr",88:"Ra",89:"Ac",90:"Th",
    91:"Pa",92:"U",93:"Np",94:"Pu",95:"Am",96:"Cm",97:"Bk",98:"Cf",99:"Es",100:"Fm",
    101:"Md",102:"No",103:"Lr",104:"Rf",105:"Db",106:"Sg",107:"Bh",108:"Hs",109:"Mt",
    110:"Ds",111:"Rg",112:"Cn",113:"Nh",114:"Fl",115:"Mc",116:"Lv",117:"Ts",118:"Og"
}

ORIENTATION_HEADERS = (
    "Standard orientation:",
    "Input orientation:",
    "Z-Matrix orientation:"
)

def read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        return f.read()

def find_optimization_completed_indices(text: str) -> List[int]:
    """
    Find indices where Gaussian reports successful optimization.
    Looks for common markers like 'Optimization completed'.
    """
    indices = [m.start() for m in re.finditer(r"Optimization completed\.?", text)]
    # Also include 'Stationary point found' for TS/TS-like workflows
    indices += [m.start() for m in re.finditer(r"Stationary point found", text)]
    indices.sort()
    return indices

def iter_orientation_blocks(text: str) -> List[Tuple[str, int, int, List[str]]]:
    """
    Return a list of geometry blocks found. Each entry:
    (header_name, start_index, end_index, lines_of_block)
    The lines_of_block are the table rows only (no header lines).
    """
    blocks = []
    # Pattern: header line, dashed line, column header line, dashed line, then table, ending at dashed line
    # Example lines after header:
    #  ---------------------------------------------------------------------
    #  Center     Atomic      Atomic             Coordinates (Angstroms)
    #  Number     Number       Type             X           Y           Z
    #  ---------------------------------------------------------------------
    #       1          6           0        -0.000000    0.000000    0.000000
    #  ...
    #  ---------------------------------------------------------------------
    header_regex = re.compile(rf"^(\s*({'|'.join(map(re.escape, ORIENTATION_HEADERS))}))\s*$", re.MULTILINE)

    for m in header_regex.finditer(text):
        header = m.group(2)
        start = m.end()
        # Find the dashed separators
        dashed = list(re.finditer(r"^\s*-{5,}\s*$", text[start:], re.MULTILINE))
        if len(dashed) < 2:
            # Not a well-formed block; skip
            continue
        # The first dashed line is right after header; the second dashed line after column names; the third ends the table
        # We need to be careful: typical structure has two dashed lines before table and one after.
        # Find the third dashed line relative to 'start'
        table_start_rel = None
        table_end_rel = None

        # Locate the first two dashed lines for headers
        if len(dashed) >= 2:
            # After the second dashed line, table starts
            table_start_rel = dashed[1].end()
        # Now find the next dashed line after table_start_rel to bound the table
        if table_start_rel is not None:
            next_dashed_after_table = re.search(r"^\s*-{5,}\s*$", text[start + table_start_rel:], re.MULTILINE)
            if next_dashed_after_table:
                table_end_rel = table_start_rel + next_dashed_after_table.start()

        if table_start_rel is None or table_end_rel is None:
            continue

        table_text = text[start + table_start_rel : start + table_end_rel]
        lines = [ln.rstrip() for ln in table_text.splitlines() if ln.strip()]
        blocks.append((header, m.start(), start + table_end_rel, lines))

    return blocks

def parse_block_lines(lines: List[str]) -> List[Tuple[str, float, float, float]]:
    """
    Parse lines of an orientation block.
    Expected columns: Center Number, Atomic Number, Atomic Type, X, Y, Z
    Returns list of tuples: (symbol, x, y, z)
    """
    geom = []
    for ln in lines:
        # Split respecting possible multiple spaces
        parts = ln.split()
        if len(parts) < 6:
            continue
        try:
            atno = int(parts[1])
            x, y, z = map(float, parts[-3:])
        except ValueError:
            continue
        symbol = PT.get(atno, "X")
        geom.append((symbol, x, y, z))
    return geom

def pick_final_geometry(text: str) -> Tuple[List[Tuple[str, float, float, float]], str]:
    """
    Choose the best geometry:
    1) If there is a completion marker, pick the last geometry block AFTER that marker (prefer 'Standard orientation').
    2) Else, pick the last geometry block in the file (prefer 'Standard orientation').
    Returns (geometry, status_comment).
    """
    blocks = iter_orientation_blocks(text)
    if not blocks:
        raise RuntimeError("No orientation blocks found in the output file.")

    completion_points = find_optimization_completed_indices(text)
    preferred = "Standard orientation:"

    def filter_blocks_after(idx: int) -> List[Tuple[str, int, int, List[str]]]:
        return [b for b in blocks if b[1] > idx]

    chosen = None
    status = "Final geometry (last available block)."

    if completion_points:
        # Take the last completion marker
        last_done = completion_points[-1]
        candidates = filter_blocks_after(last_done)
        if candidates:
            # Prefer Standard orientation, else the last one
            std = [b for b in candidates if b[0] == preferred]
            chosen = std[-1] if std else candidates[-1]
            status = "Optimized geometry (Optimization completed)."

    if chosen is None:
        # No completion or no blocks after completion; take the last block, prefer 'Standard orientation'
        std = [b for b in blocks if b[0] == preferred]
        chosen = std[-1] if std else blocks[-1]

    geom = parse_block_lines(chosen[3])
    if not geom:
        raise RuntimeError("Found an orientation block but failed to parse coordinates.")
    return geom, status

def write_xyz(geom: List[Tuple[str, float, float, float]], out_path: Path, comment: str):
    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"{len(geom)}\n")
        f.write(comment.strip() + "\n")
        for sym, x, y, z in geom:
            f.write(f"{sym:2s}  {x: .8f}  {y: .8f}  {z: .8f}\n")

def main():
    ap = argparse.ArgumentParser(description="Extract final/optimized geometry from Gaussian log to XYZ.")
    ap.add_argument("input", help="Gaussian .log/.out file")
    ap.add_argument("-o", "--output", help="Output .xyz path (default: input basename + _opt.xyz)")
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"Input file not found: {in_path}")

    text = read_text(in_path)
    geom, status = pick_final_geometry(text)

    out_path = Path(args.output) if args.output else in_path.with_suffix("").with_name(in_path.stem + "_opt.xyz")
    comment = f"{in_path.name} — {status}"
    write_xyz(geom, out_path, comment)

    print(f"Wrote {out_path} ({len(geom)} atoms)\nComment: {comment}")

if __name__ == "__main__":
    main()

