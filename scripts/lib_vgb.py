#!/usr/bin/env python3
"""Shared helpers for the Python stages. Imported, not run.

Holds the things more than one stage needs: reading project.conf, writing TSV
atomically, the crystallisation-additive list, receptor cleaning, and the
symmetry-corrected RMSD.

The RMSD section is the part worth reading. It is where a careless
implementation produces a flattering number, and it is commented accordingly.
"""

from __future__ import annotations

import csv
import os
import shutil
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# project.conf
# ---------------------------------------------------------------------------


def load_conf(path: str | Path) -> dict:
    """Read the shell project.conf into a dict of strings.

    project.conf is sourced by the bash stages, so it is shell syntax rather
    than ini or yaml. Only KEY=VALUE lines are taken, inline comments are cut,
    surrounding quotes are stripped. Nothing in it needs shell expansion, and
    keeping it as one file readable by both bash and python is worth more than
    the elegance of a second format.
    """
    conf: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        if not key.replace("_", "").isalnum():
            continue
        conf[key] = val.split("#", 1)[0].strip().strip('"').strip("'")
    return conf


def conf_int(conf: dict, key: str, default: int) -> int:
    try:
        return int(conf[key])
    except (KeyError, ValueError):
        return default


def conf_float(conf: dict, key: str, default: float) -> float:
    try:
        return float(conf[key])
    except (KeyError, ValueError):
        return default


# ---------------------------------------------------------------------------
# TSV
# ---------------------------------------------------------------------------


def read_tsv(path: str | Path) -> list[dict]:
    with open(path, newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if not (r.get(next(iter(r))) or "").startswith("#")]


def write_tsv(path: str | Path, columns: list[str], rows) -> None:
    """Write a TSV via a temporary file and one rename.

    A stage killed halfway must not leave a truncated results table that the
    next stage would read as complete.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                          lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    shutil.move(str(tmp), str(path))


def append_tsv(path: str | Path, columns: list[str], rows) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.is_file() and path.stat().st_size > 0
    with path.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                          lineterminator="\n", extrasaction="ignore")
        if not exists:
            w.writeheader()
        for r in rows:
            w.writerow(r)


# ---------------------------------------------------------------------------
# Receptor cleaning
# ---------------------------------------------------------------------------

# Components treated as crystallisation additives, cryoprotectants, buffers or
# ions rather than part of the binding site. They are stripped from the receptor
# before protonation.
#
# This list is a judgement call. Stripping something that is genuinely part of
# the site makes the pocket bigger than it is; keeping a glycerol that happens
# to sit in the pocket puts a spurious obstacle in it. The list is kept short
# and the count of atoms removed per complex is recorded, so the effect is
# visible rather than assumed.
CRYSTALLISATION_ADDITIVES = {
    # water and unassigned density
    "HOH", "DOD", "UNL", "UNX", "UNK",
    # polyols, cryoprotectants, detergents in small amounts
    "GOL", "EDO", "PEG", "PGE", "PG4", "PG0", "P6G", "1PE", "2PE", "MPD",
    "MRD", "BOG", "LDA", "C8E", "OCT", "HEZ", "PIN",
    # buffers and reductants
    "TRS", "BME", "DTT", "DTU", "MES", "EPE", "IMD", "BTB", "TAR", "MLI",
    "CIT", "FLC", "SIN", "ACT", "ACY", "FMT", "OXL", "MLA",
    # solvents
    "DMS", "DMF", "ACN", "MOH", "EOH", "IPA", "THJ",
    # monatomic ions and simple inorganics
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2",
    "CU", "CU1", "NI", "CO", "CD", "HG", "BR", "IOD", "F", "NO3", "CO3",
    "AZI", "CN", "SCN", "WO4", "MOO", "VO4", "BEF", "ALF", "NH4", "EDT",
    "PER", "SEK", "BR3", "YB", "SR", "CS", "RB", "LI", "BA", "PB", "AU",
}

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # protonation-state and modification variants pdb2pqr and others emit
    "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN", "TYM", "ARN",
    "MSE", "SEC", "PYL",
    # nucleic acids, for the handful of complexes with a nucleic chain
    "A", "C", "G", "U", "T", "DA", "DC", "DG", "DT", "DU", "I", "N",
}


def clean_receptor_pdb(src: Path, dst: Path, drop_ccd: str | None = None) -> dict:
    """Write a receptor PDB with waters, additives and stray ligand atoms gone.

    drop_ccd names the ligand of interest. Six of the 308 PoseBusters receptor
    files still contain one to four atoms of it: the authors' PyMOL selection
    left a fragment behind. One atom is not the answer sitting in the pocket,
    but it is a real atom in the grid where the ligand should be going, so it is
    removed and the count is recorded.

    Returns counts so the caller can put them in a table.
    """
    kept, n_water, n_additive, n_stray, n_cofactor_atoms = [], 0, 0, 0, 0
    cofactor_names: set[str] = set()
    drop = (drop_ccd or "").strip().upper()

    for line in src.read_text(errors="replace").splitlines():
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM"):
            resn = line[17:20].strip().upper()
            if rec == "HETATM":
                if resn == drop and drop:
                    n_stray += 1
                    continue
                if resn in ("HOH", "DOD"):
                    n_water += 1
                    continue
                if resn in CRYSTALLISATION_ADDITIVES:
                    n_additive += 1
                    continue
                if resn not in STANDARD_RESIDUES:
                    n_cofactor_atoms += 1
                    cofactor_names.add(resn)
            kept.append(line)
        # Everything else is dropped, including TER and END. This is not
        # tidiness. The archive files carry bare "TER" records with no serial
        # number or residue name, and once atoms have been deleted around them
        # pdb2pqr builds a residue with an empty atom list from one and dies
        # with an IndexError inside Residue.__init__. Six of six Astex receptors
        # failed that way before these two lines changed. REMARK and HEADER
        # records are dropped because pdb2pqr cannot carry them through anyway,
        # and CONECT records are dropped because ones pointing at deleted atoms
        # are worse than none.

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(kept) + "\nEND\n")
    return {
        "atoms_kept": len(kept),
        "waters_removed": n_water,
        "additive_atoms_removed": n_additive,
        "stray_ligand_atoms_removed": n_stray,
        "cofactor_atoms_kept": n_cofactor_atoms,
        "cofactors_kept": ",".join(sorted(cofactor_names)),
    }


# ---------------------------------------------------------------------------
# RMSD
# ---------------------------------------------------------------------------
#
# A naive atom-index RMSD is wrong for any ligand with a symmetric substructure.
# Rotate a para-substituted benzene by 180 degrees about its own axis and every
# atom lands on a chemically identical partner: the pose is the same pose, and an
# index-matched calculation reports two or three Angstrom of error. The same goes
# for a carboxylate's two oxygens, a tert-butyl's three methyls, a symmetric
# biaryl, a phosphate.
#
# The fix is to minimise over the automorphism group of the molecular graph.
# There is a second trap underneath that one, and it is worth spelling out
# because the PoseBusters paper's own wording sends you into it.
#
# RDKit has two symmetry-aware functions and they answer different questions.
#   CalcRMS      symmetry-corrected, computed in place. The molecules are left
#                where they are. This is the docking question: how far is the
#                pose from the crystal pose.
#   GetBestRMS   symmetry-corrected after an optimal superposition of the two
#                molecules. This is the conformer question: is the shape right.
#                A pose in completely the wrong subsite can score well here.
# The paper's methods say GetBestRMS. The PoseBusters code does not use it for
# the success threshold: posebusters/modules/rmsd.py computes both and applies
# rmsd_within_threshold to the CalcRMS value, reporting the GetBestRMS value
# separately as kabsch_rmsd. The in-place number is the right one and it is the
# one this pipeline calls rmsd.
#
# Three implementations are recorded per pose:
#
#   posebusters.robust_rmsd   the primary number. Calling their function rather
#                             than CalcRMS directly is deliberate: robust_rmsd
#                             retries with charges stripped, with tautomers
#                             canonicalised and with bond orders reassigned from
#                             the reference when the first attempt fails to find
#                             a substructure match. Docked poses read back from
#                             PDBQT frequently need that, and a home-rolled
#                             CalcRMS call would return NaN where theirs returns
#                             a number, which would quietly drop complexes.
#   spyrmsd                   an independent graph-isomorphism implementation.
#                             Not a tie-breaker, because there is nothing to
#                             break a tie against. It is there so a disagreement
#                             between two libraries appears as a column in the
#                             results table rather than never being looked for.
#   naive index RMSD          no symmetry correction, recorded only so the
#                             README can state how often correction changes the
#                             verdict at 2 Angstrom.
#
# All on heavy atoms only. Hydrogen positions in a crystal structure are placed
# by the refinement program rather than observed, so including them would
# measure the refinement program's guess.


def load_heavy(path: str | Path, sanitize: bool = True):
    """Read the first molecule from an SDF, hydrogens removed."""
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(path), removeHs=True, sanitize=sanitize)
    for mol in supplier:
        if mol is not None:
            return mol
    return None


def symm_rmsd_posebusters(probe, ref) -> tuple[float | None, float | None]:
    """Return (rmsd, kabsch_rmsd) using PoseBusters' own implementation.

    rmsd is in place and is the success metric. kabsch_rmsd superposes first and
    is reported alongside, because a pose with a low kabsch and a high in-place
    RMSD has the right conformation in the wrong place, which is a search
    failure rather than a conformer failure.
    """
    from rdkit import Chem

    try:
        from posebusters.modules.rmsd import robust_rmsd
    except Exception:
        return None, None

    def one(kabsch: bool):
        try:
            import math

            v = robust_rmsd(Chem.Mol(ref), Chem.Mol(probe), heavy_only=True,
                            kabsch=kabsch)
            v = float(v)
            return None if math.isnan(v) else v
        except Exception:
            return None

    return one(False), one(True)


def symm_rmsd_rdkit(probe, ref) -> float | None:
    """Symmetry-corrected heavy-atom RMSD via RDKit CalcRMS, computed in place.

    Kept as the plain-RDKit reading of the same quantity. The copies are made
    because the RDKit alignment functions are documented to be free to modify
    their arguments.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    try:
        return float(rdMolAlign.CalcRMS(Chem.Mol(probe), Chem.Mol(ref)))
    except Exception:
        return None


def symm_rmsd_spyrmsd(probe, ref) -> float | None:
    """Symmetry-corrected heavy-atom RMSD, spyrmsd, graph isomorphism."""
    try:
        import numpy as np
        from spyrmsd import rmsd as spy_rmsd
        from spyrmsd.molecule import Molecule
    except Exception:
        return None
    try:
        pm = Molecule.from_rdkit(probe)
        rm = Molecule.from_rdkit(ref)
        pm.strip()
        rm.strip()
        val = spy_rmsd.symmrmsd(
            pm.coordinates, rm.coordinates,
            pm.atomicnums, rm.atomicnums,
            pm.adjacency_matrix, rm.adjacency_matrix,
            minimize=False,
        )
        return float(np.asarray(val).ravel()[0])
    except Exception:
        return None


def first_match_rmsd(probe, ref) -> float | None:
    """RMSD using one arbitrary valid atom correspondence instead of the best.

    This is the number the symmetry-corrected RMSD is compared against, and
    getting the comparison right took two attempts.

    The first attempt matched atoms by index. That is not a symmetry question at
    all: a pose that has been through a PDBQT file comes back in a different
    atom order from the crystal SDF, so index matching compares a carbon with a
    chlorine and reports eight Angstrom on a pose that is actually within one.
    It measures file format handling, not symmetry.

    What isolates symmetry is to take a single valid graph match, the first one
    RDKit returns, and use that correspondence. The symmetry-corrected value
    minimises over every match; this one takes whichever came out first. The
    difference between the two is the cost of not correcting for symmetry, and
    for a molecule with no symmetric substructure there is exactly one match and
    the two numbers are identical.

    Never used as a result. It exists so the README can state how often
    symmetry correction changes the verdict at the 2 Angstrom line.
    """
    try:
        import numpy as np

        match = probe.GetSubstructMatch(ref)
        if not match or len(match) != ref.GetNumAtoms():
            return None
        a = probe.GetConformer().GetPositions()
        b = ref.GetConformer().GetPositions()
        diff = a[list(match)] - b
        return float(np.sqrt((diff ** 2).sum(axis=1).mean()))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# misc
# ---------------------------------------------------------------------------


def eprint(*args) -> None:
    print(*args, file=sys.stderr, flush=True)


def which_or_die(name: str, hint: str = "") -> str:
    p = shutil.which(name)
    if not p:
        eprint(f"[error] {name} not found on PATH. {hint}")
        sys.exit(3)
    return p


def du_kb(path: str | Path) -> int:
    """Apparent size of a directory tree in kilobytes."""
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total // 1024
