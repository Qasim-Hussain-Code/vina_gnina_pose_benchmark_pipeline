#!/usr/bin/env python3
"""Prepare receptors and ligands for every arm.

Receptors
---------
The archive's protein file is stripped of waters, of crystallisation additives,
and of any stray atoms of the ligand of interest. Six of the 308 PoseBusters
receptor files contain one to four atoms of their own ligand, left behind by the
authors' PyMOL selection; one atom is not the answer sitting in the pocket, but
it is a real atom in the grid where the ligand is meant to go, so it comes out
and the count goes in the table.

Protonation is pdb2pqr with PROPKA at pH 7.4, AMBER naming. The reason for
pdb2pqr rather than reduce is that PROPKA predicts a pKa per titratable group
and sets histidine tautomers accordingly, where reduce optimises the hydrogen
network without changing formal states. Neither is right; pdb2pqr at least
records what it decided.

Cofactors are kept. This matters more than it sounds: 105 of the 308 PoseBusters
receptors and 29 of the 85 Astex receptors carry a non-additive cofactor, most
often HEM, FAD, NAD(P), FMN or an iron-sulfur cluster. pdb2pqr's PQR output
drops all of them because AMBER has no parameters for them, and a heme site
docked without its heme is a large empty cavity that a ligand will happily fall
into. The route used here takes pdb2pqr's --pdb-output instead of its PQR, which
keeps the cofactor records, and types the result with Open Babel. Meeko would
type the protein better and cannot type HEM at all; the note in
prepare_receptor has the argument and what the choice costs.

Ligands
-------
Two versions per complex, because one of the four questions this repository asks
is how much of the quoted success rate comes from the starting conformer.

  refconf  the crystal ligand's own coordinates, hydrogens added in place.
           The search starts from the answer. This is the flattering arm and it
           is reported so the gap can be measured.
  genconf  the same molecule rebuilt from its SMILES with ETKDGv3 and an MMFF
           relaxation, seeded from project.conf. No alignment to the crystal
           pose at any point: aligning it would put the information back.

The SMILES is taken from the crystal SDF rather than from the PDB chemical
component dictionary, so the generated molecule is guaranteed to be the same
molecule as the reference, with the same formal charges and the same formula.
A mismatch there would make the RMSD undefined and the PoseBusters formula check
fail for reasons that have nothing to do with docking.

Cross-docking receptors
-----------------------
The partner structure is fetched from the RCSB, cleaned the same way, and has
its own ligand removed. The reference pose for RMSD is the query ligand's
crystal pose carried into the partner's frame: the two proteins are superposed
on matched C-alpha atoms and the transform is applied to the ligand. That is
where the ligand sits when this protein binds it, expressed in the coordinate
system of the receptor being docked into.

Usage
-----
    python scripts/04_prepare.py --config project.conf
    python scripts/04_prepare.py --config project.conf --dataset astex --jobs 8
    python scripts/04_prepare.py --config project.conf --what crossdock --force
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import io
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_vgb as L  # noqa: E402

RCSB_CIF = "https://files.rcsb.org/download/{pdb}.cif.gz"

RECEPTOR_COLUMNS = [
    "scope", "key", "source_pdb", "status", "reason", "atoms_kept",
    "waters_removed", "additive_atoms_removed", "stray_ligand_atoms_removed",
    "cofactor_atoms_kept", "cofactors_kept", "protonation", "typing",
    "typing_note", "pdbqt_atoms", "pdbqt_bytes", "elapsed_s", "recorded",
]
LIGAND_COLUMNS = [
    "dataset", "complex_id", "conformer", "status", "reason", "smiles",
    "heavy_atoms", "rotatable_bonds", "seed", "embed_attempts",
    "conformer_rmsd_to_crystal", "pdbqt_bytes", "elapsed_s", "recorded",
]
ALIGN_COLUMNS = [
    "pair_id", "query_pdb_id", "receptor_pdb_id", "status", "reason",
    "chains_matched", "ca_pairs", "ca_rmsd_after_superposition",
    "sequence_identity_of_matched", "partner_ligand_atoms",
    "partner_ligand_center_x", "partner_ligand_center_y",
    "partner_ligand_center_z", "recorded",
]


# ---------------------------------------------------------------------------
# receptor
# ---------------------------------------------------------------------------
def prepare_receptor(src_pdb: Path, out_dir: Path, drop_ccd: str | None,
                     tools: dict, force: bool) -> dict:
    """Clean, protonate and type one receptor. Returns a table row."""
    t0 = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    pdbqt = out_dir / "receptor.pdbqt"
    row = {"status": "ok", "reason": "", "protonation": "pdb2pqr_propka_ph7.4",
           "typing": "openbabel_gasteiger_rigid",
           "recorded": datetime.datetime.now().astimezone().isoformat(timespec="seconds")}

    if pdbqt.is_file() and pdbqt.stat().st_size > 0 and not force:
        row.update({"status": "cached", "pdbqt_bytes": pdbqt.stat().st_size,
                    "elapsed_s": 0.0})
        return row

    clean = out_dir / "clean.pdb"
    try:
        counts = L.clean_receptor_pdb(src_pdb, clean, drop_ccd)
    except Exception as exc:
        row.update({"status": "failed", "reason": f"cleaning failed: {exc!r}"[:200],
                    "elapsed_s": round(time.time() - t0, 2)})
        return row
    row.update(counts)
    if counts["atoms_kept"] < 100:
        row.update({"status": "failed",
                    "reason": f"only {counts['atoms_kept']} atoms survived cleaning",
                    "elapsed_s": round(time.time() - t0, 2)})
        return row

    # --- protonation ------------------------------------------------------
    prot = out_dir / "protonated.pdb"
    pqr = out_dir / "receptor.pqr"
    cmd = [tools["pdb2pqr"], "--ff=AMBER", "--with-ph=7.4",
           "--titration-state-method=propka", "--keep-chain", "--drop-water",
           f"--pdb-output={prot}", str(clean), str(pqr)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0 or not prot.is_file() or prot.stat().st_size == 0:
        # Open Babel as the fallback protonator. It does not predict pKa values,
        # so a complex that lands here has a worse protonation than the rest and
        # the table says so rather than hiding it.
        prot = out_dir / "protonated_obabel.pdb"
        p2 = subprocess.run([tools["obabel"], str(clean), "-O", str(prot), "-p", "7.4"],
                            capture_output=True, text=True, timeout=1800)
        if p2.returncode != 0 or not prot.is_file() or prot.stat().st_size == 0:
            tail = (p.stderr or p.stdout or "").strip().splitlines()[-1:] or [""]
            row.update({"status": "failed",
                        "reason": f"pdb2pqr and obabel both failed: {tail[0]}"[:200],
                        "elapsed_s": round(time.time() - t0, 2)})
            return row
        row["protonation"] = "obabel_ph7.4_fallback"

    # --- typing and charges ------------------------------------------------
    # Open Babel, not Meeko, and the reason is cofactors.
    #
    # Meeko's mk_prepare_receptor is the better tool for the polymer: it types
    # residues from curated templates rather than perceiving them. It cannot
    # type a cofactor it has no template for. Asked to handle HEM it tries to
    # build a template from the chemical component dictionary, gets None back,
    # and raises AttributeError inside chemtempgen.export_chem_templates_to_json.
    # HEM appears in 12 PoseBusters receptors and 5 Astex ones, and the same
    # path fails for other metal-containing components.
    #
    # The options were: delete the cofactors Meeko cannot type, hand-write
    # templates for the forty-odd cofactors in this set, or type the whole
    # receptor with Open Babel. Deleting them is wrong, because a heme site
    # without its heme is a large empty cavity that a ligand falls into and
    # scores well in. Hand-writing templates is a week of work and would still
    # leave the next dataset broken. Open Babel types every receptor the same
    # way, which also removes the confound of one protocol for the two thirds
    # without cofactors and another for the third with them.
    #
    # What this costs: Open Babel perceives atom types and assigns Gasteiger
    # charges rather than reading them from residue templates, so the protein
    # typing is cruder than Meeko's. For Vina and Vinardo that matters less than
    # it sounds, because both scoring functions distinguish only a handful of
    # receptor atom types plus hydrogen-bond donor and acceptor flags. It would
    # matter more for a force-field score. This is a fixed choice, not an
    # explored one, and the Limitations section says so.
    #
    # -xr writes a rigid receptor with no torsion tree. Non-polar hydrogens are
    # merged into their carbons, which is the united-atom convention Vina
    # expects; the polar hydrogens pdb2pqr placed are kept.
    #
    # --partialcharge gasteiger is not optional and its absence is silent. Open
    # Babel writes a PDBQT with a charge column of +0.000 for every atom unless
    # asked for a charge model, and nothing complains: the file parses, Vina
    # reads it, and the docking runs. The first receptors prepared here had
    # 11,465 atoms of zero charge each. It happens not to change the Vina or
    # Vinardo result, because neither scoring function has an electrostatic term,
    # but GNINA's CNN sees the receptor as typed atoms with charges and an
    # AutoDock4 rescoring would be nonsense. A zero-charge receptor is wrong even
    # where it is harmless.
    def obabel_pdbqt(extra: list[str]):
        p = subprocess.run([tools["obabel"], str(prot), "-O", str(pdbqt), "-xr"] + extra,
                           capture_output=True, text=True, timeout=3600)
        ok = p.returncode == 0 and pdbqt.is_file() and pdbqt.stat().st_size > 0
        return ok, p

    ok, p3 = obabel_pdbqt(["--partialcharge", "gasteiger"])
    if not ok:
        # Gasteiger assignment needs the molecule kekulized, and Open Babel
        # cannot kekulize 27 of the 308 PoseBusters receptors: the ones with a
        # heme or a similar aromatic metal-containing cofactor. Without the
        # charge model the same file converts with only a warning, so asking for
        # charges turns a warning into a fatal error and loses the complex.
        #
        # Dropping those 27 would be a silent exclusion of exactly the receptors
        # with the most interesting cofactors. Instead the conversion is retried
        # without a charge model and the typing column records which of the two
        # ran, so the difference is a column rather than a missing row. Neither
        # Vina nor Vinardo has an electrostatic term, so for those two the
        # retried receptors are identical to the rest.
        ok, p3 = obabel_pdbqt([])
        row["typing"] = "openbabel_nocharge_rigid"
        row["typing_note"] = ("gasteiger failed, kekulization: "
                              + ((p3.stderr or p3.stdout or "").strip().splitlines() or [""])[-1])[:160]
    if not ok:
        tail = (p3.stderr or p3.stdout or "").strip().splitlines()[-1:] or [""]
        row.update({"status": "failed",
                    "reason": f"obabel receptor conversion failed: {tail[0]}"[:200],
                    "elapsed_s": round(time.time() - t0, 2)})
        return row

    # Header records Vina cannot read. When pdb2pqr fails and Open Babel does
    # the protonation, the intermediate PDB carries Open Babel's own COMPND
    # (the input file path) and AUTHOR records, and the typing step copies them
    # into the PDBQT. Vina 1.2.7 stops at the COMPND line with a parse error;
    # smina skips it. Every receptor that took the fallback therefore failed
    # Vina in every arm while docking normally under Vinardo: 7 Astex, 5
    # PoseBusters and 6 cross-docking receptors in the first full run, all
    # logged in results/excluded.tsv with the COMPND line as the reason. The
    # records carry nothing the docking needs, so they are dropped here.
    lines = pdbqt.read_text().splitlines()
    kept = [ln for ln in lines if not ln.startswith(("COMPND", "AUTHOR"))]
    if len(kept) != len(lines):
        pdbqt.write_text("\n".join(kept) + "\n")
        row["header_records_dropped"] = len(lines) - len(kept)

    n_rec_atoms = sum(1 for line in pdbqt.read_text().splitlines()
                      if line.startswith("ATOM") or line.startswith("HETATM"))
    row["pdbqt_atoms"] = n_rec_atoms
    if n_rec_atoms < 100:
        row.update({"status": "failed",
                    "reason": f"receptor PDBQT has only {n_rec_atoms} atoms",
                    "elapsed_s": round(time.time() - t0, 2)})
        return row

    # The protonated PDB and the PQR go now rather than at the end of the stage:
    # across 643 receptors they are the largest thing written here and nothing
    # downstream reads them. clean.pdb stays, because 05_define_boxes.py runs
    # fpocket on it and fpocket wants a PDB rather than a PDBQT.
    for f in (pqr, out_dir / "protonated.pdb", out_dir / "protonated_obabel.pdb",
              out_dir / "receptor.log"):
        try:
            if f.is_file():
                f.unlink()
        except OSError:
            pass

    row.update({"pdbqt_bytes": pdbqt.stat().st_size,
                "elapsed_s": round(time.time() - t0, 2)})
    return row


# ---------------------------------------------------------------------------
# ligand
# ---------------------------------------------------------------------------
def _write_pdbqt(molh, path: Path) -> tuple[bool, str]:
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    prep = MoleculePreparation()
    setups = prep.prepare(molh)
    if not setups:
        return False, "Meeko returned no setup"
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return False, str(err)[:200]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(s)
    return True, ""


def prepare_ligand(complex_id: str, ligand_sdf: Path, out_dir: Path,
                   conformer: str, seed: int, force: bool) -> dict:
    """Prepare one ligand in one conformer flavour. Returns a table row."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdMolAlign

    RDLogger.DisableLog("rdApp.*")
    t0 = time.time()
    pdbqt = out_dir / f"{complex_id}.pdbqt"
    row = {"dataset": "", "complex_id": complex_id, "conformer": conformer,
           "status": "ok", "reason": "", "seed": seed if conformer == "genconf" else "",
           "embed_attempts": "", "conformer_rmsd_to_crystal": "",
           "recorded": datetime.datetime.now().astimezone().isoformat(timespec="seconds")}

    if pdbqt.is_file() and pdbqt.stat().st_size > 0 and not force:
        row.update({"status": "cached", "pdbqt_bytes": pdbqt.stat().st_size,
                    "elapsed_s": 0.0})
        return row

    ref = L.load_heavy(ligand_sdf)
    if ref is None:
        row.update({"status": "failed", "reason": "crystal SDF would not parse",
                    "elapsed_s": round(time.time() - t0, 2)})
        return row

    smiles = Chem.MolToSmiles(ref)
    row.update({"smiles": smiles, "heavy_atoms": ref.GetNumAtoms(),
                "rotatable_bonds": Descriptors.NumRotatableBonds(ref)})

    if conformer == "refconf":
        # Crystal coordinates, hydrogens added onto them. addCoords=True places
        # the hydrogens geometrically; their positions are not observed in the
        # crystal either way, and the RMSD is heavy-atom only.
        molh = Chem.AddHs(ref, addCoords=True)
    else:
        # Rebuilt from SMILES. ETKDGv3 with the seed from project.conf, then an
        # MMFF relaxation, falling back to UFF for the elements MMFF has no
        # parameters for. useRandomCoords is the documented retry for molecules
        # where distance geometry fails from the default start.
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if mol is None:
            row.update({"status": "failed", "reason": "SMILES would not parse back",
                        "elapsed_s": round(time.time() - t0, 2)})
            return row
        cid, attempts = -1, 0
        for attempt in range(3):
            attempts = attempt + 1
            ps = AllChem.ETKDGv3()
            ps.randomSeed = seed + attempt
            ps.useRandomCoords = attempt > 0
            ps.maxIterations = 1000 * (attempt + 1)
            cid = AllChem.EmbedMolecule(mol, ps)
            if cid >= 0:
                break
        row["embed_attempts"] = attempts
        if cid < 0:
            row.update({"status": "failed",
                        "reason": "ETKDGv3 failed to embed after 3 attempts",
                        "elapsed_s": round(time.time() - t0, 2)})
            return row
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            ff = AllChem.MMFFGetMoleculeForceField(mol, props) if props else None
            if ff is None:
                ff = AllChem.UFFGetMoleculeForceField(mol)
            if ff is not None:
                ff.Minimize(maxIts=2000)
        except Exception:
            pass    # an unrelaxed conformer is still a conformer; the internal
                    # RMSD column will show it is a poor one
        molh = mol
        try:
            # How far the generated conformer is from the crystal conformer as a
            # shape, after optimal superposition. This is the conformer
            # question, not the docking question, and it sets a floor on what
            # the genconf arm could possibly achieve.
            rms = rdMolAlign.GetBestRMS(Chem.Mol(Chem.RemoveHs(molh)), Chem.Mol(ref))
            row["conformer_rmsd_to_crystal"] = round(float(rms), 3)
        except Exception:
            pass

    ok, err = _write_pdbqt(molh, pdbqt)
    if not ok:
        row.update({"status": "failed", "reason": f"Meeko: {err}",
                    "elapsed_s": round(time.time() - t0, 2)})
        return row
    row.update({"pdbqt_bytes": pdbqt.stat().st_size,
                "elapsed_s": round(time.time() - t0, 2)})
    return row


# ---------------------------------------------------------------------------
# cross-docking: fetch a partner structure and carry the ligand into its frame
# ---------------------------------------------------------------------------
def fetch_cif(pdb_id: str, dest: Path, retries: int = 5) -> bool:
    if dest.is_file() and dest.stat().st_size > 0:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = RCSB_CIF.format(pdb=pdb_id.upper())
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "vina_gnina_pose_benchmark_pipeline"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                raw = fh.read()
            with gzip.open(io.BytesIO(raw)) as gz:
                dest.write_bytes(gz.read())
            return True
        except (urllib.error.URLError, OSError, EOFError):
            time.sleep(1.5 * (attempt + 1))
    return False


def cif_to_receptor_pdb(cif: Path, out_pdb: Path, drop_ccds: set[str]) -> dict:
    """Write a PDB receptor from an mmCIF, dropping named components.

    Returns the centroid of the dropped ligand as well, because that is where
    05_define_boxes.py puts the cross-docking box. The partner's own bound
    ligand is the only thing an experimenter would have to go on: they have a
    holo structure of the target with some other compound in it, and they know
    where that compound sat. Using the transformed reference ligand instead
    would hand the arm the answer.

    Uses the first model only. A cross-docking receptor taken from an NMR
    ensemble would otherwise contribute twenty overlapping copies of itself,
    and the resolution filter in 03 already excludes NMR entries, so this is
    belt and braces.
    """
    import numpy as np
    from Bio.PDB import MMCIFParser, PDBIO, Select

    class Keep(Select):
        def accept_model(self, model):
            return model.id == 0

        def accept_residue(self, residue):
            name = residue.get_resname().strip().upper()
            if name in drop_ccds:
                return False
            if name in ("HOH", "DOD"):
                return False
            return name not in L.CRYSTALLISATION_ADDITIVES or name in L.STANDARD_RESIDUES

        def accept_atom(self, atom):
            # Alternate locations: keep the first, which is the highest
            # occupancy by PDB convention.
            return (not atom.is_disordered()) or atom.get_altloc() in (" ", "A")

        def accept_chain(self, chain):
            # Biopython's PDBIO cannot write a chain identifier longer than one
            # character, and mmCIF entries for large assemblies use identifiers
            # like 'AAA'. Six of the 253 cross-docking partners are lost this
            # way. Remapping the identifiers would keep them and is not done
            # here; the failures are recorded in crossdock_align.tsv with the
            # exception text rather than dropped quietly.
            return True

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("x", str(cif))

    # Centroid of the largest copy of the dropped ligand, taken before the file
    # is written. Several entries hold the same ligand in two or four chains;
    # the copies sit in different pockets, so averaging over all of them would
    # put the box between them, in the protein core. The copy with the most
    # atoms is taken, and where copies tie the first in file order wins.
    best_coords = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip().upper() in drop_ccds:
                    coords = [a.get_coord() for a in res
                              if a.element.strip().upper() not in ("H", "D")]
                    if len(coords) > len(best_coords):
                        best_coords = coords
        break   # first model only

    io_ = PDBIO()
    io_.set_structure(structure)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io_.save(str(out_pdb), select=Keep())
    n = sum(1 for line in out_pdb.read_text().splitlines()
            if line[:6] in ("ATOM  ", "HETATM"))
    info = {"atoms_written": n, "partner_ligand_atoms": len(best_coords)}
    if best_coords:
        c = np.asarray(best_coords, dtype=float).mean(axis=0)
        info.update({"partner_ligand_center_x": round(float(c[0]), 3),
                     "partner_ligand_center_y": round(float(c[1]), 3),
                     "partner_ligand_center_z": round(float(c[2]), 3)})
    return info


def superpose_query_onto_receptor(query_pdb: Path, receptor_pdb: Path):
    """Return (rotation, translation, stats) mapping query frame to receptor frame.

    Chains are matched by sequence rather than by chain id or residue number.
    Two entries of the same protein routinely disagree about both: a construct
    numbered from the full-length sequence against one numbered from the
    expressed fragment gives a superposition on the wrong residues, silently,
    and the cross-docking reference pose would then be wrong by however far the
    offset is. Aligning the sequences costs a few milliseconds and removes that
    failure mode.
    """
    import numpy as np
    from Bio import Align
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa
    from Bio.Data.IUPACData import protein_letters_3to1

    def chains_with_ca(path: Path):
        s = PDBParser(QUIET=True).get_structure("x", str(path))
        model = next(iter(s))
        out = {}
        for ch in model:
            seq, cas = [], []
            for res in ch:
                if not is_aa(res, standard=False):
                    continue
                if "CA" not in res:
                    continue
                three = res.get_resname().strip().capitalize()
                one = protein_letters_3to1.get(three, "X")
                seq.append(one)
                cas.append(res["CA"].get_coord())
            if len(cas) >= 20:
                out[ch.id] = ("".join(seq), np.asarray(cas, dtype=float))
        return out

    qc = chains_with_ca(query_pdb)
    rc = chains_with_ca(receptor_pdb)
    if not qc or not rc:
        return None, None, {"reason": "no chain with at least 20 C-alpha atoms"}

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

    used_receptor: set[str] = set()
    qpts, rpts, matched_chains, n_identical, n_aligned = [], [], [], 0, 0

    for qid, (qseq, qca) in sorted(qc.items()):
        best = None
        for rid, (rseq, rca) in sorted(rc.items()):
            if rid in used_receptor:
                continue
            try:
                aln = aligner.align(qseq, rseq)[0]
            except Exception:
                continue
            qi, ri = aln.aligned[0], aln.aligned[1]
            pairs = sum(b - a for a, b in qi)
            if best is None or pairs > best[0]:
                best = (pairs, rid, aln)
        if best is None or best[0] < 20:
            continue
        pairs, rid, aln = best
        used_receptor.add(rid)
        matched_chains.append(f"{qid}:{rid}")
        rseq, rca = rc[rid]
        for (qa, qb), (ra, rb) in zip(aln.aligned[0], aln.aligned[1]):
            for k in range(qb - qa):
                i, j = qa + k, ra + k
                if i >= len(qca) or j >= len(rca):
                    continue
                qpts.append(qca[i])
                rpts.append(rca[j])
                n_aligned += 1
                if qseq[i] == rseq[j]:
                    n_identical += 1

    if len(qpts) < 20:
        return None, None, {"reason": f"only {len(qpts)} C-alpha pairs matched"}

    P = np.asarray(qpts)
    Q = np.asarray(rpts)
    # Kabsch. The ligand has to move with the query protein, so the transform is
    # the one that carries query coordinates onto receptor coordinates.
    pc, qc_ = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc_)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = qc_ - R @ pc
    rms = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(axis=1).mean()))
    stats = {
        "chains_matched": ";".join(matched_chains),
        "ca_pairs": len(qpts),
        "ca_rmsd_after_superposition": round(rms, 3),
        "sequence_identity_of_matched": round(n_identical / max(n_aligned, 1), 4),
        "reason": "",
    }
    return R, t, stats


def transform_sdf(src_sdf: Path, dst_sdf: Path, R, t) -> bool:
    import numpy as np
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    mol = L.load_heavy(src_sdf)
    if mol is None:
        return False
    conf = mol.GetConformer()
    pos = conf.GetPositions()
    new = (R @ pos.T).T + t
    for i in range(mol.GetNumAtoms()):
        x, y, z = new[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    dst_sdf.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(dst_sdf))
    w.write(mol)
    w.close()
    return True


# ---------------------------------------------------------------------------
# workers
# ---------------------------------------------------------------------------
_CTX: dict = {}


def _init(ctx):
    global _CTX
    _CTX = ctx


def _receptor_job(item):
    scope, key, src, drop_ccd = item
    out_dir = Path(_CTX["prepared"]) / "receptors" / scope / key
    row = prepare_receptor(Path(src), out_dir, drop_ccd, _CTX["tools"], _CTX["force"])
    row.update({"scope": scope, "key": key, "source_pdb": src})
    return row


def _ligand_job(item):
    dataset, complex_id, sdf, conformer = item
    out_dir = Path(_CTX["prepared"]) / "ligands" / conformer / dataset
    row = prepare_ligand(complex_id, Path(sdf), out_dir, conformer,
                         _CTX["seed"], _CTX["force"])
    row["dataset"] = dataset
    return row


def _crossdock_job(item):
    pair = item
    data = Path(_CTX["data"])
    prepared = Path(_CTX["prepared"])
    pid = pair["pair_id"]
    rid = pair["receptor_pdb_id"].upper()
    row = {"pair_id": pid, "query_pdb_id": pair["query_pdb_id"],
           "receptor_pdb_id": rid, "status": "ok", "reason": "",
           "recorded": datetime.datetime.now().astimezone().isoformat(timespec="seconds")}

    cif = data / "crossdock" / "cif" / f"{rid}.cif"
    if not fetch_cif(rid, cif):
        row.update({"status": "failed", "reason": "RCSB mmCIF download failed"})
        return row

    # The partner's own ligand comes out: we are docking a different molecule
    # into that pocket. Its position is still used, in 05, to place the box.
    recpdb = data / "crossdock" / "receptor_pdb" / f"{rid}.pdb"
    try:
        info = cif_to_receptor_pdb(cif, recpdb, {pair["receptor_ligand_ccd"].upper()})
    except Exception as exc:
        row.update({"status": "failed", "reason": f"mmCIF parse failed: {exc!r}"[:200]})
        return row
    row.update(info)
    if not info.get("partner_ligand_atoms"):
        row.update({"status": "failed",
                    "reason": f"ligand {pair['receptor_ligand_ccd']} not found in "
                              f"{rid}; the search API and the coordinate file disagree"})
        return row

    query_pdb = Path(pair["query_protein_pdb"])
    try:
        R, t, stats = superpose_query_onto_receptor(query_pdb, recpdb)
    except Exception as exc:
        row.update({"status": "failed", "reason": f"superposition raised: {exc!r}"[:200]})
        return row
    row.update(stats)
    if R is None:
        row["status"] = "failed"
        return row
    # A superposition that lands above 3 Angstrom on matched C-alpha atoms is
    # not the same fold in the same conformation, and a reference pose derived
    # from it would be wrong by that much before docking even starts. Three is
    # arbitrary; it is roughly where a domain has moved rather than a loop.
    #
    # This cut is the least comfortable decision in the pipeline and it biases
    # the cross-docking arm. It removed 39 of 253 pairs, and they are not a
    # random 39: they are the pairs whose two structures differ most, which is
    # to say the hardest and most interesting cross-docking cases. Excluding
    # them makes arm A4 easier than the population it is meant to represent, so
    # the A4 number is an optimistic estimate of cross-docking performance
    # rather than a neutral one.
    #
    # The fix is a local superposition on binding-site residues instead of a
    # global one on every matched C-alpha. That would give a reliable reference
    # pose even where a distant domain has swung, and would keep most of those
    # 39. It is the obvious next improvement and it is not done here. The count
    # and the per-pair C-alpha RMSD are in
    # results/preparation/crossdock_align.tsv so the size of the effect is
    # visible rather than implied.
    if stats["ca_rmsd_after_superposition"] > 3.0:
        row.update({"status": "failed",
                    "reason": f"C-alpha RMSD {stats['ca_rmsd_after_superposition']} A "
                              f"after superposition exceeds the 3.0 A limit"})
        return row

    ref_out = prepared / "reference" / "crossdock" / f"{pid}.sdf"
    if not transform_sdf(Path(pair["query_ligand_sdf"]), ref_out, R, t):
        row.update({"status": "failed", "reason": "reference ligand would not transform"})
        return row

    out_dir = prepared / "receptors" / "crossdock" / rid
    rrow = prepare_receptor(recpdb, out_dir, None, _CTX["tools"], _CTX["force"])
    if rrow["status"] == "failed":
        row.update({"status": "failed", "reason": "receptor prep: " + rrow["reason"]})
    return row


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--dataset", default="", help="posebusters, astex, or empty for both")
    ap.add_argument("--what", default="all",
                    choices=["all", "receptors", "ligands", "crossdock"])
    ap.add_argument("--jobs", type=int, default=0, help="0 means take JOBS from project.conf")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    config_dir = Path(conf["CONFIG_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    prepared = data_dir / "prepared"
    seed = L.conf_int(conf, "SEED", 20260922)
    jobs = args.jobs or L.conf_int(conf, "JOBS", 1)

    tools = {
        "pdb2pqr": L.which_or_die("pdb2pqr", "conda env vgb_bench"),
        "obabel": L.which_or_die("obabel", "conda env vgb_bench"),
        # mk_prepare_receptor is no longer called for receptors (see the note in
        # prepare_receptor) but its presence is still checked, because its
        # absence means the environment is not the one this pipeline was built
        # against and the ligand side uses the same Meeko install.
        "mk_receptor": L.which_or_die("mk_prepare_receptor.py", "conda env vgb_bench"),
    }
    ctx = {"prepared": str(prepared), "data": str(data_dir), "tools": tools,
           "seed": seed, "force": args.force}

    datasets = [args.dataset] if args.dataset else ["posebusters", "astex"]
    manifests = {}
    for ds in datasets:
        p = config_dir / f"dataset_{ds}.tsv"
        if not p.is_file():
            L.eprint(f"[error] {p} missing; run 02_fetch_benchmarks.sh first")
            return 1
        rows = L.read_tsv(p)
        manifests[ds] = rows[:args.limit] if args.limit else rows

    out = results_dir / "preparation"
    out.mkdir(parents=True, exist_ok=True)
    t_stage = time.time()

    def run_pool(fn, items, label):
        if not items:
            return []
        if jobs > 1:
            with mp.Pool(jobs, initializer=_init, initargs=(ctx,)) as pool:
                got = []
                for i, r in enumerate(pool.imap_unordered(fn, items, chunksize=1), 1):
                    got.append(r)
                    if i % 50 == 0 or i == len(items):
                        print(f"  {label}: {i}/{len(items)}", flush=True)
                return got
        _init(ctx)
        got = []
        for i, it in enumerate(items, 1):
            got.append(fn(it))
            if i % 50 == 0 or i == len(items):
                print(f"  {label}: {i}/{len(items)}", flush=True)
        return got

    # --- receptors ---------------------------------------------------------
    if args.what in ("all", "receptors"):
        items = []
        for ds, rows in manifests.items():
            for r in rows:
                items.append((ds, r["complex_id"], r["protein_pdb"], r["ccd_id"]))
        print(f"[04_prepare] {len(items)} receptors, {jobs} job(s)")
        rrows = run_pool(_receptor_job, items, "receptors")
        L.write_tsv(out / "receptor_prep.tsv", RECEPTOR_COLUMNS, rrows)
        bad = [r for r in rrows if r["status"] == "failed"]
        print(f"[04_prepare] receptors: {len(rrows) - len(bad)} ok, {len(bad)} failed")
        fb = [r for r in rrows if r.get("protonation", "").startswith("obabel")]
        if fb:
            print(f"[04_prepare] {len(fb)} receptors fell back to the obabel protonator")

    # --- ligands -----------------------------------------------------------
    if args.what in ("all", "ligands"):
        items = []
        for ds, rows in manifests.items():
            for r in rows:
                for conformer in ("refconf", "genconf"):
                    items.append((ds, r["complex_id"], r["ligand_sdf"], conformer))
        print(f"[04_prepare] {len(items)} ligand preparations, {jobs} job(s)")
        lrows = run_pool(_ligand_job, items, "ligands")
        L.write_tsv(out / "ligand_prep.tsv", LIGAND_COLUMNS, lrows)
        bad = [r for r in lrows if r["status"] == "failed"]
        print(f"[04_prepare] ligands: {len(lrows) - len(bad)} ok, {len(bad)} failed")
        gen = [r for r in lrows if r["conformer"] == "genconf"
               and isinstance(r.get("conformer_rmsd_to_crystal"), float)]
        if gen:
            vals = sorted(r["conformer_rmsd_to_crystal"] for r in gen)
            mid = vals[len(vals) // 2]
            print(f"[04_prepare] generated conformer vs crystal conformer, "
                  f"median best-fit RMSD {mid:.2f} A over {len(vals)} ligands")

        # The reference poses every arm scores against: the crystal ligand,
        # heavy atoms, copied once so later stages never touch the archive.
        for ds, rows in manifests.items():
            dest = prepared / "reference" / ds
            dest.mkdir(parents=True, exist_ok=True)
            for r in rows:
                d = dest / f"{r['complex_id']}.sdf"
                if not d.is_file() or args.force:
                    shutil.copyfile(r["ligand_sdf"], d)

    # --- cross-docking -----------------------------------------------------
    if args.what in ("all", "crossdock"):
        pairs_path = config_dir / "crossdock_pairs.tsv"
        if not pairs_path.is_file():
            print("[04_prepare] no crossdock_pairs.tsv; run 03_build_crossdock_set.py "
                  "first. Skipping the cross-docking arm.")
        else:
            pairs = L.read_tsv(pairs_path)
            if args.limit:
                pairs = pairs[:args.limit]
            by_id = {r["complex_id"]: r for rows in manifests.values() for r in rows}
            items, missing = [], 0
            for p in pairs:
                q = by_id.get(p["query_complex_id"])
                if q is None:
                    missing += 1
                    continue
                p = dict(p)
                p["query_protein_pdb"] = q["protein_pdb"]
                p["query_ligand_sdf"] = q["ligand_sdf"]
                items.append(p)
            print(f"[04_prepare] {len(items)} cross-docking pairs, {jobs} job(s)"
                  + (f" ({missing} had no manifest row)" if missing else ""))
            crows = run_pool(_crossdock_job, items, "crossdock")
            L.write_tsv(out / "crossdock_align.tsv", ALIGN_COLUMNS, crows)
            bad = [r for r in crows if r["status"] == "failed"]
            print(f"[04_prepare] cross-docking: {len(crows) - len(bad)} prepared, "
                  f"{len(bad)} failed")
            if bad:
                from collections import Counter
                for reason, n in Counter(r["reason"].split(":")[0] for r in bad).most_common(6):
                    print(f"    {n:4d}  {reason}")

    print(f"[04_prepare] stage took {time.time() - t_stage:.0f} s; "
          f"prepared tree is {L.du_kb(prepared) / 1024:.0f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
