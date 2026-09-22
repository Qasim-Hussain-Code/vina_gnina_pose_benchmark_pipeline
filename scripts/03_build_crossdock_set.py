#!/usr/bin/env python3
"""Build the cross-docking set by querying the RCSB, not by downloading a list.

Cross-docking is what virtual screening actually does. You have one structure of
a target, solved with whatever ligand happened to crystallise, and you dock a
library of other compounds into it. The receptor conformation you have is the one
that some other molecule induced. Self-docking never tests that, because the
pocket you dock into was shaped by the very ligand you are putting back.

So: for each complex in the primary set, find another PDB entry of the same
protein with a different ligand bound, and dock the primary ligand into that
other receptor. The RMSD reference stays the primary complex's own crystal pose,
after superposing the two receptors, because that is where the ligand actually
sits when this protein binds it.

How a pair is built
-------------------
1. Ask the RCSB Data API for the query entry's polymer entities and read the
   UniProt accession out of the SIFTS cross-reference, plus the entity's
   membership of the 95 per cent sequence-identity cluster that the RCSB
   computes. Using their clustering rather than running an alignment here means
   the identity cutoff is theirs, is applied the same way to every entry, and
   does not depend on my choice of alignment parameters.
2. Ask the Search API for every X-ray entry carrying that UniProt accession at
   the resolution cutoff in project.conf.
3. Keep candidates that share the query's 95 per cent cluster, so "same protein"
   is a sequence statement and not a name match.
4. Drop candidates whose bound ligands are only crystallisation additives. An
   apo structure with three glycerols in it is not a holo structure of anything,
   and the box for the cross-docking arm is placed on the partner's own ligand.
5. Drop candidates bound to the same chemical component as the query, because
   that is self-docking wearing a different PDB code.
6. Of what survives, take the highest-resolution partner. One partner per query,
   set by XDOCK_MAX_PARTNERS, so that no single well-studied target contributes
   fifty pairs and dominates the arm.

Every complex that fails any of these steps is written to the exclusion table
with which step rejected it. The counts in the README come from that table.

Usage
-----
    python scripts/03_build_crossdock_set.py --config project.conf
    python scripts/03_build_crossdock_set.py --config project.conf --force
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

DATA_GRAPHQL = "https://data.rcsb.org/graphql"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_CIF = "https://files.rcsb.org/download/{pdb_id}.cif.gz"

# Components that are almost always crystallisation additives, cryoprotectants,
# buffers or ions rather than the thing the structure is about. A candidate whose
# only heteroatoms come from this list is treated as having no ligand of interest.
#
# This list is a judgement call and it is the weakest link in the cross-docking
# set. It is deliberately short and conservative: including something here that
# is genuinely a ligand loses a pair, which is visible in the counts, while
# leaving a true additive out gives a pair whose box is centred on a glycerol,
# which is not visible in the counts at all. I would rather lose pairs.
CRYSTALLISATION_ADDITIVES = {
    # polyols and cryoprotectants
    "GOL", "EDO", "PEG", "PGE", "PG4", "PG0", "P6G", "1PE", "2PE", "MPD",
    "TRS", "BME", "DTT", "DTU", "MES", "EPE", "IMD", "BTB", "TAR", "MLI",
    "DMS", "DMF", "ACN", "MOH", "EOH", "IPA", "ACT", "ACY", "FMT", "OXL",
    "CIT", "FLC", "SIN", "SUC", "MAL", "GLC", "BGC", "FUC", "XYP",
    # ions and simple inorganics
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2",
    "CU", "CU1", "NI", "CO", "CD", "HG", "BR", "IOD", "F", "NO3", "CO3",
    "AZI", "CN", "SCN", "WO4", "MOO", "VO4", "BEF", "ALF", "SF4", "FES",
    # unassigned density, solvent, free amino acids and chelators
    "UNL", "UNX", "UNK", "HOH", "DOD", "NH4", "GLY", "ALA", "SER", "EDT",
}

# A component with fewer heavy atoms than this is not a drug-like ligand and is
# not worth centring a 25 Angstrom box on. Ten is arbitrary. It sits just below
# the smallest ligands in the PoseBusters set and just above most buffers.
MIN_LIGAND_HEAVY_ATOMS = 10


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------
class Rcsb:
    """Thin client with retries and a polite delay.

    The RCSB asks for reasonable request rates and does not publish a hard
    limit. 0.12 s between calls with 8 retries on backoff has run the whole 308
    without a single 429 here; a tighter loop started collecting them.
    """

    def __init__(self, delay: float = 0.12, retries: int = 8, cache: Path | None = None):
        self.delay = delay
        self.retries = retries
        self.cache = cache
        if cache:
            cache.mkdir(parents=True, exist_ok=True)
        self.calls = 0

    def _post(self, url: str, payload: dict) -> dict:
        body = json.dumps(payload).encode()
        last = None
        for attempt in range(self.retries):
            try:
                req = urllib.request.Request(
                    url, data=body,
                    headers={"Content-Type": "application/json",
                             "User-Agent": "vina_gnina_pose_benchmark_pipeline"},
                )
                with urllib.request.urlopen(req, timeout=90) as fh:
                    self.calls += 1
                    time.sleep(self.delay)
                    return json.load(fh)
            except urllib.error.HTTPError as exc:
                # 204 from the Search API means no hits, which is an answer.
                if exc.code == 204:
                    return {}
                last = exc
                if exc.code in (429, 500, 502, 503, 504):
                    time.sleep(1.5 * (attempt + 1))
                    continue
                raise
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
                last = exc
                time.sleep(1.5 * (attempt + 1))
        raise RuntimeError(f"RCSB request failed after {self.retries} tries: {last}")

    def graphql(self, query: str, variables: dict, cache_key: str | None = None) -> dict:
        if self.cache and cache_key:
            p = self.cache / f"{cache_key}.json"
            if p.is_file():
                return json.loads(p.read_text())
        out = self._post(DATA_GRAPHQL, {"query": query, "variables": variables})
        if self.cache and cache_key:
            (self.cache / f"{cache_key}.json").write_text(json.dumps(out))
        return out

    def search(self, payload: dict, cache_key: str | None = None) -> dict:
        if self.cache and cache_key:
            p = self.cache / f"{cache_key}.json"
            if p.is_file():
                return json.loads(p.read_text())
        out = self._post(SEARCH_URL, payload)
        if self.cache and cache_key:
            (self.cache / f"{cache_key}.json").write_text(json.dumps(out))
        return out


ENTRY_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_accession_info { initial_release_date deposit_date }
    rcsb_entry_info { resolution_combined experimental_method }
    polymer_entities {
      rcsb_id
      entity_poly { rcsb_sample_sequence_length }
      rcsb_polymer_entity_container_identifiers {
        auth_asym_ids
        reference_sequence_identifiers { database_accession database_name }
      }
      rcsb_polymer_entity_group_membership { group_id similarity_cutoff aggregation_method }
    }
    nonpolymer_entities {
      rcsb_id
      nonpolymer_comp {
        chem_comp { id formula_weight name }
        rcsb_chem_comp_descriptor { SMILES }
      }
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}
"""


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def heavy_atom_count(formula_weight, smiles: str | None) -> int:
    """Rough heavy-atom count from SMILES, falling back to weight.

    RDKit is available in this environment but importing it here just to count
    atoms in a few thousand short strings costs more than it saves, and the
    number only feeds a coarse size filter. Non-hydrogen element symbols are
    counted directly; if there is no SMILES, molecular weight over 16 is used as
    a stand-in, 16 being roughly the mass per heavy atom of an organic molecule.
    """
    if smiles:
        n = 0
        i = 0
        while i < len(smiles):
            c = smiles[i]
            if c.isalpha():
                if c in "cnops":            # aromatic lower case, all heavy
                    n += 1
                elif c == "H":
                    pass
                elif c.isupper():
                    # Two-letter element symbols such as Cl, Br, Si, Se.
                    if i + 1 < len(smiles) and smiles[i + 1].islower() and \
                            smiles[i:i + 2] in ("Cl", "Br", "Si", "Se", "Na", "Mg",
                                                "Ca", "Fe", "Zn", "Mn", "Cu", "Co",
                                                "Ni", "Al", "As", "Sn", "Li"):
                        n += 1
                        i += 1
                    else:
                        n += 1
            i += 1
        return n
    try:
        return int(float(formula_weight) / 16.0)
    except (TypeError, ValueError):
        return 0


def interesting_ligands(entry: dict) -> list[dict]:
    """Non-polymer components of an entry that could be a ligand of interest."""
    out = []
    for npe in entry.get("nonpolymer_entities") or []:
        comp = (npe.get("nonpolymer_comp") or {})
        cc = comp.get("chem_comp") or {}
        cid = cc.get("id")
        if not cid or cid in CRYSTALLISATION_ADDITIVES:
            continue
        smiles = (comp.get("rcsb_chem_comp_descriptor") or {}).get("SMILES")
        n_heavy = heavy_atom_count(cc.get("formula_weight"), smiles)
        if n_heavy < MIN_LIGAND_HEAVY_ATOMS:
            continue
        chains = (npe.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []
        out.append({
            "ccd_id": cid,
            "name": cc.get("name") or "",
            "formula_weight": cc.get("formula_weight"),
            "heavy_atoms": n_heavy,
            "auth_asym_ids": chains,
        })
    return out


def uniprot_and_cluster(entry: dict, cutoff: float = 95.0):
    """Return {accession: cluster_id} for each polymer entity carrying one."""
    out = {}
    for pe in entry.get("polymer_entities") or []:
        ids = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        refs = ids.get("reference_sequence_identifiers") or []
        accs = [r["database_accession"] for r in refs
                if r.get("database_name") == "UniProt" and r.get("database_accession")]
        if not accs:
            continue
        cluster = None
        for g in pe.get("rcsb_polymer_entity_group_membership") or []:
            if g.get("aggregation_method") == "sequence_identity" and \
                    abs(float(g.get("similarity_cutoff", -1)) - cutoff) < 1e-6:
                cluster = g.get("group_id")
                break
        for a in accs:
            out[a] = cluster
    return out


def resolution_of(entry: dict):
    r = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
    if isinstance(r, list) and r:
        return float(r[0])
    if isinstance(r, (int, float)):
        return float(r)
    return None


def search_by_uniprot(client: Rcsb, accession: str, max_res: float) -> list[str]:
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers"
                                 ".reference_sequence_identifiers.database_accession",
                    "operator": "exact_match", "value": accession}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal", "value": max_res}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.experimental_method",
                    "operator": "exact_match", "value": "X-ray"}},
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 200},
                            "results_content_type": ["experimental"]},
    }
    out = client.search(payload, cache_key=f"search_{accession}_{max_res}")
    return [r["identifier"] for r in out.get("result_set", [])]


def load_conf(path: Path) -> dict:
    """Read the shell project.conf into a dict.

    project.conf is sourced by the bash stages, so it is shell syntax rather
    than ini or yaml. Only KEY=VALUE lines are read and inline comments are
    stripped; nothing here needs shell expansion.
    """
    conf = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        if not key.replace("_", "").isalnum():
            continue
        val = val.split("#", 1)[0].strip().strip('"').strip("'")
        conf[key] = val
    return conf


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path, help="project.conf")
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if config/crossdock_pairs.tsv exists")
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N query complexes (for a smoke test)")
    args = ap.parse_args()

    conf = load_conf(args.config)
    config_dir = Path(conf["CONFIG_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    max_res = float(conf.get("XDOCK_MAX_RESOLUTION", 2.5))
    min_ident = float(conf.get("XDOCK_MIN_SEQ_IDENTITY", 0.95)) * 100.0
    max_partners = int(conf.get("XDOCK_MAX_PARTNERS", 1))

    out_pairs = config_dir / "crossdock_pairs.tsv"
    out_excl = config_dir / "crossdock_excluded.tsv"
    if out_pairs.is_file() and not args.force:
        print(f"[03_build_crossdock_set] {out_pairs.name} exists; skipping "
              f"(--force to rebuild).")
        return 0

    manifest = config_dir / "dataset_posebusters.tsv"
    if not manifest.is_file():
        print("[03_build_crossdock_set] config/dataset_posebusters.tsv missing; "
              "run 02_fetch_benchmarks.sh first.", file=sys.stderr)
        return 1
    with manifest.open() as fh:
        queries = list(csv.DictReader(fh, delimiter="\t"))
    if args.limit:
        queries = queries[:args.limit]

    cache = data_dir / "cache" / "rcsb"
    client = Rcsb(cache=cache)
    stamp = datetime.date.today().isoformat()

    # --- entry records for the query complexes -------------------------------
    print(f"[03_build_crossdock_set] fetching entry records for "
          f"{len(queries)} query complexes")
    entries: dict[str, dict] = {}
    q_pdbs = sorted({q["pdb_id"].upper() for q in queries})
    for batch in chunked(q_pdbs, 50):
        data = client.graphql(ENTRY_QUERY, {"ids": batch},
                              cache_key="entry_" + "_".join(batch[:1]) + f"_{len(batch)}")
        for e in (data.get("data", {}).get("entries") or []):
            if e:
                entries[e["rcsb_id"].upper()] = e

    excluded: list[tuple[str, str]] = []
    # accession -> candidate entry ids, so one search serves every query that
    # maps to the same protein
    candidates_by_acc: dict[str, list[str]] = {}

    stage_counts = Counter()
    pair_rows = []

    for i, q in enumerate(queries, 1):
        cid, pdb, ccd = q["complex_id"], q["pdb_id"].upper(), q["ccd_id"].upper()
        if i % 25 == 0 or i == len(queries):
            print(f"  {i}/{len(queries)} queries, {len(pair_rows)} pairs so far, "
                  f"{client.calls} API calls")

        qe = entries.get(pdb)
        if qe is None:
            excluded.append((cid, "query entry not returned by the RCSB Data API; "
                                  "probably obsoleted since the benchmark was built"))
            stage_counts["no_query_entry"] += 1
            continue

        q_map = uniprot_and_cluster(qe, min_ident)
        if not q_map:
            excluded.append((cid, "no UniProt cross-reference on any polymer entity, "
                                  "so same-protein cannot be established"))
            stage_counts["no_uniprot"] += 1
            continue

        # Gather candidates for every accession this entry carries. Multi-chain
        # complexes with two different proteins give two accessions; a partner
        # sharing either is still a structure of a protein the ligand binds.
        cand_ids: set[str] = set()
        for acc in q_map:
            if acc not in candidates_by_acc:
                candidates_by_acc[acc] = search_by_uniprot(client, acc, max_res)
            cand_ids.update(candidates_by_acc[acc])
        cand_ids.discard(pdb)
        if not cand_ids:
            excluded.append((cid, f"no other X-ray entry at {max_res} A or better "
                                  f"shares a UniProt accession with this one"))
            stage_counts["no_same_protein_entry"] += 1
            continue

        # Fetch the candidates we have not seen yet.
        todo = sorted(c.upper() for c in cand_ids if c.upper() not in entries)
        for batch in chunked(todo, 50):
            data = client.graphql(ENTRY_QUERY, {"ids": batch},
                                  cache_key="entry_" + batch[0] + f"_{len(batch)}")
            for e in (data.get("data", {}).get("entries") or []):
                if e:
                    entries[e["rcsb_id"].upper()] = e

        q_clusters = {c for c in q_map.values() if c}
        viable = []
        reasons = Counter()
        for c in sorted(cand_ids):
            ce = entries.get(c.upper())
            if ce is None:
                reasons["candidate_not_returned"] += 1
                continue
            res = resolution_of(ce)
            if res is None or res > max_res:
                reasons["resolution"] += 1
                continue
            # Sequence identity, via the RCSB's own 95 per cent clustering.
            c_map = uniprot_and_cluster(ce, min_ident)
            c_clusters = {x for x in c_map.values() if x}
            if q_clusters and not (q_clusters & c_clusters):
                reasons["below_identity_cutoff"] += 1
                continue
            ligs = interesting_ligands(ce)
            if not ligs:
                reasons["no_ligand_of_interest"] += 1
                continue
            other = [L for L in ligs if L["ccd_id"].upper() != ccd]
            if not other:
                reasons["same_ligand_only"] += 1
                continue
            # Largest ligand in the partner, as the one most likely to be the
            # thing the structure is about and the one whose centroid is the
            # most defensible box centre.
            other.sort(key=lambda L: (-L["heavy_atoms"], L["ccd_id"]))
            viable.append((res, c.upper(), other[0]))

        if not viable:
            worst = reasons.most_common(1)[0][0] if reasons else "no candidates survived"
            detail = ", ".join(f"{k}={v}" for k, v in sorted(reasons.items()))
            excluded.append((cid, f"no partner survived filtering ({detail})"))
            stage_counts[f"no_partner_{worst}"] += 1
            continue

        viable.sort(key=lambda t: (t[0], t[1]))      # best resolution first
        for res, partner, lig in viable[:max_partners]:
            pe = entries[partner]
            pair_rows.append({
                "pair_id": f"{cid}__into__{partner}",
                "query_complex_id": cid,
                "query_pdb_id": pdb,
                "query_ccd_id": ccd,
                "query_resolution": f"{resolution_of(qe):.2f}" if resolution_of(qe) else "NA",
                "receptor_pdb_id": partner,
                "receptor_resolution": f"{res:.2f}",
                "receptor_ligand_ccd": lig["ccd_id"],
                "receptor_ligand_heavy_atoms": lig["heavy_atoms"],
                "receptor_ligand_chains": ",".join(lig["auth_asym_ids"]),
                "uniprot": ",".join(sorted(q_map)),
                "identity_cutoff_percent": f"{min_ident:.0f}",
                "receptor_release_date": (pe.get("rcsb_accession_info") or {}).get("initial_release_date", "NA")[:10],
            })
        stage_counts["paired"] += 1

    # --- write ---------------------------------------------------------------
    cols = ["pair_id", "query_complex_id", "query_pdb_id", "query_ccd_id",
            "query_resolution", "receptor_pdb_id", "receptor_resolution",
            "receptor_ligand_ccd", "receptor_ligand_heavy_atoms",
            "receptor_ligand_chains", "uniprot", "identity_cutoff_percent",
            "receptor_release_date"]
    tmp = out_pairs.with_suffix(".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in pair_rows:
            w.writerow(r)
    shutil.move(str(tmp), str(out_pairs))

    with out_excl.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["dataset", "complex_id", "stage", "reason", "recorded"])
        for cid, why in excluded:
            w.writerow(["crossdock", cid, "03_build_crossdock_set", why, stamp])

    summary = results_dir / "crossdock_set_summary.tsv"
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["item", "value"])
        w.writerow(["queries_considered", len(queries)])
        w.writerow(["pairs_built", len(pair_rows)])
        w.writerow(["queries_paired", stage_counts["paired"]])
        w.writerow(["queries_excluded", len(excluded)])
        w.writerow(["max_resolution_angstrom", max_res])
        w.writerow(["min_sequence_identity_percent", f"{min_ident:.0f}"])
        w.writerow(["max_partners_per_query", max_partners])
        w.writerow(["min_ligand_heavy_atoms", MIN_LIGAND_HEAVY_ATOMS])
        w.writerow(["additives_excluded_count", len(CRYSTALLISATION_ADDITIVES)])
        w.writerow(["rcsb_api_calls", client.calls])
        w.writerow(["built", datetime.datetime.now().astimezone().isoformat(timespec="seconds")])
        for k, v in sorted(stage_counts.items()):
            if k != "paired":
                w.writerow([f"excluded_{k}", v])

    # Distinct receptors, because that is the number that sets the download and
    # preparation cost of the arm.
    n_receptors = len({r["receptor_pdb_id"] for r in pair_rows})
    print(f"[03_build_crossdock_set] {len(pair_rows)} pairs from {len(queries)} "
          f"queries, {n_receptors} distinct receptors to prepare")
    print(f"[03_build_crossdock_set] {len(excluded)} queries excluded; reasons in "
          f"{out_excl}")
    for k, v in sorted(stage_counts.items()):
        print(f"    {k:34s} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
