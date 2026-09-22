# What belongs in data/

Nothing here is tracked by git except this file. Everything else is either
downloaded from a source recorded in `config/sources.tsv` or written by a script
in `scripts/`, so it is all reproducible from a clean clone.

On the machine this was developed on, `data/` is not inside the repository at
all. The repository sits on a filesystem with 8 GB free and the analysis needs
more than that once the GNINA binary and the benchmark archive are unpacked, so
`scripts/00_configure.sh --data-dir /home/qasim/vgb_bench_data` was used and
`DATA_DIR` in `project.conf` points there. Check `project.conf` for where it
actually is on any given machine. If you pass no `--data-dir`, it is this
directory.

## Layout

```
data/
  raw/                         the benchmark archive as downloaded, plus the
                               Zenodo record JSON and the pinned 308-member id
                               list. Verified against the md5 the Zenodo API
                               publishes, not against a checksum pasted into a
                               script.
  benchmarks/
    posebusters_benchmark_set/ 428 complex folders as the authors distribute
                               them. The 308 used are named in
                               config/dataset_posebusters.tsv; the other 120
                               are listed in config/dataset_excluded.tsv with
                               the reason.
    astex_diverse_set/         85 complex folders.
  tools/
    gnina                      the 2.1 GB static release binary. Version and
                               sha256 in results/environment/versions.tsv.
    gnina_release.json         the GitHub release record it came from.
  cache/rcsb/                  cached RCSB Search and Data API responses, so
                               re-running 03_build_crossdock_set.py does not
                               repeat 725 API calls. Safe to delete; it will be
                               refetched.
  crossdock/
    cif/                       partner structures as mmCIF, one per receptor.
    receptor_pdb/              the same, cleaned, with the partner's own ligand
                               removed. Used for docking and for the validity
                               checks on the cross-docking arm.
  prepared/
    receptors/<scope>/<key>/   receptor.pdbqt and clean.pdb per receptor.
                               clean.pdb is kept because fpocket reads it in
                               stage 5; the protonated intermediate and the PQR
                               are deleted as soon as the PDBQT exists.
    ligands/refconf/<dataset>/ ligand PDBQT built on the crystal coordinates.
    ligands/genconf/<dataset>/ ligand PDBQT built from SMILES with ETKDGv3 and
                               the seed in project.conf.
    reference/<dataset>/       the crystal ligand each pose is measured against.
    reference/crossdock/       the query ligand carried into the partner's frame
                               by the superposition in stage 4.
  poses/<arm>/<method>/<dataset>/<key>.sdf.gz
                               the top N poses per run, N from TOP_N_POSES.
                               Gzipped because this is the only thing here that
                               grows with the size of the experiment.
  work/                        per-run scratch. Every docking run writes its
                               PDBQT output here and it is deleted as soon as
                               the poses are in the archive and the score is in
                               the results table. A killed run leaves a
                               directory behind and the next run removes it.
```

## Size

Measured on this machine after a full run of five arms, three methods and both
datasets. The current figures are in `logs/*.resources.tsv`, which is what the
README quotes; the breakdown below is for planning.

| What | Size | Re-fetchable |
|---|---|---|
| GNINA binary | 2.1 GB | yes, from the GitHub release |
| benchmark archive, downloaded | 54 MB | yes, from Zenodo record 8278563 |
| benchmark archive, unpacked | 211 MB | yes |
| prepared receptors | around 600 MB for 634 receptors | yes, from stage 4 |
| cross-docking mmCIF files | a few MB per entry | yes, from the RCSB |
| pose archives | roughly 25 kB per run | yes, from stage 6 |

The disk gate in `scripts/06_dock.sh` measures the per-run cost on the first ten
complexes of an arm, projects the total for the arms still to come, and refuses
to start if it will not fit in `DISK_GB`. That is the number to trust over this
table.

## Deleting it

Safe to remove entirely between analyses:

```bash
source project.conf
rm -rf "$DATA_DIR"
```

Everything in `results/`, `figures/` and `logs/` survives, and those are what the
README and the report are built from. Re-running `bash run_all.sh` rebuilds
`data/` from scratch; the stage stamps in `logs/` will need removing first if you
want the earlier stages to run again rather than skip.
