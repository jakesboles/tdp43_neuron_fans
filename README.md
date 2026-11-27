# tdp43_neuron_fans
Long-read (ONT) and short-read (10x Genomics 3' OCM) single-nucleus RNA-seq of
neuronal nuclei sorted from fresh-frozen human brain into TDP-43+ and TDP-43-
populations.

## ONT pre-processing

Storage on Quest:

| Location | Holds |
| --- | --- |
| `/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans` (large, transient) | `pod5/`, `basecalled/`, `cellranger/`, and all epi2me data (`epi2me_outs/`, `epi2me_runs/`, `epi2me_work/`, `epi2me_store/`), `nf_work/` |
| `/projects/b1169/boles/tdp43_neuron_fans` (kept) | `logs/` (SLURM and Nextflow logs), `nf_outs/`, processed objects you copy over |

Sample sheets in this repo (submit every job **from the repo root**):

| File | What it holds |
| --- | --- |
| `config/samples.tsv` | one row per sequencing run (folder under `pod5/`) and the 10x pool it belongs to |
| `config/pools.tsv` | one row per 10x pool, with `expected_cells` |

1. `jobs/01_basecalling.sh`: one GPU task per run (row of `samples.tsv`).
   It writes unaligned BAM to `basecalled/<run>.bam` and skips runs already done.
   If a task times out, resubmit the same array index: it resumes from the
   partial BAM (`dorado --resume-from`) instead of starting over.
2. `jobs/02_epi2me_single_cell.sh`: one task per pool (row of `pools.tsv`).
   It runs epi2me-labs/wf-single-cell (`--bam`) on the BAMs of every run in the pool.

```bash
mkdir -p /projects/b1169/boles/tdp43_neuron_fans/logs
jid=$(sbatch --parsable jobs/01_basecalling.sh)          # or --array=5-12 for specific rows
sbatch --dependency=afterok:${jid} jobs/02_epi2me_single_cell.sh
```

`jobs/nf.sh <pool>` runs the separate 10x_LR_preprocessing Nextflow pipeline,
which does its own basecalling, on one pool. It is an alternative to steps
1-2 and is not needed for them.

## Short-read scRNA-seq (10x 3' OCM)

R scripts in `scripts/sc_scripts/`. Every script runs from the project root,
`/projects/b1169/boles/tdp43_neuron_fans`, and reads Cell Ranger output for
every pool with a `cellranger_dir` in `config/pools.tsv`. Pools still set to
`NA` are skipped.

| Output | Location |
| --- | --- |
| Data objects (BPCells matrices, metadata and reduction RDS files) | `data/sc_data/<script>/` |
| Tables and PNG plots | `results/sc_results/<script>/[<cell type>/]` |
| External reference tables (GEO, Novartis) | `external_data/` |

Only `01_qc.R` writes a counts matrix (`data/sc_data/01_qc/bpcells`, real
counts for QC-passing cells). Later scripts rebuild their object from it,
subset to their own cells, and re-normalize, plus their own `metadata.rds`
and reductions.

1. `01_qc.R`: build, QC, and filter all pools
2. `02_doubletfinder.R`: SLURM array, one task per sample
   (`jobs/sc_jobs/02_doubletfinder.sh`; set `--array` to the number of samples)
3. `03_integration.R`: Harmony integration, clustering, markers
4. `04_annotation.R`: interactive labeling; cluster-to-label maps must be redone
   whenever clustering changes
5. `05_tdp43_de_{pseudobulk,mast,wilcox}.R`: TDPneg vs TDPpos within a cell type
6. `06_liu_de.R`, `07_meta_de.R`, `summary_plots.R`,
   `tdp43neg_vs_novartis_data.R`: comparisons and figures
