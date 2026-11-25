# tdp43_neuron_fans
Long-read (ONT) and short-read (10x Genomics 3' OCM) single-nucleus RNA-seq of
neuronal nuclei sorted from fresh-frozen human brain into TDP-43+ and TDP-43-
populations.

## ONT pre-processing

Storage on Quest:

| Location | Holds |
| --- | --- |
| `/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans` (large, transient) | `pod5/`, `basecalled/`, Nextflow work folders |
| `/projects/b1169/boles/tdp43_neuron_fans` (kept) | `logs/`, `epi2me_outs/`, `nf_outs/`, Nextflow launch folders |

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
