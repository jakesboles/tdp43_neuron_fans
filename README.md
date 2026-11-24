# tdp43_neuron_fans
Long-read (ONT) and short-read (10x Genomics 3' OCM) single-nucleus RNA-seq of
neuronal nuclei sorted from fresh-frozen human brain into TDP-43+ and TDP-43-
populations.

## ONT pre-processing

Storage on Quest:

| Location | Holds |
| --- | --- |
| `/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans` (scratch) | `pod5/`, `basecalled/`, nextflow work dirs |
| `/projects/b1169/boles/tdp43_neuron_fans` (kept) | `epi2me_outs/`, `nf_outs/`, nextflow launch dirs and logs |

Sample sheets in this repo (submit every job **from the repo root**):

| File | What it holds |
| --- | --- |
| `config/samples.tsv` | one row per sequencing run (folder under `pod5/`) and the 10x pool it belongs to |
| `config/pools.tsv` | one row per 10x pool, with `expected_cells` |

1. `jobs/00_make_basecall_chunks.sh` (login node, no sbatch): splits each run's
   POD5s into chunks of N files (default 10) so every basecalling job fits in
   the 48 h GPU limit. It prints the `--array` range for step 1.
2. `jobs/01_basecalling.sh`: one GPU task per chunk. It writes plain FASTQ to
   `basecalled/<run>/<run>_partNNN.fastq`.
3. `jobs/02_epi2me_single_cell.sh`: one task per pool. It runs
   epi2me-labs/wf-single-cell on every chunk from every run in the pool.

```bash
mkdir -p /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/logs /projects/b1169/boles/tdp43_neuron_fans/logs
bash jobs/00_make_basecall_chunks.sh 10        # prints e.g. sbatch --array=1-9 ...
jid=$(sbatch --parsable --array=1-9 jobs/01_basecalling.sh)
sbatch --dependency=afterok:${jid} --array=1 jobs/02_epi2me_single_cell.sh   # pools whose runs are all basecalled
```

When more POD5s arrive, run step 0 again (not while 01 or 02 jobs are
running). Finished chunks are skipped, and it prints the range still to do.
Only chunk a run once its transfer is complete.

`jobs/nf.sh <pool>` runs the separate 10x_LR_preprocessing Nextflow pipeline,
which does its own basecalling, on one pool. It is an alternative to steps
0-2 and is not needed for them.
