# tdp43_neuron_fans
Long-read (ONT) and short-read (10x Genomics 3' OCM) single-nucleus RNA-seq of
neuronal nuclei sorted from fresh-frozen human brain into TDP-43+ and TDP-43-
populations.

## ONT pre-processing

Project directory on Quest: `/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans`

| File | What it holds |
| --- | --- |
| `config/samples.tsv` | one row per sequencing run (folder under `pod5/`) and the 10x pool it belongs to |
| `config/pools.tsv` | one row per 10x pool, with `expected_cells` for wf-single-cell |

1. `jobs/01_basecalling.sh`: SLURM array, one task per row of `samples.tsv`.
   Dorado basecalls each run separately to `basecalled/<run>.fastq.gz`.
2. `jobs/02_epi2me_single_cell.sh`: SLURM array, one task per row of
   `pools.tsv`. Combines the FASTQs of every run in a pool and runs
   epi2me-labs/wf-single-cell once per pool, writing to `epi2me_outs/<pool>/`.

```bash
mkdir -p /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/logs
cd /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans
jid=$(sbatch --parsable jobs/01_basecalling.sh)
sbatch --dependency=afterok:${jid} jobs/02_epi2me_single_cell.sh
```

If you add or remove runs or pools, update the two TSVs and the
`#SBATCH --array` range of each script to match.
