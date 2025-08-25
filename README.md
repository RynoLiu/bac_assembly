# Tool for quick assemble microbial shotgun sequence

Usage:

```
bash assembly.sh \
-i $fastq_dir \
-o $output_dir
```

Parameters:

| Parameters  | Description  |
|     ---     |     ---      |
| -i, --input | Path to the input sequence directory. |
| -o, --outdir| Output directory. |
| -a, --assembly-options | Assembly options. (default: --megahit) |
| -b, --binning-options | Binning options. (default: --metabat2 --maxbin2 --concoct) |
| -r,  --refine-options | Bin refine options. (default: -c 50 -x 10) |
| -t,  --threads | Threads. (default: 1) |
| -m, --memory| memory in GB. (default: 24) |
| --clean | Clean mode. Remove all work files. (If your hard drive has less space) |
| -h,  --help| Display this help message. |
