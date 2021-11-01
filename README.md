# midnax
[![DOI](https://zenodo.org/badge/404365756.svg)](https://zenodo.org/badge/latestdoi/404365756)
MItochondrial DNA eXtracter. With the help of blast extracts contig with mtDNA inside. Needs closest reference mtDNA to run.
## Install, set and run
midnax is available in conda, to install and set is use following commands:
1) Download midnax in separate conda environment: `conda create -n midnax_env -c conda-forge -c bioconda -c zilov midnax`
2) Activate the environment: `conda activate midnax_env`
3) To run midnax on your assembly use the following command:
   ```
   pannopi -a /path/to/assembly.fasta -m /path/to/mtdna.fasta -t 32 -o /path/to/outdir
   ```
## Options
```
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        path to genome asssembly in FASTA format
  -m MTDNA, --mtdna MTDNA
                        path to reference mtDNA in FASTA format
  -p PREFIX, --prefix PREFIX
                        prefix for output files
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -t THREADS, --threads THREADS
                        number of threads [default == 8]
  -d, --debug           do not run, print commands only
```
## Citation
Please cite midnax if you use it in your project!