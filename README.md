# midnax
MItochondrial DNA eXtracter. With the help of blast extracts contig with mtDNA inside. Needs closest reference mtDNA to run.
## Install, set and run
midnax is available in conda, to install and set is use following commands:
1) Download midnax in separate conda environment: `conda create -n midnax_env -c conda-forge -c bioconda -c zilov midnax`
2) Activate the environment: `conda activate midnax_env`
1) To run midnax on your assembly use the following command:
   ```
   pannopi -a /path/to/assembly.fasta -m /path/to/mtdna.fasta -t 32 -o /path/to/outdir
   ```
