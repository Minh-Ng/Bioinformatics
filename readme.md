Put the inputs (reference fasta and rna-seq read fastq) inside of the inputs folder
Run the following command in order to give read/execution permissions to others
	`chmod 755 pipeline.sh`
Run `./pipeline.sh –c` in order to clear the outputs of the previous run
Run `./pipeline.sh` or `sh pipeline.sh` in order to run the pipeline.
The fastqc outputs can be seen in the `fastqc_results` folder. 
The cufflinks results can be seen in the `cufflinks_results` folder.
The star results will be outputted in the folder that pipeline.sh resides.
This pipeline can take a long time as STAR can take a long time to run especially on large sets of data.

