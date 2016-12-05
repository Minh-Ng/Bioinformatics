The following need to be installed in order to work properly:
* [fastqc] - fastq quality checking tool
* [STAR] - rna-seq alignment tool
* [cufflinks] - outputs gene structure information from rna-seq reads
* [samtools] - sort SAM file and output to BAM or SAM

Useful:
* [seqtk] - convert fastq into fasta
* [cuffdiff] - compare cufflinks output
 
## Running the Pipeline
Put the inputs (reference fasta and rna-seq read fastq) inside of the inputs folder  
Run `chmod 755 pipeline.sh` in order to give read/execution permissions to others  
Run `./pipeline.sh –c` in order to clear the outputs of the previous run.  
Run `./pipeline.sh` or `sh pipeline.sh` in order to run the pipeline.  

## Results
The fastqc outputs can be seen in the `fastqc_results` folder.   
The cufflinks results can be seen in the `cufflinks_results` folder.  
The star results will be outputted in the folder that pipeline.sh resides.  
This pipeline can take a long time as STAR can take a long time to run especially on large sets of data.

[fastqc]: http://www.bioinformatics.babraham.ac.uk/projects/download.html
[STAR]: https://github.com/alexdobin/STAR
[cufflinks]: http://cole-trapnell-lab.github.io/cufflinks/install/
[samtools]: http://samtools.sourceforge.net/
[seqtk]: https://github.com/lh3/seqtk
[cuffdiff]: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/