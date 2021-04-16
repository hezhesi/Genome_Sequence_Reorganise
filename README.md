# Genome Sequence Reorganise

`Genome_Sequence_Reorganise.R` is an R script, to generate the revised genome sequence based on a standardized formatted table for correcting misplaced segments of genomes.

After GOGGs analysis and when anomalies were detected, we manually generated a standardized formatted table for correcting misplaced segments of genomes, with correct placement based on where in the genome the pattern of alleles present in each line of the mapping population best matches. Fine placement (i.e., between precisely which two gene models the segment should be inserted) was determined to preserve collinearity with the A. thaliana genome, based on top BLAST similarity to A. thaliana CDS gene models. 

These tables (exmample shown as in `Chiifu-NI100-TO1000_edits.xlsx`)and the original genome resource files were then taken into `Genome_Sequence_Reorganise.R`, to generate the revised genome sequence and annotation files. 
