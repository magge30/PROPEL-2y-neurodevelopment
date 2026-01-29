# PROPEL - 2y neurodevelopment 

Scripts to produce the analysis and figures for the PROPEL 2y neurodevelopment paper: <br>
Dynamics of early gut microbiota maturation in extremely preterm infants are associated with neurodevelopment at 2 years <br>
preprint: ttps://doi.org/10.1101/2025.03.17.25324095

## Data
- Sequences available at European Nucleotide Archive : PRJEB36531. <br>
- metadata, OTU tabel and ASV table in .csv format <br>
- phyloseq obejct (data.rds) integrating the .cvs data ready to run with the Rscirpts  <br>

## Bioinformatics <br>
From Raw sequence data to ASV table (bbduk, dada2 and ASV table preprocessin) <br>
Scripts: https://github.com/magge30/PROPEL-ELBW-16S <br>

## Biostatistics: <br>
Generalised R scripts for the statisical anlysis:
- SparseMCMM <br>
- Alpha-diversity <br>
- Beta-diversity <br>
- PLS-DA <br>
- Longitudinal analys <br>
