### summary statistics from mapping - (from bam files) 


# Summary statistics of variant calling - chr 22 (form vcf files) [script.R](statsRscript/SNstatsVCF.R)

### [Number of MNPs](img/SN-numberofMNPs.png) number of rows with a MNP, such as CC>TT
### [Number of SNPs](img/SN-numberofSNPs.png) number of rows with a SNP
### [Number of Indels](img/SN-numberofindels.png) number of rows with an indel
### [Number of multiallelic SNPs sites](img/SN-numberofmultiallelicSNPsites.png) number of rows with multiple alternate alleles, all SNPs
### [Number of multiallelic sites](img/SN-numberofmultiallelicsites.png) number of rows with multiple alternate alleles
### [Number of no-ALTs](img/SN-numberofno-ALTs.png) reference-only sites, ALT is either "." or identical to REF
### [Number of Others](img/SN-numberofothers.png) number of rows with other type, for example a symbolic allele or a complex substitution, such as ACT>TCGA
### [Number of Records](img/SN-numberofrecords.png) number of data rows in the VCF

# Substitution types - chr 22 (from vcf files) [script.R](statsRscript/STstatsVCF.R)

### [Substitution types](img/ST-Substitutiontypes.png) 

# Transitions/transversions - chr 22 (from vcf files) [script.R](statsRscript/TSTVstatsVCF.R)

### [Transitions](img/transition.png)
### [Transversion](img/transversion.png)
### [TSTV ratio](img/tstvratio.png)

# Depth - chr 22 (from vcf files) [script.R](statsRscript/DEPTHstatsVCF.R)

### [Depth](img/depth-quality1.png)  number of reads, on average, thet are likely to be aligned at a given reference base position.
 