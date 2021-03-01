#### Required tools and programming language
- bwa (http://bio-bwa.sourceforge.net/bwa.shtml)
- sambamba (https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
- freebayes (https://github.com/freebayes/freebayes)
- samtools/htslib (https://github.com/samtools/htslib)
- bcftools(http://samtools.github.io/bcftools/bcftools.html)
- Vt (https://genome.sph.umich.edu/wiki/Vt)
- vep (https://www.ensembl.org/info/docs/tools/vep/index.html)
- vcftools (https://vcftools.github.io/index.html)
- Python3
- R
#### 1. Align reads to reference genome 
```
bwa mem -t 24 -R "@RG\tID:$id\tSM:$id" reference.GRCh38.fa sample_R1.fastq.gz sample_R2.fastq.gz | samtools view -b - > sample.raw.bam && touch sample.align_ok
```
- sort bed file

```
sambamba sort -t 16 -m 32G --tmpdir /scratch -o sample.raw.sorted.bam sample.raw.bam
```
- remove PCR duplicates
```
sambamba markdup -t 8 --tmpdir ~/tmp --overflow-list-size 500000 sample.raw.sorted.bam sample.bam
```
#### 2. Variant Calling by sample and chromosome
```
freebayes -f reference.GRCh38.fa -r $chr -g 500 -b sample.bam > sample.$chr.fb.vcf && touch sample.$chr.fb_ok
```
-bgzip and tabix
```
bgzip sample.$chr.fb.vcf
```
```
tabix -p vcf sample.$chr.fb.vcf.gz
```
- filter vcf for quality > 20
```
bcftools filter -i "QUAL>20"  sample.$chr.fb.vcf.gz -O z -o sample.$chr.fb.filtered.vcf.gz && touch sample.$chr_filter_ok
```
- tabix
- normalize filtered vcf
```
vt normalize -n -r GRCh38.fa sample.$chr.filtered.fb.vcf.gz > sample.$chr.fb.norm.vcf && touch sample.$chr_norm_ok
```
- bgzip and tabix
- decompose normalized vcf
```
vt  decompose_blocksub sample.$chr.fb.norm.vcf.gz > sample.$chr.fb.decomp.vcf && touch sample.$chr_decomp_ok
```
-bgzip and tabix

#### 3. Merge all samples
```
bcftools merge -O z -0 -o testdata/sam/samples.chr22.vcf.gz sample1.$chr.fb.decomp.vcf.gz  sample2.$chr.fb.decomp.vcf.gz sample3.$chr.fb.decomp.vcf.gz sampleX.$chr.fb.decomp.vcf.gz 
```
- bgzip and tabix

#### 4.  Annotate variants using `V.E.P`
merged vcf per chromosome. output is a table

testdata/sam/samples.chr22.vcf.gz -> testdata/sam/samples.chr22.vep.tsv
testdata/con/controls.chr22.vcf -> testdata/con/controls.chr22.vep.tsv  
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/data/

```
vep --af_1kg --af_gnomad --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --fork 4 --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --dir_cache /data/biocontainers/vepcache --offline --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,AF,AFR_AF,AMR_AF,ASN_AF,EUR_AF,EAS_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,CADD_RAW,CADD_PHRED" --force_overwrite --variant_class -i allsamples.chrx.vcf --plugin CADD,/CADD/whole_genome_SNVs.tsv.gz -o  allsamples.$chr.vep.tsv && touch tabOk/allsamples.$chr.table_ok
```

#### 4.1  Remove Header of file VEP
for i in {1..22} X ; do grep -v "##" ${i}.vep.tsv > ${i}.vep.noheader.tsv; done
for i in {1..22} X; do sed -i '/Uploaded_variatio/s/^#//g' ${i}.vep.noheader.tsv ; done

#### 5. Parsing VEP annotations with [grep.csq.sql.py](https://github.com/ezcn/grep/tree/master/6_grepPipeline/scr/grep.csq.sql.py)
 per chromosome
 ```
 Options :
 -db = path to new sql db
 -i = path to input vep.noheader.tsv table file
 -g = path to gene list file
 -pli = path to pLI score table
 -o = path to output file
 ```
 
```
command line example
python3 /grepPipeline/scr/grep.csq.sql.py -db samples.chr22.db -i vep/samples.chr22.vep.tsv -g testdata/db/all_geneList.tsv -pli testdata/db/pLIscore.tsv -fathmmcoding testdata/db/fathmm_xf_coding.chr22.tsv -fathmmnc testdata/db/fathmm_xf_noncoding.$(chr).tsv -o samples.chr22.sql.csq
```

#### 6. Extract individual's allele count from vcf 
1. Extract counts with vcftools 
```
vcftools --gzvcf samples.chr22.vcf.gz --out $(chr)/$(id).chr22_counts  --counts --indv $(id)
```
2. Change file formatting with [altCounts.py](https://github.com/ezcn/grep/tree/master/6_grepPipeline/scr/altCounts.py)
```
-i = path to input file
-id = sample ID
python3 altCounts.py -i /$(chr)/$(id).$(chr)_counts.frq.count -id $(id)
```

#### 7. Filter variable sites with [grep.filter.all.slq.py](https://github.com/ezcn/grep/tree/master/6_grepPipeline/scr/grep.filter.all.slq.py)
per chromosome 
```
Options :
-dbC = path to control db 
-dbS = path to control db
-cl = list of controls IDs
-sl = list of sampes IDs
-ff = threshold for first frequency definition
-f = threshold for second frequency definition
-type = feature_type(genic , intergenic, regularoty)
-r = value for filtering rare variants (not equal to)
-pli = threshold for  pLI score
-cadd = treshold for CADD score
-g = number of gene lists
-n = number of control individual to sample
-i = number of iterations
-pathTodirCtrl = path to directory with controls allele counts files
-ctgn = path to control genes file
-gtd = path to discarded genes file
-pathTodir = path to directory with samples allele counts files
-chrom = chromosome
-ac =  allele count treshold
-gt = threshold for gene frequency in control population
-o =  path to output file
 ```

```
command line example : 
python3 scr/grep.filter.all.slq.py -db controls.chr22.db -dbS samples.chr22.db -cl testdata/con/controls_id.txt -sl testdata/sam/samples_id.txt -ff 0.01 -f 0.05 -type Transcript -r false -pli 0.7 -cadd 0.5 -g 2 -n 10 -i 100 -pathTodirCtrl /control/alleleCount/chr22 -ctgn controlGenes.chr22.tsv -gtd GenesToDiscardHgdp.chr22.tsv -pathTodir /samples/counts/chr22 -chrom chr22 -ac 1 -gt 0.5 -o samples.chr22.filtered.tsv
```
- concat all chromosomes -> allSamples.filtered.tsv

#### 8. Format results with [grep.curation_nolow.R](https://github.com/ezcn/grep/tree/master/6_grepPipeline/scr/grep.curation_nolow.R) 
merge all chromosomes samples files
```
Options : 
filtered samples file
plot and tables names prefix

command line example :
Rscript grepPipeline/scr/grep.curation_nolow.R allSamples.filtered.tsv Grep_allsamples
```
