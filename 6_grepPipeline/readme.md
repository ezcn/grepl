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
bcftools merge -O z -0 -o samples.vcf.gz sample1.$chr.fb.decomp.vcf.gz  sample2.$chr.fb.decomp.vcf.gz sample3.$chr.fb.decomp.vcf.gz sampleX.$chr.fb.decomp.vcf.gz 
```
- bgzip and tabix

#### 4.  Annotate variants using `V.E.P`
merged vcf per chromosome. output is a table

allsamples.$chr.vcf -> allsamples.$chr.vep.tsv
all.control.chrx.vcf -> allcontrols.$chr.vep.tsv  
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/data/

```
vep --af_1kg --af_gnomad --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --fork 4 --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --dir_cache /data/biocontainers/vepcache --offline --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,AF,AFR_AF,AMR_AF,ASN_AF,EUR_AF,EAS_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,CADD_RAW,CADD_PHRED" --force_overwrite --variant_class -i allsamples.chrx.vcf --plugin CADD,/CADD/whole_genome_SNVs.tsv.gz -o  allsamples.$chr.vep.tsv && touch tabOk/allsamples.$chr.table_ok
```

#### 4.1  Remove Header of file VEP
for i in {1..22} X ; do grep -v "##" ${i}.vep.tsv > ${i}.vep.noheader.tsv; done
for i in {1..22} X; do sed -i '/Uploaded_variatio/s/^#//g' ${i}.vep.noheader.tsv ; done

#### 5. Format data base files per chromosome 
- merge gene lists: merge_allGenes.py  
Input gene list: tab-delimited two column file with ENSGENEID, category. Genes in more than one category can be repeated 
Output: merged file 
- fathmm 

#### 6. Parsing VEP annotations with [grep.csq.sql.py](https://github.com/ezcn/grep/tree/master/grepPipeline/scr/grep.csq.sql.py)` 
 per chromosome
```
python3 /grepPipeline/scr/grep.csq.sql.py -db /db/sample.chrx.db -i vep/all.sample.chrx.vep.tsv -g testdata/db/all_geneList.tsv -pli testdata/db/pLIscore.tsv -fathmmcoding /lustre/home/enza/fathmm_database/coding/fathmm_xf_coding.$(chr).tsv -fathmmnc /lustre/home/enza/fathmm_database/noncoding/fathmm_xf_noncoding.$(chr).tsv -o /lustre/home/enza/grep/csq/all.sample.chrx.sql.csq
```

#### 7. Extract individual's allele count from vcf 
1. Extract counts with vcftools 
```
vcftools --gzvcf all.sample.chrx.vcf.gz --out $(chr)/$(id).$(chr)_counts  --counts --indv $(id)
```
2. Change file formatting with `altCounts.py`
```
python3 grepPipeline/scr/altCounts.py -i /$(chr)/$(id).$(chr)_counts.frq.count -id $(id)
```

#### 8. Filter variable sites with [grep.filter.all.slq.py](https://github.com/ezcn/grep/tree/master/grepPipeline/scr/grep.filter.all.slq.py)`
per chromosome 
```
python3 grep.filter.all.slq.py -db /control/db/hgdp.$(chr).full.db -dbS /samples/db/allGrep.$(chr).db -cl control.list.txt -sl samples.list.txt -ff 0.01 -f 0.05 -type Transcript -r false -pli 0.7 -cadd 0.5 -g 2 -n 10 -i 100 -pathTodirCtrl /control/alleleCount/$(chr) -ctgn controlGenes.$(chr).tsv -gtd GenesToDiscardHgdp.$(chr).tsv -pathTodir /samples/counts/$(chr) -chrom $(chr) -ac 1 -gt 0.5 -o samples.$(chr).filtered.tsv
```
- concat all chromosomes -> allSamples.filtered.tsv

#### 9. Format results with [grep.curation_nolow.R](https://github.com/ezcn/grep/tree/master/grepPipeline/scr/grep.curation_nolow.R) 

```
Rscript grepPipeline/scr/grep.curation_nolow.R allSamples.filtered.tsv Grep_allsamples
```
