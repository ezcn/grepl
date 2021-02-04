library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify two arguments: <inputFile> <uniqueCodeForResults>', call.=FALSE)
} 

code=args[2]
mydata=read.table(args[1], header=T, sep='\t')
#myf <- mydata %>% filter(pLIscore >= 0.9 & caddPercent >=0.9 | sumGene >=2) %>% filter(GrandMean <= 0.05) %>% filter(IMPACT != "LOW")
myf <- mydata %>% filter(GrandMean <= 0.05) %>% filter(IMPACT != "LOW")
myf$impact= factor(myf$IMPACT,levels=c( 'MODERATE','HIGH'))
#maxv <- myf %>% group_by(index_x, SYMBOL) %>% count() %>% filter(n<=args[3])
maxv <- myf %>% select(index_x, SYMBOL) %>% distinct() %>% group_by( SYMBOL) %>% tally() %>% filter(n<=5)
myd <-merge(myf, maxv)
#myd<-myf
write.table(myd, "totally_filtered_variants.tsv", sep = "\t", quote = F, row.names=F, col.names=T)

#### summary stats 
sink(paste(code, '.output.txt', sep=''))
totnumberofuniqevariants=myd %>% select(index_x) %>% distinct() %>% tally() 
print(paste(code, 'Tot Number of Uniqe Variants', totnumberofuniqevariants, sep='\t'))
totalnumberuniqueNOB=myd %>% filter(rare05=='NOB') %>% select (index_x) %>% distinct() %>% tally()
print(paste(code, 'Tot Number of NOB',totalnumberuniqueNOB, sep='\t'))

impdf<-myd %>% select (index_x , IMPACT ) %>% distinct() %>% group_by(IMPACT) %>% tally()
print(impdf)


avsdSites<- myd %>% select(index_x, ID) %>% distinct() %>% group_by(ID) %>% tally() %>% summarize(mean(n), sd(n))
print(paste(code, 'Average number of sites per sample', avsdSites[1], sep='\t'))
print(paste(code, 'Sd number of sites per sample', avsdSites[2], sep='\t'))
totnumberofuniqegenes=myd %>% select(Gene) %>% distinct() %>% tally()
print(paste(code, 'Total number of Unique Genes', totnumberofuniqegenes, sep='\t'))
totnumberofuniqetranscripts=myd %>% select(Feature) %>% distinct() %>% tally()
print(paste(code, 'Total number of Unique Transcripts', totnumberofuniqetranscripts, sep ='\t'))
avesdGenes<- myd %>% select(Gene, ID) %>% distinct() %>% group_by(ID) %>% tally() %>% summarize(mean(n), sd(n))
print(paste(code, 'Average number of genes per sample', avesdGenes[1], sep='\t'))
print(paste(code, 'Sd number of genes per sample', avesdGenes[2], sep='\t'))
avesdTranscripts<- myd %>% select(Feature, ID) %>% distinct() %>% group_by(ID) %>% tally() %>% summarize(mean(n), sd(n))
print(paste(code, 'Average number of transcripts per sample', avesdTranscripts[1], sep='\t'))
print(paste(code, 'Sd number of transcripts per sample', avesdTranscripts[2], sep='\t'))
sink()


####### 1. remove genes with too many variants 
### Multiple variants in the same gene

print('multiple variants same gene')
myd %>% select(index_x,ID, Gene,IMPACT ) %>% distinct() %>% group_by( ID, Gene,IMPACT) %>% tally() %>% filter(n>1) %>% spread(IMPACT, n) %>% write.table(paste(code, '.multiplevariants_pergene.tsv'), sep="\t", quote=F, row.names=F, col.names=T)

mycol=c("#49BEB7", "#FF5959")
myd %>% select(index_x,ID, SYMBOL,IMPACT) %>% distinct() %>% group_by( ID, SYMBOL,IMPACT) %>% tally() %>% filter(n>1) %>% ggplot(aes(reorder(SYMBOL, n) , n , fill=IMPACT ))  +  geom_bar(stat='identity')+coord_flip() +theme_minimal() +scale_fill_manual(values=mycol)
ggsave('variantspergenemore2.png')

#myd %>% select(index_x,sample, gene_symbol, impact) %>% distinct() %>% group_by( gene_symbol, impact) %>% tally() %>% spread(impact, n) %>% summarize(tot=sum(MODERATE, HIGH, na.rm=T )) %>% ggplot(aes(tot)) +geom_density()+theme_minimal()+xlab('Number of variants per gene') + geom_vline( xintercept=quantile(tot, .95)) 


myg =myd %>% select(index_x,ID, Gene,IMPACT) %>% distinct() %>% group_by( Gene,IMPACT) %>% tally() %>% spread(IMPACT, n) %>% summarize(tot=sum(MODERATE, HIGH, na.rm=T )) 

dens <- density(myg$tot)
df <- data.frame(x=dens$x, y=dens$y)
probs <- c(0,  0.95, 0.99, 1)
quantiles <- quantile(myg$tot, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))
ggplot(df, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant))
#scale_x_continuous(breaks=quantiles) 
#+ scale_fill_brewer(palette="Reds", guide="none") +theme_minimal()+labs(x='Number of variants per gene' , y='Density')
ggsave('variantspergenedensity.png')



####### 3. Results 
#######
print('~~ number')
#myd %>% select(sample, index_x, impact ) %>% distinct() %>% group_by (sample, impact) %>% tally() %>% spread(impact, n)
mycol=c( "#FF5959","#49BEB7")
number=myd %>% select(ID, index_x,IMPACT ) %>% distinct() %>% group_by (ID, IMPACT) %>% tally() 
#number$impact= factor(number$impact,levels=c('HIGH','MODERATE'))
ggplot(number, aes(as.factor(ID), n, fill=IMPACT)) + geom_bar(stat='identity') + ggtitle(paste (code,'- Unique variants per sample', sep=' ')) + theme_minimal() + scale_fill_manual(values=mycol)+coord_flip() +xlab('') + ylab('Number of unique variants') 
ggsave(paste(code, '_number.png', sep=''))

######################
print('sumgene')
myd %>% select(ID, Gene, sumGene)%>% distinct () %>% group_by(ID, sumGene)  %>% tally() %>% ggplot (aes(as.factor(ID), n, fill=sumGene)) + geom_bar(stat='identity')+ theme_minimal()  +coord_flip() +scale_fill_gradient(low = "#132B43",  high = "#56B1F7") +labs(fill = "Number of lists", y='Number of unique genes', x='' )  #+ ggtitle(paste(code,'- Genes in lists', sep=' '))             
ggsave(paste (code, '_sumGene.png', sep=''))  


#####################
print('most_severe_consequence')
csqtype=myd %>% select(ID, index_x, IMPACT, Consequence) %>% distinct() %>% group_by (ID, IMPACT, Consequence) %>% tally() 
#number$impact= factor(number$impact,levels=c('HIGH','MODERATE'))
nbcolors=length(levels(myd$Consequence)) 
mycolors <- colorRampPalette(brewer.pal(8,'Set2'))(nbcolors)
ggplot(csqtype, aes(as.factor(ID), n, fill=Consequence)) + geom_bar(stat='identity')  + theme_bw()+ facet_wrap (IMPACT ~ ., scales='free')+ scale_fill_manual(values = mycolors) +coord_flip() +labs(fill = "Consequences", y='Number of unique variants', x='' ) #+ ggtitle(paste(code, '- most_severe_consequence', sep=' '))
ggsave(paste(code, '_most_severe_consequence.png', sep=''))



######
print('geneTranscripts')
genes =myd %>% select(ID, SYMBOL ) %>% group_by(ID) %>% distinct() %>% tally()
transcripts = myd %>% select(ID, Feature ) %>% group_by(ID) %>% distinct() %>% tally()
geneTranscripts=full_join(genes, transcripts, by='ID')
colnames(geneTranscripts)=c('sample','genes','transcripts')
geneTranscripts %>% gather ('Feature_type','n', c(genes, transcripts)) %>% ggplot(aes(as.factor(sample) , n , color=Feature_type) ) + geom_point() + ggtitle(paste (code, '- Genes and transcripts', sep=' ')) + theme_bw() + scale_color_nejm()
ggsave(paste (code, '_geneTranscripts.png', sep='')) 


########################
print('cadd pli ')
myd %>% ggplot(aes(CADD_RAW, color=IMPACT)) + geom_density() + ggtitle(paste(code,'- CADD', sep=' ')) + theme_bw() + scale_color_nejm()
ggsave(paste(code, '_CADD.png', sep=''))

myd %>% ggplot(aes(pLIscore, color=IMPACT)) + geom_density() + ggtitle(paste(code,'- pLI', sep=' ')) + theme_bw() + scale_color_nejm()
ggsave(paste(code, '_pLI.png', sep=''))


############ 
print('plicadsumgene')
myd %>% ggplot(aes(pLIscore, CADD_RAW, color=as.factor(sumGene))) + geom_point() + ggtitle(paste(code, '- pLI CADD sumGene', sep=' ')) + theme_bw() + scale_color_futurama()
ggsave(paste(code, '_pLiCADDsumGene.png', sep='')) 


##########################
print('plicadd')
#mycol=c('#ffbd69', '#fe346e', '#b21f66', '#381460') 
mycol=c('#FFBD69', '#2EF276','#FE346E', '#B21F66', '#381460')
datawithpliandcadd=myd %>% select (pLIscore, CADD_RAW, sumGene, IMPACT) %>% distinct()%>% filter (!is.na(pLIscore) & !is.na(CADD_RAW)) %>% tally() 
myd %>% select (ID, pLIscore, CADD_RAW, sumGene, IMPACT) %>% distinct() %>% ggplot(aes(pLIscore, CADD_RAW, color=as.factor(sumGene), alpha =0.7 )) + geom_point() + ggtitle(paste (code, 'pLI, CADD, gene lists \n# unique sites with data =',datawithpliandcadd, sep='')) + theme_bw() + scale_color_manual (values=mycol) + facet_wrap(ID ~ .)
ggsave(paste(code, '_plicadd.png', sep='')) 



### Aggregate analysis 
print('aggregate')
mycol=c( "#FF5959", "#49BEB7")

myd %>% select (IMPACT, SYMBOL, ID) %>%distinct() %>% group_by( IMPACT, SYMBOL)  %>% tally() %>% filter(n>1) %>% ggplot( aes(reorder(SYMBOL, n), n, fill=IMPACT)) + geom_bar(stat='identity', position='dodge') + xlab('genes') + ylab('number of samples') + scale_fill_manual(values=mycol) + coord_flip() + theme_minimal() + facet_wrap (IMPACT ~ ., scales ='free') + ggtitle(paste(code, '- Genes shared among samples', sep=' ')) 
ggsave(paste(code, '_aggregate.png', sep='')) 


myd %>% select (index_x, IMPACT, SYMBOL, ID) %>%distinct() %>% group_by( IMPACT, SYMBOL)  %>% tally() %>% filter(n>1) %>% ggplot( aes(reorder(SYMBOL, n), n, fill=IMPACT)) + geom_bar(stat='identity', position='dodge') + xlab('genes') + ylab('number of samples') + scale_fill_manual(values=mycol) + coord_flip() + theme_minimal() + facet_wrap (IMPACT ~ ., scales ='free') + ggtitle(paste(code, '- Genes shared among samples', sep=' ')) 
ggsave(paste(code, '_aggregateindex_X.png', sep='')) 

#######
print('aggregate impact alt count')
sh<-myd %>% select (index_x, IMPACT, SYMBOL, ID, ALTcount) %>%distinct() %>% group_by( IMPACT, SYMBOL)  %>% count() %>% rename(sharedVariants = n) %>% filter(sharedVariants>1)
mydsh <- merge(myd,sh)
mydsh$ID <- as.factor(mydsh$ID)
#mycol = c("#e2979c","#e7305b", "#f4ebc1", "#a0c1b8")
mycol = c("#e2979c", "#f4ebc1", "#a0c1b8")
mydsh %>% select(ID, SYMBOL, ALTcount, IMPACT)  %>% mutate(colorTile = ifelse (IMPACT == "HIGH" & ALTcount == 2 ,"highHom" , ifelse(IMPACT == "HIGH" & ALTcount ==1 , "highHet", ifelse(IMPACT == "MODERATE" & ALTcount == 2, "ModHom", "ModHet")))) %>% ggplot(aes(ID, SYMBOL, fill = as.factor(colorTile))) + geom_tile() + scale_fill_manual(values= mycol) + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1)) + xlab("") + ylab("") + theme(legend.title=element_blank())
ggsave(paste(code, '_aggregate_impact_count.png', sep=''),width = 10, height = 15)

##### not shared genes
mycol = c("#e2979c","#e7305b", "#f4ebc1", "#a0c1b8")
Nsh<-myd %>% select (index_x, IMPACT, SYMBOL, ID, ALTcount) %>%distinct() %>% group_by( IMPACT, SYMBOL)  %>% count() %>% rename(notShared = n) %>% filter(notShared==1)
mydnsh<-merge(myd,Nsh)
mydnsh %>% select(index_x, IMPACT, SYMBOL, ID, ALTcount) %>% distinct() %>% mutate(colorTile = ifelse (IMPACT == "HIGH" & ALTcount == 2 ,"highHom" , ifelse(IMPACT =="HIGH" & ALTcount ==1 , "highHet", ifelse(IMPACT == "MODERATE" & ALTcount == 2, "ModHom", "ModHet")))) %>% ggplot(aes(SYMBOL, fill= as.factor(colorTile))) +geom_bar() + coord_flip() + facet_wrap( .~ ID ,scale= "free", ncol = 10) + scale_fill_manual(values = mycol) + ylim(0,3) + theme_bw() + xlab("") + ylab("") + theme(legend.title=element_blank())
ggsave(paste(code, '_individual_impact_count.png', sep='') , width = 15, height = 20)

#### all genes per sample
mycol = c("#e2979c","#e7305b", "#f4ebc1", "#a0c1b8")
allg<-myd %>% select (index_x, IMPACT, SYMBOL, ID, ALTcount) %>%distinct() %>% group_by( IMPACT, SYMBOL)  %>% count() %>% rename(allGenes = n)
mydallg<-merge(myd,allg)
mydallg %>% select(index_x, IMPACT, SYMBOL, ID, ALTcount) %>% distinct() %>% mutate(colorTile = ifelse (IMPACT == "HIGH" & ALTcount == 2 ,"highHom" , ifelse(IMPACT =="HIGH" & ALTcount ==1 , "highHet", ifelse(IMPACT == "MODERATE" & ALTcount == 2, "ModHom", "ModHet")))) %>% ggplot(aes(SYMBOL, fill= as.factor(colorTile))) +geom_bar() + coord_flip() + facet_wrap( .~ ID ,scale= "free", ncol = 10) + scale_fill_manual(values = mycol) + ylim(0,3) + theme_bw() + xlab("") + ylab("") + theme(legend.title=element_blank())
ggsave(paste(code, '_allGenes_impact_count.png', sep='') , width = 15, height = 20)