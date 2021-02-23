library(ggplot2)
library(ggsci) 
library(tidyr)
library(dplyr)
library(ggrepel)

idlist=read.table('ID_category.tsv', header=F, col.names=c('id', 'category'), sep='\t')


# homozygosity by descent, runs of homozygosity 
temphom =read.table("memorial_exome.hbd", header=F , sep="\t",col.names=c('id1', 'hapid1' ,'id2', 'hapid2' , 'chr', 'start', 'end', 'lod', 'cM') )

hbd=merge(temphom, idlist, by.x='id1', by.y='id')

hbd %>% mutate(lenbp= end-start )%>%  group_by(id1, category) %>% summarize (totROH=sum(lenbp)) %>% ggplot(aes(id1,totROH/1000000 , fill=category))+geom_bar(stat='identity')+coord_flip()+theme_minimal()  +facet_wrap(category ~ ., scales='free_y' ) + ylab('genome in ROH (Mbp)') +xlab('')
ggsave('overview.png')



hbd$lenMb =(hbd$end-hbd$start)/1000000
hbd$rohtype =ifelse(hbd$lenMb<=0.5, 'Short (<= 0.5Mb)', ifelse(hbd$lenMb>5, 'Long (>5 Mb)', 'Intermediate (0.5Mb > ROH =< 5Mb)') ) 
 


##~~~  sum vs number of roh
sumNb= hbd %>% filter(category=='PREIMPLANTATION DEVELOPMENT ARREST') %>%  group_by(id1, rohtype ) %>% summarize ( sumROHs = sum(lenMb) , nbROHs= length(lenMb)) 
sumNb$rohtype <- factor(sumNb$rohtype, levels=c('Short (<= 0.5Mb)', 'Intermediate (0.5Mb > ROH =< 5Mb)', 'Long (>5 Mb)'))
ggplot(sumNb, aes(sumROHs , nbROHs , color = rohtype)  ) +geom_point() + theme_bw() +geom_text_repel(aes(label=id1), nudge_x = 0.25, nudge_y = 0.25) +xlab('Sum of ROHs (Mb)') + ylab('Number of ROHs') + scale_color_aaas()+ theme(legend.position="none")  + facet_wrap(rohtype ~ . , scales='free')
ggsave('sumvsnumber.png')

# FROH genomic inbreedding coefficient PMID:18760389  doi: 10.1016/j.ajhg.2008.08.007
#FROHis the fraction of each genome in ROH >1.5 Mb. doi.org/10.1038/s41467-019-12283

########################

#identity by descent 

tempibd=read.table('memorial_exome.ibd',header=F , sep="\t",col.names=c('id1', 'hapid1' ,'id2', 'hapid2' , 'chr', 'start', 'end', 'lod', 'cM') )

tempibd1=merge(tempibd, idlist, by.x='id1', by.y='id')
ibd=merge(temibd1, idlist, by.x='id2', by.y='id')
ibd$lenMb =(hbd$end-hbd$start)/10000000

ibd %>% filter (category.x=='PREIMPLANTATION DEVELOPMENT ARREST' & category.y=='PREIMPLANTATION DEVELOPMENT ARREST') %>% mutate (pair=paste(id1, hapid1, id2,  hapid2, sep='_')) %>% ggplot(aes(pair, lenMb, fill=lod)) + geom_bar(stat='identity') + facet_wrap(chr ~  . , scales ='free_y') +coord_flip() + scale_fill_gsea() + theme_bw()
ggsave('ibd.png')


 m4 %>% group_by(gene_symbol) %>% summarize(tot=sum(HIGH,LOW,MODERATE, na.rm=T  ))%>% filter (tot>20 ) %>% ggplot(aes(reorder(gene_symbol, tot ), tot)) +geom_bar(stat='identity')+coord_flip()


 m4 %>% group_by(gene_symbol) %>% summarize(tot=sum(HIGH,LOW,MODERATE, na.rm=T  ))%>% filter (tot>20 ) %>% ggplot(aes(reorder(gene_symbol, tot ), tot)) +geom_bar(stat='identity')+coord_flip()
