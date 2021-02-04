library(dplyr)
library(argparse)
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--countfile", action="store")
parser$add_argument("--numberCases", action="store")
parser$add_argument("--numberControls", action="store_true")
args <- parser$parse_args()
# lp95obs --> observed -log10(p-value) @95% of all genes
# lp95exp --> expected -log10(p-value) @95% of all genes
# lp0obs --> always zero, observed p-value of the gene with the max expected -log10(p-value) among genes with p-value==1
# lp0exp --> expected p-value of the gene with the max expected -log10(p-value) among genes with p-value==1
# INPUT: table og genes and counts
# ca1Het --> count of cases with at least one heterozygous qualifying variant in the gene
# co1Het --> count of controls with at least one heterozygous qualifying variant in the gene
# ca1Hom --> count of cases with at least one homozygous qualifying varinat in the gene
# co1Hom --> count of controls with at least one homozygous qualifying variant in the gene 

totcase=numberCases
totcon=numberControls
# read count file header: gene, ca1Het, ca1Hom, co1Het, co1Hom
pdaHet<-provapda %>% select(SYMBOL, ID) %>% distinct() %>%  group_by(SYMBOL) %>% tally() %>% rename(ca1Het = n)
pdaHom<-provapda %>% filter(ALTcount == 2) %>% select(SYMBOL, ID) %>% distinct() %>%  group_by(SYMBOL) %>% tally() %>% rename(ca1Hom = n)
pdaAll<- merge(pdaHet, pdaHom, by = "SYMBOL", all = T) %>% replace(is.na(.), 0)

myd=read.table(countfile, header=T )


######### pvalues recessive and dominant (non vanno bene, spostare nelle righe successive al posto giusto)
#mymatDom<-matrix(c(ca1Het, numberCases-ca1Het, co1Het, numberControls-co1Het), ncol=2, nrow=2)
#mymatRec<-matrix(c(ca1Hom, numberCases-ca1Hom, co1Hom, numberControls-co1Hom), ncol=2, nrow=2)

#fisher 
myd<- myd %>%  rowwise() %>% mutate(fpdom=fisher.test(matrix(c(ca1Het, numberCases-ca1Het, co1Het, numberControls-co1Het), ncol=2, nrow=2))$p.value, fstatdom=fisher.test(matrix(c(ca1Het, numberCases-ca1Het, co1Het, numberControls-co1Het), ncol=2, nrow=2))$estimate, fprec=fisher.test(matrix(c(ca1Hom, numberCases-ca1Hom, co1Hom, numberControls-co1Hom), ncol=2, nrow=2))$p.value, fstatrec=fisher.test(matrix(c(ca1Hom, numberCases-ca1Hom, co1Hom, numberControls-co1Hom), ncol=2, nrow=2))$estimate, log10fpdom=-log10(fpdom),  log10fprec=-log10(fprec) )
#chisquare
myd<- myd %>%  rowwise() %>% mutate(cpdom=chisq.test(matrix(c(ca1Het, numberCases-ca1Het, co1Het, numberControls-co1Het), ncol=2, nrow=2))$p.value, cstatdom=chisq.test(matrix(c(ca1Het, numberCases-ca1Het, co1Het, numberControls-co1Het), ncol=2, nrow=2))$statistic, cprec=chisq.test(matrix(c(ca1Hom, numberCases-ca1Hom, co1Hom, numberControls-co1Hom), ncol=2, nrow=2))$p.value, cstatrec=chisq.test(matrix(c(ca1Hom, numberCases-ca1Hom, co1Hom, numberControls-co1Hom), ncol=2, nrow=2))$statistic, log10cpdom=-log10(cpdom),  log10cprec=-log10(cprec) )  

	
######## calulate expected pval
#dominant fisher 
myd<-myd[order(myd$log10fpdom),]
myd$log10fpdomexp<-sort(-log10(runif(nrow(myd))))  ## 'nrows' random values from uniform distribution, sorted

#dominant chisq 
myd<-myd[order(myd$log10cpdom),]
myd$log10cpdomexp<-sort(-log10(runif(nrow(myd))))  ## 'nrows' random values from uniform distribution, sorted

#recessive fisher 
myd<-myd[order(myd$log10fprec),]
myd$log10fprecexp<-sort(-log10(runif(nrow(myd))))

#recessive chisq 
myd<-myd[order(myd$log10cprec),]
myd$log10cprecexp<-sort(-log10(runif(nrow(myd))))


######## calculate lambda 
#lambda95 fisher dominant  (fd)
fdlp95obs <- quantile( subset(myd, log10fpdom!=0)$log10fpdom, 0.95) #p95
fdlp95exp <- quantile( subset(myd, log10fpdom!=0)$log10fpdomexp, 0.95) #u95
fdlp0exp <- maxlogp1 <- max(subset(myd, fpdom==1)$log10fpdomexp) #u0
fdlp0obs <-0 #p0

fdlambda95 <-(fdlp95obs[[1]]-fdlp0obs) / (fdlp95exp[[1]]-fdlp0exp) #slope

#lambda95 fisher recessive 
frlp95obs <- quantile( subset(myd, log10fprec)!=0)$log10fprec, 0.95) #p95
frlp95exp <- quantile( subset(myd, log10fprec!=0)$log10fprecexp, 0.95) #u95
frlp0exp <- maxlogp1 <- max(subset(myd, fprec==1)$log10fprecexp) #u0
frlp0obs <-0 #p0

frlambda95 <-(frlp95obs[[1]]-frlp0obs) / (frlp95exp[[1]]-frlp0exp) #slope

#lambda95 chisq dominant 
cpp95obs <- quantile( subset(myd, log10cpdom!=0)$log10cpdom, 0.95) #p95
cpp95exp <- quantile( subset(myd, log10cpdom!=0)$log10cpdomexp, 0.95) #u95
cpp0exp <- maxlogp1 <- max(subset(myd, cpdom==1)$log10cpdomexp) #u0
cpp0obs <-0 #p0

cplambda95 <-(cpp95obs[[1]]-cpp0obs) / (cpp95exp[[1]]-cpp0exp) #slope
#lambda95 chisq recessive 


#lambda fisher dominant 


> provapda %>% select(SYMBOL, ID) %>% distinct() %>%  group_by(SYMBOL) %>% tally() 



######## PLOT fisher dominant
yint<-fdlp95obs[[1]]-fdlambda95*fdlp95exp[[1]]

#yint<-flp95obs[[1]]-lambda95*lp95exp[[1]]
maxp<-ceiling(max(myd$log10fpdom))
myd %>% ggplot (aes(log10fpdomexp, log10fpdom) )+ geom_point() + geom_abline(intercept=yint, slope=fdlambda95, linetype=3, color='grey')+ geom_abline(intercept=0, slope=1, color='grey') + xlim(0, maxp)+ylim(0, maxp) + theme_minimal() +labs(x='expected -log10(p-value)' , y='observed -log10(p-value)' , title= paste('Dominant test - fdlambda95 =' ,round(fdlambda95, 3)  , sep= ' ' )  )
ggsave('dominant.png')

#### chisq dominant
yint<-cpp95obs[[1]]-cplambda95*cpp95exp[[1]]
maxp<-ceiling(max(myd$log10cpdom))
myd %>% ggplot (aes(log10cpdomexp, log10cpdom) )+ geom_point() + geom_abline(intercept=yint, slope=cplambda95, linetype=3, color='grey')+ geom_abline(intercept=0, slope=1, color='grey') + xlim(0, maxp)+ylim(0, maxp) + theme_minimal() +labs(x='expected -log10(p-value)' , y='observed -log10(p-value)' , title= paste('Dominant test - cplambda95 =' ,round(cplambda95, 3)  , sep= ' ' )  )

#recessive
#calulate expected pval
myd<-myd[order(myd$log10prec),]
myd$log10prec_exp<-sort(-log10(runif(nrow(myd)))) ## 'nrows' random values from uniform distribution, sorted
#calculate lambda95
rlp95obs <- quantile( subset(myd, log10prec!=0)$log10prec, 0.95) #p95
rlp95exp <- quantile( subset(myd, log10prec!=0)$log10prec_exp, 0.95) #u95
rlp0exp <- maxlogp1 <- max(subset(myd, prec==1)$log10prec_exp) #u0
rlp0obs <-0 #p0
rlambda95 <-(rlp95obs[[1]]-rlp0obs) / (rlp95exp[[1]]-rlp0exp) #slope
ryint<-rlp95obs[[1]]-rlambda95*rlp95exp[[1]]
rmaxp<-ceiling(max(myd$log10prec))
myd %>% ggplot (aes(log10prec_exp, log10prec) )+ geom_point() + geom_abline(intercept=ryint, slope=rlambda95, linetype=3, color='grey')+ geom_abline(intercept=0, slope=1, color='grey') + xlim(0, rmaxp)+ylim(0, rmaxp) + theme_minimal() +labs(x='expected -log10(p-value)' , y='observed -log10(p-value)' , title= paste('Recessive test - rlambda95 =' ,round(rlambda95, 3)  , sep= ' ' )  )
ggsave('dominant.png')



samplingXsquare<-rchisq(10520,df=1)
pchisq(samplingXsquare, df=1)

expPX<-pchisq(samplingXsquare, df=1)
myd<-myd[order(myd$log10cpdom),]
myd$log10cpdomexp<-sort(-log10(expPX))
