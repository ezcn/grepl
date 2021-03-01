import sqlite3
import pandas as pd
import re, sys, argparse, random, glob, subprocess

def main():
	parser = argparse.ArgumentParser()
	###~~~~ new sql db name
	parser.add_argument("-dbC", help="path to control db ",type=str,required=True)
	parser.add_argument("-dbS", help="path to control db ",type=str,required=True)
	###~~~ input files
	#parser.add_argument("-scsq", help="path to input grep csq file ",type=str,required=True)
	#parser.add_argument("-ccsq", help="path to input hgdp csq file ",type=str,required=True)
	parser.add_argument("-cl", help="list of control samples id ",type=str,required=True)
	parser.add_argument("-sl", help="list of sampes id ",type=str,required=True)
	###~~~ filtering arguments
	parser.add_argument("-ff", help="threshold for first rare frequency definition (01) ", type=float,required=True)
	parser.add_argument("-f", help="threshold for second rare frequency definition (05)", type=float,required=True)
	parser.add_argument("-type", help=" feature_type(genic , intergenic, regularoty) ", type=str,required=True)
	parser.add_argument("-r", help=" rare variant not equal to (05)", type=str,required=True)
	#parser.add_argument("-rr", help=" rare variant for second threshold not equal to (05)", type=str,required=True)
	parser.add_argument("-pli", help="threshold for  pLI score  ",type=float, required= True)
	parser.add_argument("-cadd", help=" treshold for CADD score  ",type=float , required= True)
	parser.add_argument("-g", help=" number of gene lists  ",type=float , required= True)
	###~~~ HGDP arguments
	parser.add_argument("-n", help="number of individual to sample ",type=int,required=True)
	parser.add_argument("-i", help="number of iterations ",type=int,required=True)
	parser.add_argument("-pathTodirCtrl", help="path to file with control allele counts directory  ",type=str, required= True)
	parser.add_argument("-ctgn", help="path to control genes file",type=str,required=True)
	parser.add_argument("-gtd", help="path to genes to discard file",type=str,required=True)

	###~~~ Grep arguments
	parser.add_argument("-pathTodir", help="path to file with sample allele counts directory  ",type=str, required= True)
	parser.add_argument("-chrom", help="chromosome name  ",type=str, required= True)
	parser.add_argument("-ac", help=" allele count >= of   ",type=int , required= True)
	#parser.add_argument("-ctrlgen", help="path to hgdp genes to discard file ",type=str,required=True)
	parser.add_argument("-gt", help="threshold for excluding genes ",type=float,required=True)
	#parser.add_argument("-v", help=" path to file for variants per gene count ",type=str , required= True)
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()

	##### open hgdp db
	conn = sqlite3.connect(args.dbC)  ### create new sql db
	c = conn.cursor()
	##### remove tables pli, fathmm and gene list
	c.execute("DROP TABLE IF EXISTS pLItable;")
	c.execute("DROP TABLE IF EXISTS geneList;")
	c.execute("DROP TABLE IF EXISTS fatmTab;")
	c.execute("DROP TABLE IF EXISTS fatmNCtab;")
	########## Replace - with 0 in frequencies columns
	c.execute("UPDATE indexTable SET AFR_AF = REPLACE (AFR_AF, '-', 0) , AMR_AF = REPLACE (AMR_AF, '-', 0) , ASN_AF = REPLACE (ASN_AF, '-', 0) , EUR_AF = REPLACE (EUR_AF, '-', 0) , EAS_AF = REPLACE (EAS_AF, '-', 0) , SAS_AF = REPLACE (SAS_AF, '-', 0) , AA_AF = REPLACE (AA_AF, '-', 0) , EA_AF = REPLACE (EA_AF, '-', 0) , gnomAD_AF = REPLACE (gnomAD_AF, '-', 0), gnomAD_AFR_AF = REPLACE (gnomAD_AFR_AF, '-', 0) ,  gnomAD_AMR_AF = REPLACE (gnomAD_AMR_AF, '-', 0) , gnomAD_ASJ_AF = REPLACE (gnomAD_ASJ_AF, '-', 0), gnomAD_EAS_AF = REPLACE (gnomAD_EAS_AF, '-', 0), gnomAD_FIN_AF = REPLACE (gnomAD_FIN_AF, '-', 0), gnomAD_NFE_AF = REPLACE (gnomAD_NFE_AF, '-', 0), gnomAD_OTH_AF = REPLACE (gnomAD_OTH_AF, '-', 0), gnomAD_SAS_AF = REPLACE (gnomAD_SAS_AF, '-', 0),  CADD_RAW = REPLACE (CADD_RAW, '-', 0);")
	####### add info to hgdp table
	thr = args.ff
	c.execute("DROP TABLE IF EXISTS freqTable1Ctr;")
	c.execute('CREATE TABLE freqTable1Ctr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real);')

	c.execute("INSERT INTO freqTable1Ctr SELECT *, CASE WHEN AFR_AF =0 AND AMR_AF =0 AND ASN_AF =0 AND EUR_AF =0 AND EAS_AF =0 AND SAS_AF =0 AND AA_AF =0 AND EA_AF =0 AND gnomAD_AFR_AF =0 AND gnomAD_AMR_AF =0 AND gnomAD_ASJ_AF =0 AND gnomAD_EAS_AF =0 AND gnomAD_FIN_AF =0 AND gnomAD_NFE_AF =0 AND gnomAD_OTH_AF =0 AND gnomAD_SAS_AF =0 THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM indexTable;" ,(thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,))

	conn.commit()

	thr2 = args.f
	c.execute("DROP TABLE IF EXISTS freqTableCtr;")
	c.execute('CREATE TABLE freqTableCtr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real, rare05 real);')

	c.execute("INSERT INTO freqTableCtr SELECT *, CASE WHEN AFR_AF =0 AND AMR_AF =0 AND ASN_AF =0 AND EUR_AF =0 AND EAS_AF =0 AND SAS_AF =0 AND AA_AF =0 AND EA_AF =0 AND gnomAD_AFR_AF =0 AND gnomAD_AMR_AF =0 AND gnomAD_ASJ_AF =0 AND gnomAD_EAS_AF =0 AND gnomAD_FIN_AF =0 AND gnomAD_NFE_AF =0 AND gnomAD_OTH_AF =0 AND gnomAD_SAS_AF =0 THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM freqTable1Ctr;" ,(thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,))

	conn.commit()
	###### create sumGene column in control samples df
	c.execute("DROP TABLE IF EXISTS sumgeneTabCtr;")
	c.execute('CREATE TABLE sumgeneTabCtr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real, rare05 real, sumGene real);')

	c.execute("INSERT INTO sumgeneTabCtr SELECT Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,MAX_AF,CADD_RAW, CADD_PHRED, pLIscore, EmbryoDev, DDD, Lethal, Essential, Misc, index_x text, rare01,rare05, (EmbryoDev + DDD + Lethal + Essential + Misc ) FROM freqTableCtr;")

	conn.commit()

	###### create CADD percentile column in control samples df
	c.execute("DROP TABLE IF EXISTS finalTableCtr;")
	c.execute("CREATE TABLE finalTableCtr (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real, rare05 real, sumGene real, caddPercent real);")


	c.execute("INSERT INTO finalTableCtr SELECT * , ROUND(PERCENT_RANK() OVER (ORDER BY CADD_RAW),3) FROM sumgeneTabCtr WHERE CADD_RAW != '-';")

	conn.commit()

	######## filter control samples 
	typ = args.type ; rareThresh = args.r  ; pliscore = args.pli ; caddscore = args.cadd ; numgene = args.g

	#c.execute("CREATE TABLE filtro AS SELECT * FROM finalTable WHERE IMPACT != 'MODIFIER' AND feature_type = ? AND rare != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?);" , (typ,rareThresh, pliscore, caddscore, numgene,))
	query = "SELECT * FROM finalTableCtr WHERE IMPACT != 'MODIFIER' AND Feature_type = ? AND rare05 != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?); "
	df_finalCtr = pd.read_sql_query(query,conn, params = (typ,rareThresh, pliscore, caddscore, numgene))
	#query = "SELECT * FROM filtro;"
	#df_finalHg.to_csv()
	#df_finalHg = pd.read_sql_query(query,conn)

	##~~ repeat args.i times on args.n samples from controls
	listControl = [line.rstrip('\n') for line in open(args.cl)]
	cycle=0
	while cycle < args.i:
		cycle+=1
		##~~ choose a random sample form controls of size args.n
		randomSample=random.sample(listControl, args.n)
		##~~ integrate with samples ID and csqAlele count
		control_filtered_allsamples=pd.DataFrame()
		for ss in randomSample:
			tmpdf=df_finalCtr
			tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %( args.pathTodirCtrl, ss, args.chrom) )
			if "chr" in tmpSamp.loc[0].key:
				tmpSamp.key = tmpSamp.key.str.lstrip("chr")
			tmpSamp = tmpSamp[ (tmpSamp['ALTcount']>=args.ac )]
			tmpSamp["sample"] = ss
			if "chr" in tmpdf.loc[0].index_x:
				tmpdf.index_x = tmpdf.index_x.str.lstrip("chr")
			df=tmpdf.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
			control_filtered_allsamples=pd.concat([control_filtered_allsamples, df])
		tmpg=control_filtered_allsamples[[ 'SYMBOL', 'ID']].drop_duplicates().groupby('SYMBOL').count().transform(lambda x: x /float(args.n) )
		if cycle==1:  genesPerSample=tmpg  # create genesPerSample at first cycle
		else: genesPerSample=genesPerSample.join(tmpg, on='SYMBOL', how='outer', lsuffix='_genesPerSample', rsuffix='_tmpg').fillna(0)
		
	##~~ make GrandMean over args.i and discard genes that on average shows up in args.gt individuals over args.i iterations	
	genesPerSample['GrandMean']=genesPerSample.sum(axis=1 )/float(args.i)
	genesToDiscardControl=genesPerSample[genesPerSample['GrandMean']> float(args.gt)]
	genesPerSample.to_csv(args.ctgn, sep='\t', index=False)
	genesToDiscardControl.to_csv(args.gtd, sep='\t', index=False)

	conn.close()
	#### open grep db
	conn = sqlite3.connect(args.dbS) 
	c = conn.cursor()
	#### load table genePerSample in grep db
	genesPerSample = pd.read_table(args.ctgn, sep="\t")
	genesPerSample.columns = genesPerSample.columns.str.strip()
	genesPerSample.to_sql("genesMean", conn, if_exists="replace")

	c.execute("UPDATE indexTable SET AFR_AF = REPLACE (AFR_AF, '-', 0) , AMR_AF = REPLACE (AMR_AF, '-', 0) , ASN_AF = REPLACE (ASN_AF, '-', 0) , EUR_AF = REPLACE (EUR_AF, '-', 0) , EAS_AF = REPLACE (EAS_AF, '-', 0) , SAS_AF = REPLACE (SAS_AF, '-', 0) , AA_AF = REPLACE (AA_AF, '-', 0) , EA_AF = REPLACE (EA_AF, '-', 0) , gnomAD_AF = REPLACE (gnomAD_AF, '-', 0), gnomAD_AFR_AF = REPLACE (gnomAD_AFR_AF, '-', 0) ,  gnomAD_AMR_AF = REPLACE (gnomAD_AMR_AF, '-', 0) , gnomAD_ASJ_AF = REPLACE (gnomAD_ASJ_AF, '-', 0), gnomAD_EAS_AF = REPLACE (gnomAD_EAS_AF, '-', 0), gnomAD_FIN_AF = REPLACE (gnomAD_FIN_AF, '-', 0), gnomAD_NFE_AF = REPLACE (gnomAD_NFE_AF, '-', 0), gnomAD_OTH_AF = REPLACE (gnomAD_OTH_AF, '-', 0), gnomAD_SAS_AF = REPLACE (gnomAD_SAS_AF, '-', 0),  CADD_RAW = REPLACE (CADD_RAW, '-', 0);")


	###### retrieve frequencies Grep samples
	thr = args.ff
	c.execute("DROP TABLE IF EXISTS freqTable1;")
	c.execute('CREATE TABLE freqTable1 (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real);')

	c.execute("INSERT INTO freqTable1 SELECT *, CASE WHEN AFR_AF =0 AND AMR_AF =0 AND ASN_AF =0 AND EUR_AF =0 AND EAS_AF =0 AND SAS_AF =0 AND AA_AF =0 AND EA_AF =0 AND gnomAD_AFR_AF =0 AND gnomAD_AMR_AF =0 AND gnomAD_ASJ_AF =0 AND gnomAD_EAS_AF =0 AND gnomAD_FIN_AF =0 AND gnomAD_NFE_AF =0 AND gnomAD_OTH_AF =0 AND gnomAD_SAS_AF =0 THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM indexTable;" ,(thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,))

	conn.commit()

	thr2 = args.f
	c.execute("DROP TABLE IF EXISTS freqTable;")
	c.execute('CREATE TABLE freqTable (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real, rare05 real);')

	c.execute("INSERT INTO freqTable SELECT *, CASE WHEN AFR_AF =0 AND AMR_AF =0 AND ASN_AF =0 AND EUR_AF =0 AND EAS_AF =0 AND SAS_AF =0 AND AA_AF =0 AND EA_AF =0 AND gnomAD_AFR_AF =0 AND gnomAD_AMR_AF =0 AND gnomAD_ASJ_AF =0 AND gnomAD_EAS_AF =0 AND gnomAD_FIN_AF =0 AND gnomAD_NFE_AF =0 AND gnomAD_OTH_AF =0 AND gnomAD_SAS_AF =0 THEN 'NOB' WHEN AFR_AF <=? AND AMR_AF <=? AND ASN_AF <=? AND EUR_AF <=? AND EAS_AF <=? AND SAS_AF <=? AND AA_AF <=? AND EA_AF <=? AND gnomAD_AFR_AF <=? AND gnomAD_AMR_AF <=? AND gnomAD_ASJ_AF <=? AND gnomAD_EAS_AF <=? AND gnomAD_FIN_AF <=? AND gnomAD_NFE_AF <=? AND gnomAD_OTH_AF <=? AND gnomAD_SAS_AF <=? THEN 'true' ELSE 'false' END FROM freqTable1;" ,(thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,))

	conn.commit()
	###### create sumGene column grep samples
	c.execute("DROP TABLE IF EXISTS sumgeneTab;")
	c.execute('CREATE TABLE sumgeneTab (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real, rare05 real, sumGene real);')

	c.execute("INSERT INTO sumgeneTab SELECT Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,MAX_AF,CADD_RAW, CADD_PHRED, pLIscore, EmbryoDev, DDD, Lethal, Essential, Misc, index_x text, rare01,rare05, (EmbryoDev + DDD + Lethal + Essential + Misc) FROM freqTable;")

	conn.commit()

	###### create CADD percentile column grep samples
	c.execute("DROP TABLE IF EXISTS finalTable;")
	c.execute("CREATE TABLE finalTable (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, rare01 real, rare05 real, sumGene real, caddPercent real);")


	c.execute("INSERT INTO finalTable SELECT * , ROUND(PERCENT_RANK() OVER (ORDER BY CADD_RAW),3) FROM sumgeneTab WHERE CADD_RAW != '-';")

	conn.commit()

	########## FILTERING grep samples

	typ = args.type ; rareThresh = args.r ; pliscore = args.pli ; caddscore = args.cadd ; numgene = args.g

	#c.execute("CREATE TABLE filtro AS SELECT * FROM finalTable WHERE IMPACT != 'MODIFIER' AND feature_type = ? AND rare != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?);" , (typ,rareThresh, pliscore, caddscore, numgene,))
	query = "SELECT * FROM finalTable WHERE IMPACT != 'MODIFIER' AND Feature_type = ?  AND rare05 != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGene >= ?); "
	df_final = pd.read_sql_query(query,conn, params = (typ,rareThresh, pliscore, caddscore, numgene))
	
	listSamples = [line.rstrip('\n') for line in open(args.sl)]
	ssall=pd.DataFrame()
	for ss in listSamples:
		tmpgrep=df_final
		tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %( args.pathTodir, ss, args.chrom) )
		if "chr" in tmpSamp.loc[0].key:
			tmpSamp.key = tmpSamp.key.str.lstrip("chr")
		if "chr" in tmpgrep.loc[0].index_x:
			tmpgrep.index_x = tmpgrep.index_x.str.lstrip("chr")
		tmpSamp = tmpSamp[ (tmpSamp['ALTcount']>= args.ac )]
		df=tmpgrep.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
		ssall=pd.concat([ssall, df])
		#df_final.to_csv("filtropdaslq" , sep="\t", index=False)
		

	ssall.to_sql("grepFilter", conn,if_exists="replace")
	c.execute("DROP TABLE IF EXISTS noCtrlGenes;")
	c.execute("CREATE TABLE noCtrlGenes AS SELECT grepFilter.*, CASE WHEN genesMean.GrandMean IS NULL THEN 0 ELSE genesMean.GrandMean END as GrandMean FROM grepFilter LEFT JOIN genesMean ON grepFilter.SYMBOL = genesMean.SYMBOL;")
	#hgdpMean = args.gt
	#c.execute("CREATE TABLE noCtrlGenes AS SELECT grepFilter.* FROM grepFilter LEFT JOIN genesMean ON grepFilter.SYMBOL = genesMean.SYMBOL WHERE GrandMean <= ?; " , (hgdpMean,))
	#variants = args.maxv
	#c.execute("CREATE TABLE variantsPerGene AS SELECT index_x, SYMBOL ,COUNT (DISTINCT index_x) FROM noCtrlGenes GROUP BY SYMBOL;")
	#HAVING COUNT (DISTINCT index_x) >=?;" ,(variants,))
	#query1 = "SELECT * FROM variantsPerGene;"
	#variants = pd.read_sql_query(query1,conn)
	#variants.to_csv(args.v, sep = "\t", index = False)
	#c.execute("DELETE FROM noCtrlGenes WHERE EXISTS (SELECT * FROM variantsPerGene WHERE variantsPerGene.SYMBOL = noCtrlGenes.SYMBOL);")
	conn.commit()
	query = "SELECT * FROM noCtrlGenes;"
	df_end = pd.read_sql_query(query,conn)
	df_end.to_csv(args.o, sep = "\t", na_rep= "NA", index = False)


if __name__ == "__main__":
	main()
