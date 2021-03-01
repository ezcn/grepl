import sqlite3
import pandas as pd
import re, sys, argparse

def main():
	parser = argparse.ArgumentParser()
	###~~~~ new sql db name
	parser.add_argument("-db", help="path to input vep table file ",type=str,required=True)
	###~~~ input files
	parser.add_argument("-i", help="path to input vep table file ",type=str,required=True)
	###~~~~ databases
	parser.add_argument("-g", help="path to gene file list ",required=True)
	parser.add_argument("-pli", help="path to table of pLI score  ",type=str, required= True)
	#parser.add_argument('-fathmmcoding', help="path to fathmm file",type=str, required= True)
	#parser.add_argument('-fathmmnc', help="path to fathmm coding file",type=str, required= True)
	###~~~ output file
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()


	conn = sqlite3.connect(args.db)  ### create new sql db
	c = conn.cursor()
	c.execute("DROP TABLE IF EXISTS myTable;")
	c.execute("DROP TABLE IF EXISTS pLItable;")
	c.execute("DROP TABLE IF EXISTS geneList;")
	c.execute("DROP TABLE IF EXISTS fatmTab;")
	c.execute("DROP TABLE IF EXISTS fatmNCtab;")
	###### open vep table with pandas and create table inside grep.db
	df = pd.read_table(args.i, sep="\t", index_col= "Uploaded_variation")
	df.columns = df.columns.str.strip()
	df.to_sql("myTable", conn)
	###### open pLI table and create new table inside grep.db
	df=pd.read_table(args.pli, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("pLItable", conn)
	###### open geneList table and create new table inside grep.db
	df=pd.read_table(args.g, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("geneList", conn)
	##### open fathmm coding table and create new table inside grep.db
	'''
	df=pd.read_table(args.fathmmcoding, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("fatmTab", conn)
	###### open fathmm non coding
	df=pd.read_table(args.fathmmnc, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("fatmNCtab", conn)
	'''
	
	###### join myTable and pLItab on common column creating new table

	c.execute('CREATE TABLE firstjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real);')

	c.execute('INSERT INTO firstjoin SELECT myTable.*, pLItable.pLI FROM myTable LEFT JOIN pLItable ON myTable.Feature = pLItable.transcript;')
	conn.commit()


	##### join genelist

	c.execute('CREATE TABLE genesjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real);')

	c.execute('INSERT INTO genesjoin SELECT firstjoin.*, geneList.EmbryoDev, geneList.DDD, geneList.Lethal, geneList.Essential, geneList.Misc FROM firstjoin LEFT JOIN geneList ON firstjoin.Gene = geneList.ensID;')
	conn.commit()

	###### add index_x column

	c.execute('CREATE TABLE indexTable (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text);')

	c.execute("INSERT INTO indexTable SELECT *, Location ||\':/\'|| Allele FROM genesjoin;")

	conn.commit()


	####### join fathmm score
	#c.execute('CREATE TABLE codjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real);')

	#c.execute('INSERT INTO codjoin SELECT indexTable.*, fatmTab.FathmmCoding FROM indexTable LEFT JOIN fatmTab ON indexTable.index_x = fatmTab.key;')
	#conn.commit()
	###### join fathm non coding

	#c.execute('CREATE TABLE noncodjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AF real,AFR_AF real,AMR_AF real,ASN_AF real,EUR_AF real,EAS_AF real,SAS_AF real,AA_AF real,EA_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,MAX_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, index_x text, FathmmCod real, FathmmNonCod real);')

	#c.execute('INSERT INTO noncodjoin SELECT codjoin.*, fatmNCtab.FathmmNonCod FROM codjoin LEFT JOIN fatmNCtab ON codjoin.index_x = fatmNCtab.key;')
	#conn.commit()

	
	query = "SELECT * FROM indexTable;" 
	df_final = pd.read_sql_query(query, conn).fillna(0)
	df_final.to_csv(args.o, sep="\t", index=False)

if __name__ == "__main__":
	main()

