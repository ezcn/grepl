import sqlite3
import pandas as pd
import re, sys, argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-db", help="path to csq db  ",type=str,required=True)
	parser.add_argument("-cg", help="path to control genes file ",type=str,required=True)
	parser.add_argument("-sl", help="list of sampes id ",type=str,required=True)
	parser.add_argument("-pathTodir", help="path to file with sample allele counts directory  ",type=str, required= True)
	parser.add_argument("-chrom", help="chromosome name  ",type=str, required= True)
	parser.add_argument("-ac", help=" allele count >= of   ",type=int , required= True)
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()

	#### access to db
	conn = sqlite3.connect(args.db)
	c = conn.cursor()
	c.execute("DROP TABLE IF EXISTS CandidateGenes;")
	### open file with candidate genes
	df=pd.read_table(args.cg)
	df.columns = df.columns.str.strip()
	df.to_sql("CandidateGenes", conn)

	###### find candidate genes in csq files
	query = "SELECT * FROM noncodjoin JOIN CandidateGenes ON noncodjoin.Gene = CandidateGenes.Gene;"
	dfCtr = pd.read_sql_query(query, conn)

	##### add sample and AltCounts

	def countedSC(dfCtr, listSamples, pathTodir, chrom, output):
		
		tmpPda = dfCtr
		if "chr" in tmpPda.index_x.unique()[0]:
			tmpPda.index_x = tmpPda.index_x.str.lstrip("chr")
		
		allCandidateGenes=pd.DataFrame()
		
		for sample in listSamples:
			tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %( pathTodir, sample, chrom) )
			if "chr" in tmpSamp.key.unique()[0]:
				tmpSamp.key = tmpSamp.key.str.lstrip("chr")
			df=tmpPda.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
			allCandidateGenes=pd.concat([allCandidateGenes, df])
		
		allCandidateGenes.to_csv(output, sep = "\t", na_rep= "NA", index = False)

	
	listSamples = [line.rstrip('\n') for line in open(args.sl)]
	
	try:
		countedSC(dfCtr, listSamples, args.pathTodir, args.chrom, args.o)
	except Exception as e:
		print("Candidate genes absents, skip...")
	

if __name__ == "__main__":
	main()