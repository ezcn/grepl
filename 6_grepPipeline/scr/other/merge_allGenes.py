############open all gene files

file1 = pd.read_csv('/lustrehome/silvia/genes/embryo.tsv', sep="\t")
file2 = pd.read_csv('/lustrehome/silvia/genes/DDD.tsv', sep="\t")
file3 = pd.read_csv('/lustrehome/silvia/genes/lethal.tsv', sep="\t")
file4 = pd.read_csv('/lustrehome/silvia/genes/essential.tsv', sep="\t")
file5 = pd.read_csv('/lustrehome/silvia/genes/miscarriage.tsv', sep="\t")

############ rename columns

file1.columns = [ 'ensID', 'gene_type']                                                                                                                    
file2.columns = [ 'ensID', 'gene_type']
file3.columns = [ 'ensID', 'gene_type']
file4.columns = [ 'ensID', 'gene_type']
file5.columns = [ 'ensID', 'gene_type']

############ insert new column 'value'

file1.insert(1, "value", "1")
file2.insert(1, "value", "1")
file3.insert(1, "value", "1")
file4.insert(1, "value", "1")
file5.insert(1, "value", "1")

############ reshape

df1=file1.pivot(index='ensID', columns='gene_type', values='value')
df2=file2.pivot(index='ensID', columns='gene_type', values='value')
df3=file3.pivot(index='ensID', columns='gene_type', values='value')
df4=file4.pivot(index='ensID', columns='gene_type', values='value')

############## merge all lists

list=df1.reset_index().merge(df2, on="ensID", how="outer").fillna(0)
list2=list.merge(df3, on="ensID", how="outer").fillna(0)
list3=list2.merge(df4, on="ensID", how="outer").fillna(0)
list_final=list3.merge(df5, on="ensID", how="outer").fillna(0)

############ write file tsv

list_final.to_csv("/lustrehome/silvia/genes/all_geneList.tsv", sep="\t", index=False)
