import re , numpy 
import pandas as pd 

sites = [line.rstrip('\n').split('-')[0].lstrip('chr')  for line in open('/home/enza/ezcngit/grep/grepPipeline/testdata/test.sites')]
#print (sites)

#samples  = [line.rstrip('\n') for line in open('/home/enza/ezcngit/grep/grepPipeline/testdata/test.samples')]
#count=0
#for sa in samples:
#    count+=1
#   fname='/home/enza/ezcngit/grep/grepPipeline/testdata/con/%s.chr22.tsv' % (sa)
    #df=pd.DataFrame(pd.read_csv(fname,sep="\t"))#, index_col='index_x'))
    #df[['chr','pos', 'all']] = df.index_x.str.split(expand=True)

    #df.assign(site = df['index_x'].str.split(':')[0]+':' +a.split(':')[1])           
    #if count==1 : alldf=df
    #else: alldf.merge(df)                   
#with open('/home/enza/ezcngit/grep/grepPipeline/testdata/con/all.hgdp.csq') as f:
with open('/home/enza/ezcngit/grep/grepPipeline/testdata/sam/all.sam.csq') as f:
    for line in f:
        key=line.split()[0].split(':')[0]+':' +line.split(':')[1]
        #print(key)  
        if key in sites: 
            print(line.rstrip())
