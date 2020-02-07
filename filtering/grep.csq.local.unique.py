#!/usr/bin/python3
import pandas as pd
import numpy as np
import re, sys, argparse, gzip, json
import dask.dataframe as dd
from ast import literal_eval
import ast
import decimal

###### Function START ###################################
###################################################
# using the NodeTransformer, you can also modify the nodes in the tree,
# however in this example NodeVisitor could do as we are raising exceptions
# only.
class Transformer(ast.NodeTransformer):
    ALLOWED_NAMES = set(['Decimal', 'None', 'False', 'True'])
    ALLOWED_NODE_TYPES = set([
        'Expression', # a top node for an expression
        'Tuple',      # makes a tuple
        'Call',       # a function call (hint, Decimal())
        'Name',       # an identifier...
        'Load',       # loads a value of a variable with given identifier
        'Str',        # a string literal

        'Num',        # allow numbers too
        'List',       # and list literals
        'Dict',       # and dicts...
    ])

    def visit_Name(self, node):
        if not node.id in self.ALLOWED_NAMES:
            raise RuntimeError("Name access to %s is not allowed" % node.id)

        # traverse to child nodes
        return self.generic_visit(node)

    def generic_visit(self, node):
        nodetype = type(node).__name__
        if nodetype not in self.ALLOWED_NODE_TYPES:
            raise RuntimeError("Invalid expression: %s not allowed" % nodetype)

        return ast.NodeTransformer.generic_visit(self, node)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getInfoFromVepLocally (jsonWithVEPannotations):
        vepInfo={}; vepInfoCommon={}    
        #~~ filename is the json output of VEP runned locally 
        #~~ filename is parsed by line and by alternate allele  
        with open(jsonWithVEPannotations, 'r') as f:
                for line in f:
                        info=json.loads(line) #~~ make a dictionary out of a string 
                        #~~ derive the key chr:pos:/alt from the "input" element and process each alternate allele 
                        locusdata=info['input'].split(); altAlleles=locusdata[4].split(","); mychr=locusdata[0]; mypos=locusdata[1]
                        for altAl in altAlleles:
                                mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
                                #print(mykey) 
                                most=info["most_severe_consequence"]
                                vepInfoCommon[mykey]={}
                                vepInfoCommon[mykey]['most_severe_consequence']=most 
                                #~~  check if the most sequence is in a transcript
                                if 'transcript_consequences' in info:
                                        for tc in info['transcript_consequences']:
                                                if most in  tc['consequence_terms'] :
                                                        tcTranscript=tc['transcript_id']
                                                        vepInfo[(mykey, tcTranscript)]={}
                                                        tcAllele=tc['variant_allele']
                                                        if tcAllele ==altAl :
                                                                vepInfo[(mykey, tcTranscript)]['csqAllele']=tcAllele
                                                                vepInfo[(mykey, tcTranscript)]['gene_id']=tc['gene_id']
                                                                vepInfo[(mykey, tcTranscript)]['gene_symbol']= tc['gene_symbol']
                                                                vepInfo[(mykey, tcTranscript)]['impact']=tc['impact']
                                                                vepInfo[(mykey, tcTranscript)]['key']=mykey
                                                                vepInfo[(mykey, tcTranscript)]['element_id']=tc['transcript_id']
                                #~~ check if the most severe consequence is in a regulatory feature  
                                elif 'regulatory_feature_consequences' in info:
                                        for rf  in info['regulatory_feature_consequences']:
                                                if most in  rf['consequence_terms']:
                                                        rfRegulatory=rf['regulatory_feature_id']
                                                        vepInfo[(mykey, rfRegulatory)]={}
                                                        rfAllele=rf['variant_allele']
                                                        if rfAllele ==altAl :
                                                                vepInfo[(mykey, rfRegulatory)]['csqAllele']=rfAllele
                                                                vepInfo[(mykey, rfRegulatory)]['impact']=rf['impact']
                                                                vepInfo[(mykey, rfRegulatory)]['key']=mykey
                                                                vepInfo[(mykey, rfRegulatory)]['element_id']=rf['regulatory_feature_id']
                                #~~ check if the most severe consequence is in an intergenic region 
                                elif 'intergenic_consequences' in info:
                                        for ic in info['intergenic_consequences']:
                                                if most in  ic['consequence_terms']:
                                                        vepInfo[(mykey, 'intergenic')]={}
                                                        icAllele=ic['variant_allele']
                                                        if icAllele ==altAl:    
                                                                vepInfo[(mykey, 'intergenic') ]['csqAllele']=icAllele
                                                                vepInfo[(mykey, 'intergenic') ]['key']=mykey
                                                                vepInfo[(mykey, 'intergenic') ]['element_id']='intergenic'
                                #~~ retrive info on rsID, starting position, frequencies for the csqAllele found in the previous code  
                                if "colocated_variants" in info:
                                        infoCV = info["colocated_variants"][0]
                                        #~~ adding info for starting position, this is needed for merging with CADD score.
                                        vepInfoCommon[mykey]["start"] = infoCV["start"]

                                        #~~ check if rsid is present 
                                        if "id" in infoCV: vepInfoCommon[mykey]["id"] = infoCV["id"]
                                        
                                        #~~ retrieve allelle frequencies 
                                        if 'frequencies' in infoCV:
                                                vepInfoCommon[mykey]["frequencies"]= infoCV["frequencies"]


    #print (json.dumps(info, indent=4 ) )  
    #print (json.dumps(vepInfo, indent=4) ) 
        return vepInfo, vepInfoCommon





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def VepSOTermInfo (vepinfofile): 
    """read external file with info on VEP consequences  """
    lSOTerm=[]  ### list of SOTerm 
    
    countlinesCsq= True
    for csqLine in open(vepinfofile, 'r'):
        if countlinesCsq:
            csqTitle=csqLine.rstrip().split('\t')
            countlinesCsq=False
        else:
            myRowList=csqLine.rstrip().split('\t')
            lSOTerm.append(myRowList[0])

    return  lSOTerm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def VepRankingInfo (vepinfofile): 
    """read external file with info on VEP consequences  """
    dRank={"HIGH":4, "LOW": 2, "MODERATE":3, "MODIFIER":1}
    dSOTermRank={}
    lSOTerm=[]  ### list of SOTerm ordered by severity

    countlinesCsq= True
    for csqLine in open(vepinfofile, 'r'):
        if countlinesCsq:
            csqTitle=csqLine.rstrip().split('\t')
            countlinesCsq=False
        else:
            myRowList=csqLine.rstrip().split('\t')
            dCsq= dict(zip(csqTitle, myRowList ))
            dSOTermRank[dCsq['SO term']]=dRank[dCsq['IMPACT']]
            lSOTerm.append(myRowList[0])

    #print (lSOTerm)
    lScores=list(reversed(range(len(lSOTerm)))) 
    #print (lScores) 
    dSOTermFineRank=dict(zip(lSOTerm, map(int, lScores) ))
    #print (dSOTermFineRank)
    return dSOTermFineRank
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getInfoFromVep (Position):
    """Retrieves information from Variant Effect Predictor API
    dependencies: requests, json 
    Position =  1:333333:/T (T is the alternate allele)   """
    freq_dict={}
    server="https://rest.ensembl.org"
    ext = "/vep/human/region/"+ Position +"?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        return Position
        #r.raise_for_status()
        #sys.exit()
    decoded= r.json()
    info = decoded[0]
    if "colocated_variants" in info:
        id_search = info["colocated_variants"][0]
        if "id" in id_search:
            freq_dict["id"] = id_search["id"]
        # adding info for starting position, this is needed for merging with CADD score.
        #if "start" in id_search: 
        #    freq_dict["starting_position"] = id_search["start"]   
        if 'frequencies' in id_search:
            to_parse = list(id_search["frequencies"].values())[0]
            for var in to_parse.keys():
                freq_dict[var]=to_parse[var]

    if "most_severe_consequence" in info:
        freq_dict["most_severe_consequence"]=info["most_severe_consequence"]
        most=freq_dict["most_severe_consequence"]
        if 'transcript_consequences' in info: 
            for i in info['transcript_consequences']:
                if most in  i['consequence_terms'] :
                    csqAllele=i['variant_allele']
                    freq_dict['csqAllele']=csqAllele
                    freq_dict['gene_id']=i['gene_id']
                    freq_dict['gene_symbol']=i['gene_symbol']
                else:
                    if 'regulatory_feature_consequences' in info: 
                        for r  in info['regulatory_feature_consequences']:
                            if most in  r['consequence_terms']: 
                                csqAllele=r['variant_allele']
                                freq_dict['csqAllele']=csqAllel
    return freq_dict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def csqAlleleCount (csqAllele, altAllele, altAlleleCount, ploidy): 
    #csqAllele= cons allele  from vep  
    #altAlleleCount= integer, alternate allele count 
    # ploidy = number of allels at one locus 
    
    if csqAllele == altAllele: mycsqAlleleCount = altAlleleCount      ## Consequence allele is the alternate allele 
    else: mycsqAlleleCount = ploidy-altAlleleCount    ## Consequence allele is the reference allele;
    
    return mycsqAlleleCount
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def csqAlleleFeaturesLike (csqAllele, altAllele, altAlleleCount, GL ): 
    #csqAllele= cons allele  from vep  
    #altAlleleCount= integer, alternate allele count 
    #GL, string, comma separated vector of genotpe likelihoods
    #"""

    if csqAllele == altAllele: 
        mycsqAlleleCount = altAlleleCount      ## csqAllele counts  
    else: 
        mycsqAlleleCount = 2-altAlleleCount 
    Likl=GL.split(",")[altAlleleCount] ## csqAllele Likelihood 
    return [ csqAllele, mycsqAlleleCount, Likl]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def csqAlleleFeaturesMulti(genotype, csqAllele, refAllele, altAlleles, altAlleleCount, GL): 
    allAlleles=[refAllele]+ altAlleles.split(',')
    templist=altAlleleCount.split(',')
    allCounts=[2-sum( map(int, templist) )]+templist
    dCount=dict(zip(allAlleles, allCounts)) 
    if csqAllele in dCount:  
        csqCount=int(dCount[csqAllele]) 

        """ there is no itertool or numpy therfore I HAVE TO write the following...
        it MUST be rewritten to take the upperdiagonal 
        use numpy
        """

        tempGenolist=[str(a)+str(b) for a in range (len(allAlleles)) for b in range (len(allAlleles)) ]
        """This is equivalent to: 

        tempGenolist=[]
        for a in range (len(allAlleles)):        
            for b in range (len(allAlleles)): 
                tempGenolist.append(str(a) +str(b)) 
        """
        if len(allAlleles)==2: genotypesIndexToRetain=[0,1,3]
        elif len(allAlleles)==3: genotypesIndexToRetain=[0,1,2,4,5,8]
        
        tempGenolistNoEquivalents=[i for i in tempGenolist if tempGenolist.index(i) in genotypesIndexToRetain] 
        dGL=dict(zip(tempGenolistNoEquivalents, GL.split(',')  ))    
        genotypeFormat=str(genotype[0])+str(genotype[2])
        likl=dGL[genotypeFormat]

        myout= [csqAllele, csqCount, likl]
    else: myout =False   ## no match between csq allele and alleles of teh vcf     
    return myout
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def checkFreq (listFreq, threshold): 
    # check if a variant is rare:  none of the populations in listfreq has frequency greater than threshold
    rare=True
    if len(listFreq) >0: 
        for freq in listFreq:
            if freq > threshold:
                rare=False
                break
    else: rare="NOB" #never observed 
    return(rare)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def CountCSQ_REF_ALT (csqAllele, refAllele, altAlleles, GTfields):
    """csqAllele = type: string, consequence allele  from vep  
    refAllele = type: string, reference allele from variant calling
    altAlleles = type: string, comma-separated list of alternate allele from variant calling
    nbAploidSamples = type: int, number of total Alleles
    GTfields = type: list, linesplit[9:] (take only from 9° column to the end) ["0/0:.:.:.:.:.:.:.","1/1:16:0,16:0:0:16:628:-34.365,-4.81648,0"]
    hetGenotypes = type: int, heterozygosity in samples
    countPassLine = type: int, line the are PASS in frequency analisys    """
    myres=False
    splitAltAlleles=altAlleles.split(",")
    allAlleles=[refAllele]+ splitAltAlleles
    mygstring=""; hetGenotypes=0; countPassLine=0;
    GTsplit=[i.split(":")[0] for i in GTfields]
    
    if GTsplit[0] != GTsplit[2]:
        hetGenotypes+=1

    for idx, item in enumerate(GTsplit):
        if item !='./.':#if i is not "./.":
            mygstring+=item; 
            if item[0] != item[1]: hetGenotypes+=1

    CountAlleles=[]
    for i in range(len(allAlleles)):
        #if str(i) in mygstring:
        CountAlleles.append(mygstring.count(str(i)))
    dAllele=dict(zip(allAlleles,CountAlleles))
        
    freqREF="{0:4.2f}".format(dAllele[refAllele])
        
        #### ALT freq 
    sumALT=0
    for i in splitAltAlleles:
        sumALT+=dAllele[i]
    freqALT="{0:4.2f}".format(sumALT)
    #### MAF freq
    MAF=min(freqALT,freqREF)
    #### CSQ freq
    if csqAllele in dAllele: 
        csqAllCount=dAllele[csqAllele]    
        freqCsq="{0:4.2f}".format(csqAllCount) 
            
    else: freqCsq='NA'

    #### define myres 
    myres= [freqCsq, freqREF, freqALT, MAF]
    return myres
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def AnnotateFreqCSQ_REF_ALT (csqAllele, refAllele, altAlleles, GTfields):
    """csqAllele = type: string, consequence allele  from vep  
    refAllele = type: string, reference allele from variant calling
    altAlleles = type: string, comma-separated list of alternate allele from variant calling
    nbAploidSamples = type: int, number of total Alleles
    GTfields = type: list, linesplit[9:] (take only from 9° column to the end) ["0/0:.:.:.:.:.:.:.","1/1:16:0,16:0:0:16:628:-34.365,-4.81648,0"]
    hetGenotypes = type: int, heterozygosity in samples
    countPassLine = type: int, line the are PASS in frequency analisys    """
    myres=False
    splitAltAlleles=altAlleles.split(",")
    allAlleles=[refAllele]+ splitAltAlleles
    mygstring=""; nbHaploidSamples=0; hetGenotypes=0; countPassLine=0;
    GTsplit=[i.split(":")[0] for i in GTfields]
    
    for idx, item in enumerate(GTsplit):
        if item !='./.':#if i is not "./.":
            mygstring+=item; nbHaploidSamples+=2
            if item[0]!=item[2]: hetGenotypes+=1
    CountAlleles=[]
    for i in range(len(allAlleles)):
        #if str(i) in mygstring:
        CountAlleles.append(mygstring.count(str(i)))
    dAllele=dict(zip(allAlleles,CountAlleles))
        
    if nbHaploidSamples!=0:
        #countPassLine+=1
        #### REF freq
        freqREF="{0:4.2f}".format(dAllele[refAllele]/nbHaploidSamples)
        
        #### ALT freq 
        sumALT=0
        for i in splitAltAlleles:
            sumALT+=dAllele[i]
        freqALT="{0:4.2f}".format(sumALT/nbHaploidSamples)
        
        #### MAF freq
        MAF=min(freqALT,freqREF)
        
        #### CSQ freq
        if csqAllele in dAllele: 
            csqAllCount=dAllele[csqAllele]    
            freqCsq="{0:4.2f}".format(csqAllCount/nbHaploidSamples) 
            # for obtain a '%' instead of integer
            # "{0:4.2f}%" : 
            # 4 = four number include the float ,
            # 2 = two float numbers,
            # f = for floating numbers; 
            # "{0:.0f}%" = if i want only '%' without float.
        
        else: freqCsq='NA'

        #### heterozygosity
        het= (hetGenotypes*2)/float(nbHaploidSamples)

        #### define myres 
        myres= [freqCsq, freqREF, freqALT, MAF, het]
    return myres
    #    else:
    #        myres= ['NA', freqREF, freqALT, MAF, hetGenotypes, countPassLine]
    #        return myres
    #else: 
    #    myres= ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']
    #    return myres
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### Functions END ###################################
###################################################




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="path to  input vcf file ",type=str,required=True)
    parser.add_argument("-j", help="path to  json  file ",type=str,required=True)
    parser.add_argument("-g", help="path to gene file list ",required=True)
    #parser.add_argument("-i", help="threshold for SOTerm Impact  ", type=int,required=True)
    parser.add_argument("-r", help="threshold for rare variant definition ", type=float,required=True) 
    parser.add_argument("-v", help="path to table of vep consequences  ",type=str, required= True)
    parser.add_argument("-p", help="path to table of pLI score  ",type=str, required= True)   
    parser.add_argument("-o", help="path to output file  ",type=str, required= True)
    #parser.add_argument("-st", help="path to stats file  ",type=str, required= True)
    parser.add_argument("-e", help="path to error file",type=str,required=True)
    parser.add_argument("-c", help="consequence allele count ",type=int,required=False, default=0)
    parser.add_argument("-w", help="path to  weights  file ",type=str,required=True)
    args = parser.parse_args()
        
         
    ##### 0b. read weights 
    dWeig={}
    for wline in open (args.w):
        w=wline.rstrip().split() 
        dWeig[w[1]]=int(w[2])
    #print (dWeights)
    
    ##### 1. parse vcf to produce dVcf[mykey]=[myref, myqual, dFormat["GT"]]; mykey is  1:333333:/T (T is the alternate allele) "
    
    dVcf={}
    for line in gzip.open(args.f, 'r'):
        decodedLine=line.decode()  ## line.decode() is necessary to read encoded data using gzip in python3
        if not re.match('#', decodedLine):
            linesplit=decodedLine.rstrip().split()
            mychr=linesplit[0]; mypos=linesplit[1]; myref=linesplit[3]; myalt=linesplit[4]; myqual=float(linesplit[5]); altAlleles=myalt.split(",")
            tempformattitle=linesplit[8].split(":")
            tempformatcontent=linesplit[9].split(":")
            dFormat=dict(zip(tempformattitle, tempformatcontent))

            for altAl in altAlleles:
                mykey=mychr.lstrip("chr") + ":" + mypos + ":/" + altAl
                dVcf[mykey]=[myref, myalt, myqual, dFormat["GT"]]   
        

    ##### 2. get VEP info from local json file 
    dV = getInfoFromVepLocally (args.j )   #yelds two dictionaries 
    dVepTrans=dV[0]; dVepCommon=dV[1]

    #~~csq allele features : number of allele with conseqces and genotype    
    for telem in dVepTrans:
        if 'csqAllele' in dVepTrans[telem]: 	
            vcfkey=telem[0]; csqAllele=dVepTrans[telem]['csqAllele']; altAllele=dVcf[vcfkey][1]; genotype=dVcf[vcfkey][3]
            altAlleleCount=2-genotype.count('0')
            csqCount= csqAlleleCount (csqAllele, altAllele, altAlleleCount, 2)
            dVepTrans[telem]['csqCount']=csqCount  #1 het #2 homo

    #~~check if variant is rare 
    for celem in dVepCommon:
        if 'frequencies' in dVepCommon[celem]: 
            for cAll in dVepCommon[celem]['frequencies']: 
                listFreq=dVepCommon[celem]['frequencies'][cAll].values()
                rare=checkFreq (listFreq, args.r )
                dVepCommon[celem]['rare']=rare
                dVepCommon[celem]['commonCSQall']=cAll
    
    #~~create dataframe from dV
    dfTrans = pd.DataFrame(dVepTrans).T 
    dfCommon= pd.DataFrame(dVepCommon).T
    
    
    
    #dfTrans.to_csv("/lustrehome/silvia/cicci.trans.tsv", sep="\t")#, index=True )
    #dfCommon.to_csv("/lustrehome/silvia/cicci.common.tsv", sep="\t")#, index=True )
    df=dfTrans.set_index('key').join(dfCommon, how="left"  )
    #df.to_csv("/lustrehome/silvia/cicci.ttcos.tsv", sep="\t")
	
	
    '''
                     most_severe_consequence            id csqAllele          gene_id gene_symbol gnomad_nfe gnomad_eas gnomad_asj ....
    5:64548:/A        intergenic_variant           NaN       NaN              NaN         NaN        NaN        NaN        NaN ....
    1:10623:/C     upstream_gene_variant  rs1207502826         C  ENSG00000223972     DDX11L1        NaN        NaN        NaN ....
    1:95068:/A            intron_variant   rs867757274         A  ENSG00000238009  AL627309.1        NaN        NaN        NaN ....
    1:942451:/C         missense_variant     rs6672356         C  ENSG00000187634      SAMD11     0.9999          1          1 ....
    1:1519044:/T          intron_variant   rs147532057         T  ENSG00000197785      ATAD3A        NaN        NaN        NaN ....
    ''' 

    #### 3. info form gene lists
     
    gene_list = pd.read_csv(args.g,sep="\t")
    genes=gene_list[["ensID", "final_score"]]
    df_last =df.reset_index().merge(genes,left_on="gene_id",right_on="ensID",how="left").drop("ensID", axis=1)
    df_last.rename({"final_score":"gene_score"},axis=1,inplace=True)
 
    '''
    def convert_string_to_array(x):
        tree = ast.parse(x, mode='eval')
        transformer = Transformer()
        # raises RuntimeError on invalid code
        transformer.visit(tree)
        # compile the ast into a code object
        clause = compile(tree, '<AST>', 'eval')
        # make the globals contain only the Decimal class,
        # and eval the compiled object
        return eval(clause, dict(Decimal=decimal.Decimal))
    
    try:
        df.loc[:,"gene_id"] = df.gene_id.apply(convert_string_to_array)
    except Exception as e:
        pass

    df.loc[:,"score_gene"] = df.gene_id.apply(lambda x: gene_list[gene_list.ensID==x])
    '''
    #### 3. SOScore 
    dSOTermFineRank=VepRankingInfo(args.v)
    soScore = pd.Series(dSOTermFineRank,name="soScore").to_frame().reset_index()
    df_final = df_last.reset_index().merge(soScore,left_on="most_severe_consequence",right_on="index").set_index("index_x").drop("index_y",axis=1)
    #df_final.to_csv(args.o,sep="\t",index=True)


    #### 4. CADDD 
    #final score
    #df_last.loc[:,"gpScore"] = (df_last.csqCount.astype(float) * dWeig[wCAC]) + (df_last.rare.astype(float) * dWeig[wRare]) + (df_last.soScore.astype(float) * dWeig[wRank]) + df_last.score_gene_list.astype(float)
    ##dask_df = dd.read_csv('/lustre/home/enza/CADD/SNV_ch*.tsv',sep="\t")
    #creation of the key
    ##dask_df["key"] = dask_df["#Chrom"].map(str) +":"+ dask_df["Pos"].map(str) +":"+"/"+ dask_df["Alt"]
    ##dask_df = dask_df.drop(["#Chrom","Pos","Ref","Alt","PHRED"],axis=1)
    ##dask_df = dask_df.set_index("key")
    ##merged = dd.merge(df_last, dask_df, left_index=True, right_index=True,how="left") 

    ##merged.to_csv(args.o+"*.tsv",sep="\t",index=True)
    #df_for_stats.to_csv(args.st,sep="\t",index=True)

    #### 5. pLI 
    ### load pli table
    pLI_score = pd.read_csv(args.p,sep="\t")
    pliScore=pLI_score[["transcript", "pLI"]]
    df_info = df_final.reset_index().merge(pliScore,left_on="element_id",right_on="transcript",how="left").drop("transcript", axis=1)
    df_info.to_csv(args.o,sep="\t",index=True)
    ### add pli to df 


if __name__ == "__main__":
    main()
