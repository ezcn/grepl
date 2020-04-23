#!/usr/bin/python3
import pandas as pd
import glob, argparse,  subprocess, random 


def Freq_CSQ_REF_ALT (csqAllele, refAllele, altAlleles, missing_data_format, genotypeslist):
    """csqAllele = type: string, consequence allele  from vep  
    refAllele = type: string, reference allele from variant calling
    altAlleles = type: string, comma-separated list of alternate allele from variant calling
    nbAploidSamples : calculated 
    GTfields = type: list, list of genotypes ["0/0", "0/1", ".", "./."]
    hetGenotypes = type: int, heterozygosity in samples"""
    myres=False
    listAlt=altAlleles.split(",")   
    listAll=[refAllele]+ listAlt
    stringOfGenotypes=""; nbHaploidSamples=0
    for item in (genotypeslist):    
        if item != missing_data_format: 
            stringOfGenotypes+=item; nbHaploidSamples+=2
    #return stringOfGenotypes,  nbHaploidSamples
  
    CountAlleles=[]
    for i in range(len(listAll)):  # 0 for REF, 1 for ALT1, 2 for ALT2 ...
        CountAlleles.append(stringOfGenotypes.count(str(i)))
    dAllele=dict(zip(listAll,CountAlleles))

    if nbHaploidSamples!=0:
        freqREF="{0:4.2f}".format(dAllele[refAllele]/nbHaploidSamples)
        freqAlt=[]
        for i in listAlt: 
            freqAltTemp="{0:4.2f}".format(dAllele[i]/nbHaploidSamples)
            freqAlt.append(freqAltTemp)
        #~~~ CSQ 
        if csqAllele in dAllele:
            csqAllCount=dAllele[csqAllele]
            freqCsq="{0:4.2f}".format(csqAllCount/nbHaploidSamples) 
        else: freqCsq='NA'
        myres= [freqCsq, freqREF, freqAlt]
    return myres
