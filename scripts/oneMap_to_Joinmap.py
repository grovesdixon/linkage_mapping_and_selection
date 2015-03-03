#!/usr/bin/env python
##NAMEOFTHEPROGRAM.py
##written 6/26/14 by Groves Dixon
ProgramName = 'oneMap_to_joinMap.py'
LastUpdated = '6/26/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This program takes a oneMap Input file and converts it into a 
JoinMap input file.

'''

AdditionalProgramInfo = '''
Additional Program Information:


'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = True, dest = 'input', help = 'The the OneMap input file')
parser.add_argument('-s', required = True, dest = 'samples', help = 'The the Samples file that has the sample names the match the genotypes in teh OneMap Input file. (This was output by GCF_to_OneMap)')
parser.add_argument('-p', required = False, dest = 'parents', help = 'The the Concensus parent Genotypes file. (This was also output by GCF_to_OneMap)')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-n', required = True, dest = 'name', help = 'The desired name for the data file (included in the header of the JoinMap loc file output. Probably something informative to distinguish this file from other similar ones including the parameters used to set it up)')
parser.add_argument('-popType', required = False, default = "CP", dest = 'popType', help = 'The name for the population (included in the header of the JoinMap loc file.)')
args = parser.parse_args()


#Assign Arguments
InfileName = args.input
OutfileName = args.out
sampleFileName = args.samples
Name = args.name
PopulationType = args.popType



def build_crossType_conversion_dict():
    """Sets up the cross type conversions
    to go from onemap to JoinMap"""
    conversionDict = {
        "A.1" : "<abxcd>",
        "A.2" : "<efxeg>",
        "B3.7" : "<hkxhk>",
        "D1.9" : "none",
        "D1.10" : "<lmxll>",
        "D2.14" : "none",
        "D2.15" : "<nnxnp>"
        }
    return conversionDict

def read_onemap_input(InfileName, conversionDict):
    '''Function to read the oneMap input file and 
    arrange the data into dictionaries.
    '''
    lineNumber = 0
    variantList = [] 
    crossTypeDict = {} #store the crosstypes keyed to marker names
    conCrossTypeDict = {}
    genotypeDict = {} #store the genotypes as lists keyed to marker names
    skippedCount = 0
    with open(InfileName, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split()
            if lineNumber == 1:
                sampleCount = int(line[0])
                variantCount = int(line[1])
                continue
            marker = line[0][1:]
            crossType = line[1]
            if crossType == "D2.14":
                skippedCount += 1
                continue
            if crossType == "D1.9":
                skippedCount += 1
                continue ##skip the cross types that don't work in JoinMap
            convertedCrossType = conversionDict[crossType]
            genotypes = line[2]
            genotypeList = genotypes.split(",")
            if len(genotypeList) != sampleCount:
                print "WARNING...Number of genotypes does not match sample total given in the heading of the OneMap file"
            crossTypeDict[marker] = crossType
            genotypeDict[marker] = genotypeList
            conCrossTypeDict[marker] = conversionDict[crossType]
            variantList.append(marker)
    if (len(crossTypeDict.keys()) + skippedCount) != variantCount:
        print "WARNING...The Number of Variants in the OneMap input file does not match the value given in the header"
    return conCrossTypeDict, crossTypeDict, genotypeDict, variantList

    


def build_genotype_conversion_dict():
    """Function that sets of the genotype annotation
    conversion file. Outputs a nested dictionary
    first with JoinMap style cross types as keys connecting
    to dictionaries of the genotype annotation conversions for
    that crosstype. Note JoinMap manual says order alleles 
    are given in doesn't matter. So I just mirrored Onemap.
    Also note JoinMap doesn' seem to have equivalents for
    corssTypes D1.9 or D2.14"""
    genoConvertDict = {
        "<abxcd>" : {#(cross type A.1) this are actually the same but I keyed then in just to make it simple
            'ac' : 'ac',
            'ad' : 'ad',
            'bc' : 'bc',
            'bd' : 'bd',
            '-'  : '--'
        },
        "<efxeg>" : {
            'a'  : 'ee',
            'ac' : 'eg',
            'ba' : 'fe', 
            'bc' : 'fg',
            '-'  : '--'
        },
        "<hkxhk>" : {
            'a'  : 'hh',
            'ab' : 'hk',
            'b'  : 'kk',
            '-'  : '--'
        },
        "<lmxll>" : {
            'a'  : 'll',
            'ab' : 'ml',
            '-'  : '--'
        },
        "<nnxnp>" : {
            'a'  : 'nn',
            'ab' : 'np',
            '-'  : '--'
        }
    }
    return genoConvertDict

def read_samples(sampleFileName, genotypeDict, variantList):
    """Reads in the sample names"""
    with open(sampleFileName, 'r') as infile:
        for line in infile:
            line = line.strip("\n").split()
            samples = line[2].split(",")
    if len(samples) != len(genotypeDict[variantList[1]]):
        print len(samples)
        print len(genotypeDict[variantList[1]])
        print "WARNING...The number of samples in {} does not match with the number of genotypes in each row".format(sampleFileName)
        exit()
    return samples



def convert_genotypes(variantList, conCrossTypeDict, genoConvertDict, genotypeDict):
    """Function to convert the genotype data from onemap to 
    Joinmap notation. Uses the dictionary 'conversionDict' created above"""
    joinMapGenotypeDict = {}
    for i in variantList:
        crossType = conCrossTypeDict[i]
        joinMapGenotypeDict[i] = [] #list to store the genotypes for this variant
        for geno in genotypeDict[i]:
            joinMapGeno = genoConvertDict[crossType][geno]
            # print
            # print geno
            # print joinMapGeno
            joinMapGenotypeDict[i].append(joinMapGeno)
    return joinMapGenotypeDict
        

def output(joinMapGenotypeDict, conCrossTypeDict, variantList):
    """Outputs the data as a JoinMap Input File"""
    with open(OutfileName, 'w') as out:
        date = time.strftime("%m/%d/%Y")
        timeOday = time.strftime("%H:%M:%S")
        header = ''';JoinMap Loc File Converted from {}
;Converted on {} at {}\n
name = {}
popt = {}
nloc = {}
nind = {}\n'''.format(InfileName, date, timeOday, Name, PopulationType, len(variantList), len(genotypeDict[variantList[1]]))
        print "This is the header to be used in the output file:"
        print header
        out.write(header)
        for i in variantList:
            crossType = conCrossTypeDict[i]
            genotypeList = joinMapGenotypeDict[i]
            outString1 = "\n{}\t{}\n".format(i, crossType)
            outString2 = "\t" + '\t'.join(genotypeList)
            out.write(outString1)
            out.write(outString2)
            

conversionDict = build_crossType_conversion_dict()
conCrossTypeDict, crossTypeDict, genotypeDict, variantList = read_onemap_input(InfileName, conversionDict)
genoConvertDict = build_genotype_conversion_dict()
# sampleList = read_samples(sampleFileName, genotypeDict, variantList)
joinMapGenotypeDict = convert_genotypes(variantList, conCrossTypeDict, genoConvertDict, genotypeDict)
output(joinMapGenotypeDict, conCrossTypeDict, variantList)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


