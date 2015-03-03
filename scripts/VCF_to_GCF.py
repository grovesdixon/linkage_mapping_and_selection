#!/usr/bin/env python
##VCF_to_OneMap_input.py
##written 11/6/13 by Groves Dixon
ProgramName = 'VCF_to_GCF.py'
LastUpdated = '1/24/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
#global variables are capitalized and camel case: GlobalVaraibleName
#function variables are not capitalized and camel case: functionVariableName
#function names separated by underscores: function_to_do_a_thing()


##Assign Global User Help Variables

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script is intended to convert VCF files into a simple 
variant call text file. This will allow simplified custom
filtering and conversion to other input file types such as
OneMap inputs.
'''

AdditionalProgramInfo = '''
Additional Program Information:
Importantly, this script will 'sort' genotype data,
so if heterozygotes with alleles listed in different orders
will be normalized. For example A/G and G/A will both be
recorded as A/G in the output.
'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-input', '-i', required = True, metavar = 'inputName.vcf', dest = 'In', help = 'The name of the input vcf file')
parser.add_argument('-output', '-o', required = True, metavar = 'outputName.out', dest = 'Out', help = 'The desired name for the output file')
args = parser.parse_args()

#Assign Arguments
InfileName = args.In
OutfileName = args.Out

#Assign some global Variables
RequiredHeaderParts = ['#CHROM', 'POS', 'REF', 'ALT', 'INFO', 'FORMAT']


#Functions

def read_file(InfileName, Delim):
    '''Function to read in a file as a list of lists
    '''
    lineNumber = 0
    fileList = []
    with open(InfileName, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split(Delim)
            ##skip lines starting with ##
            if line[0][0:2] == '##':
                continue
            #record the header
            if line[0][0] == '#':
                header = line
                continue
            fileList.append(line)
    return header, fileList

def parse_heading(Header, Required):
    '''Function to check the header and gather metadata.
    Looks through the header and make sure
    all the necessary header elements are present. Returns the 
    indices for the required header elements.
    '''
    headerIndices = []
    for i in Required:
        if i not in Header:
            exit('\nERROR:\nThe required component {} was not found in header. Please check input file.\n'.format(i))
        else:
            headerIndices.append(Header.index(i))
    print('\nChecking that header has necessary components...\n')
    print('\nAll necessary header components have been found\n')
    print('\nHeader indices are as follows:')
    for x in range(len(Required)):
        print('{}: {}'.format(Required[x], str(headerIndices[x])))
    chromIndex = Header.index('#CHROM')
    infoIndex = Header.index('INFO')
    posIndex = Header.index('POS')
    refIndex = Header.index('REF')
    altIndex = Header.index('ALT')
    formatIndex = Header.index('FORMAT')
    return chromIndex, infoIndex, posIndex, refIndex, altIndex, formatIndex

def get_sample_list(Header):
    '''The samples are included in the header in VCF files.
    This function pulls a list of just the samples from the header.
    Also checks the sample list to make sure parents are present.
    '''
    start = Header.index('FORMAT')+1 ##in VCF samples always start after FORMAT
    sampleList = Header[start:]
    print('\nNumber of Samples Found in VCF: {}\n'.format(len(sampleList))) 
    ##check for parents and exit if not there
    return sampleList
        
def get_variant_list(FileList, header):
    '''Function pulls two lists from the vcf filelist:
    a list of all the variants recorded as CHROM_POS
    and a coindexed list with the INFO column for each variant
    Indexing in this list will be used to keep track of which variants
    the sample genotypes are connected to.
    Also outputs the reference allele list, indexed the same as
    the variant list.
    '''
    #create two empty lists to store info in
    variantList = []
    infoList = []
    referenceAlleleList = []
    #being loop through data
    for line in FileList:
        chrom = line[ChromIndex] ; pos = line[PosIndex] ; info = line[InfoIndex] ##using indices assinged by parse_header()
        #reformat the variant as 'CHROM_POS'
        variant = '{}_{}'.format(chrom,pos)
        ##add the variant and info for this line to the two lists
        infoList.append(info)
        variantList.append(variant)
        ref = line[RefIndex]
        referenceAlleleList.append(ref)
    numberVariants = len(variantList)
    print('\n{} Variants Found in the VCF\n'.format(numberVariants))
    print
    return variantList, referenceAlleleList, infoList
    
def gather_genotypes(fileList, sampleList, variantList, header):
    '''This function iterates through the file list and records
    the genotype data for each sample in a dictionary. The genotype data for each 
    variant position is recorded as a list for the sample. The list for
    each sample coordinates with the VariantList. SNP calls 
    are stored in dictionary as A,T,G,and C. The base calls are retrieved using
    metadata gathered by parse_heading().
    
    Dictionary return format:
    {sample1Name: [['A','T'],['G','G'],['-','-']...], sample2Name: [[etc.]]}
    where the list of genotypes associated with the sampleName is equivalent 
    to the list of variant positions. *Note that the genotypes are reordered
    alphabetically so that all genotypes with the same two alleles will be
    identical lists
    '''
    #empty dictionary to store genotype data
    genotypeDict = {}
    #make a key for each sample with list to store genotype data
    for sample in sampleList:
        genotypeDict[sample] = [] ##list to store genotypes for each variant
    lineNumber = 0
    ##begin loop through each line of data
    for line in fileList:
        lineNumber += 1
        variant = VariantList[lineNumber-1]
        ##get ref and alt based on indices output by parse_heading()
        ref = line[RefIndex] #single base
        alt = line[AltIndex].split(',')##can have multiple alternates delimited by ','
        ##get the index for the genotype data from the metadata provided in FORMAT.
        genotypeIndex = line[FormatIndex].split(':').index('GT')##FORMAT is colon delimited. GT marks index for genotype.
        samplesStart = FormatIndex+1 ##samples always begin after FORMAT column in VCF
        ##get slice of the line with the genotype data
        genotypeData = line[samplesStart:]
        #set up counter to keep track of index
        sampleIndex = -1
        ##iterate through the genotype data in the line
        for i in genotypeData:
            i = i.split(':') ##the genotype entries (like the FORMAT elemnt) is colon delimited
            #keep track of index
            sampleIndex += 1
            #get sample name from sampleList using index
            sampleName = sampleList[sampleIndex]
            ##use the genotypeIndex to grab the genotype data
            genotype = i[genotypeIndex]
            
            if '|' in genotype:
                genotype = genotype.split('|') ##genotypes can be split with '|' or '/' for phased and unphased data repsectively
            elif '/' in genotype:
                genotype = genotype.split('/')
            else:
                print('/nunrecognized genotype delimiter\
                for sample {} for variant number {}'.format(sampleName, variantList[lineNumber-1]))
                continue
            #genotype is now recorded as a list of integers
            #convert the integers into A T G C or -
            convertedGenotype = [] #variable to store the genotype in converted form
            #loop through the two alleles for this genotype
            
            # print
            # print
            # print variant
            # print sampleName
            # print 'vcf genotype:'
            # print genotype
            # print 'ref:'
            # print ref
            # print 'alt:'
            # print alt
            
            ##convert the genotype into normal letters
            for allele in genotype:
                #convert for missing data
                if allele == '.': ## . marks missing data in VCF file
                    convertedGenotype.append('-') ## - marks missing data in OneMap
                #convert for present data
                else:
                    allele = int(allele)##alleles are stored as numbers in VCF. Make integer for indexing
                    #convert if it is reference allele
                    if allele == 0:
                        convertedGenotype.append(ref)## zero indicates the allele is the reference
                    #convert if it is an alternate allele 
                    else:
                        convertedGenotype.append(alt[(allele-1)]) # nonzero tells which of the alternate alleles it is (indexed from 1 in VCF so subtract for python)
            #done converting the alleles for the genotype
            #add the converted genotype to the correct sample's list in the dictionary
            convertedGenotype.sort()
            genotypeDict[sampleName].append(convertedGenotype) ##sorting here insures all genotypes with same alleles are identical
    return genotypeDict
        
def output(GenotypeDict, SampleList, VariantList, OutfileName, InfileName):
    with open(OutfileName, 'w') as out:
        out.write('#data from {} as simple genotype table'.format(InfileName))
        out.write('\n*\t' + '\t'.join(SampleList)) 
        for i in range(len(VariantList)):
            genotypeList = []
            variant = VariantList[i]
            for x in SampleList:
                genotype = GenotypeDict[x][i]
                genotype = '/'.join(genotype)
                genotypeList.append(genotype)
                variantString = '\t'.join(genotypeList)
            out.write('\n' + variant + '\t' + variantString)
                            
def output_references(variantList, refList, refOutName):    
    with open(refOutName, 'w') as out:
        out.write('variant\tref')
        for i in range(len(variantList)):
            outString = '\n{}\t{}'.format(variantList[i], refList[i])
            out.write(outString)
            

Header, FileList = read_file(InfileName, '\t')
SampleList = get_sample_list(Header)
ChromIndex, InfoIndex, PosIndex, RefIndex, AltIndex, FormatIndex = parse_heading(Header, RequiredHeaderParts)
VariantList, ReferenceAlleleList, InfoList = get_variant_list(FileList, Header)
GenotypeDict = gather_genotypes(FileList, SampleList, VariantList, Header)
# print GenotypeDict["AC1c_2_TACGTG_CATC.trim.bt2"]
# print len(GenotypeDict["AC1c_2_TACGTG_CATC.trim.bt2"])
output(GenotypeDict, SampleList, VariantList, OutfileName, InfileName)
RefOutName = InfileName[0:-4] + '_ReferenceAlleles.txt'
output_references(VariantList, ReferenceAlleleList, RefOutName)

#return time to run
Time = time.time() - Start_time
print('Time took to run: {}'.format(Time))


