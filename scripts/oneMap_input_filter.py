#!/usr/bin/env python
##oneMap_input_filter.py
##written 1/10/14 by Groves Dixon
# outline of how it works at bottom of program
ProgramName = 'oneMap_input_filter.py'
LastUpdated = '6/23/14'
#also double checked for accuracy on 5/23/14
By = 'Groves Dixon'
VersionNumber = '4.0'

##Assign Global User Help Variables

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script is intended to filter a file of genotypes
in OneMap input format. It filters based on genotype frequencies (agreement with Mendelian Expectations)
and based on absent data.
'''

AdditionalProgramInfo = '''
This is written expecting that the input file was
generated with the script 'VCF_to_OneMap_input.py'. 
This script has already filtered based on 
parental genotypes and agreement between parental replicates.

This script does not require parental genotype information
because it assumes that the cross type information included in 
the input file are accurate.

Changes as of 5/24/14:
Added a feature where the samples names for the genotypes are
retained in a separate 'header' file. This is also input
so that Mendelian filtering can be done using only controls,
or any selected subset of samples.

Changes 6/24/14:
Switched back to normal p values instead of adjusted.
'''

#import modules
import time
import argparse
from sys import argv
from sys import exit
import sys
import numpy
import scipy.stats
StartTime = time.time()

#set up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-input', '-i', required = True, metavar = 'inputName.vcf', dest = 'In', help = 'The name of the input vcf file')
# args = parser.parse_args()
parser.add_argument('-output', '-o', default='oneMap_input_filter_output.txt', metavar = 'outputName.out', dest = 'Out', help = 'The desired name for the output file')
parser.add_argument('-pCut', '-p', default=0.2, type = float, metavar = 'pvalue_threshold', dest = 'pCut', help = 'The P value threshold to use for eliminating variants based on divergence from Mendelian expectations. Default is 0.2. Note that this p value refers to the adjusted value returned after FDR correction')
parser.add_argument('-sCut', '-s', default = 0.8, type = float, metavar = 'sample_proportion_genotyped_theshold', dest = 'sCut', help = 'The cutoff for the proportion of variants that must be correctly genotyped below which a sample will be removed')
parser.add_argument('-vCut', '-v', default = 0.8, type = float, metavar = 'variant_proportion_genotyped_theshold', dest = 'vCut', help = 'The cutoff for the proportion of samples that must be genotyped for a variant for it to be kept')
parser.add_argument('-Samples', dest = 'Header', help = 'The name of a header file output from GCF_to_OneMap_v2.py. It will list the sample names in the same format as the onemap input file.')
parser.add_argument('-M', dest = 'Select', help = 'The file name of a table listing the samples to use for the Mendelian Filter')
args = parser.parse_args()
InfileName = args.In
OutfileName = args.Out
PvalueCutoff = args.pCut
SampleGenotypedCutoff = args.sCut
VariantGenotypedCutoff = args.vCut
HeaderFile = args.Header
SelectSamplesFile = args.Select


# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector
# stats = importr('stats')


#set up some global variables
CrossTypeLabels = ['A.1','A.2','B3.7', 'D1.9', 'D1.10', 'D2.14', 'D2.15']

##Assign Functions
def read_file(InfileName):
    '''function to read in the genotype file
    returns the data in three lists:
    genotypeList, crossTypeList, and variantList
    The indices of these lists correspond to the variants
    The indices of the lists nested within genotypeList correspond
    to the samples.
    '''
    genotypeList = []
    variantList = []
    crossTypeList = []
    with open (InfileName, 'r') as infile:
        for line in infile:
            if line[0] != '*':
                continue
            else:
                line = line.strip('\n').split()
                variantList.append(line[0])
                crossTypeList.append(line[1])
                genotypes = line[2].split(',')
                genotypeList.append(genotypes)
    variantCount = len(variantList)
    sampleCount = len(genotypeList[0])
    print
    print 'Number of Samples in File = {}'.format(sampleCount)
    print
    print 'Number of Variants Detected in File = {}'.format(variantCount)
    return genotypeList, variantList, crossTypeList

def build_expected_ratio_dict():
    '''function to build the dictionary connecting cross types with their expected
    ratio vectors. Data will be recorded in nested dictionaries. Top set of keys
    is the cross types. For each cross type, the next set of keys is the possible genotypes
    linked with their expected ratios based on the cross type. See OneMap tutorial for table
    of cross type notations and expected ratios.
    '''
    expRatiosDict = {'A.1' : {'ac' : .25, 'ad' : .25, 'bc' : .25, 'bd' : .25},
    'A.2' : {'a' : .25, 'ac' : .25, 'ba' : .25, 'bc' : .25},
    'B3.7' : {'a' : .25, 'ab' : .5, 'b' : .25},
    'D1.9' : {'ac' : .5, 'bc' : .5},
    'D1.10' : {'a' : .5, 'ab' : .5},
    'D2.14' : {'ac' : .5, 'bc' : .5},
    'D2.15' : {'a' : .5, 'ab' : .5}}
    return expRatiosDict

def build_counting_dict():
    '''function to build the dictionary connecting cross types with their expected
    ratio vectors. Data will be recorded in nested dictionaries. Top set of keys
    is the cross types. For each cross type, the next set of keys is the possible genotypes
    linked with their expected ratios based on the cross type. See OneMap tutorial for table
    of cross type notations and expected ratios.
    '''
    countingDict = {'A.1' : {'ac' : 0, 'ad' : 0, 'bc' : 0, 'bd' : 0},
    'A.2' : {'a' : 0, 'ac' : 0, 'ba' : 0, 'bc' : 0},
    'B3.7' : {'a' : 0, 'ab' : 0, 'b' : 0},
    'D1.9' : {'ac' : 0, 'bc' : 0},
    'D1.10' : {'a' : 0, 'ab' : 0},
    'D2.14' : {'ac' : 0, 'bc' : 0},
    'D2.15' : {'a' : 0, 'ab' : 0}}
    return countingDict

def remove_bad_genotypes(genotypeList, crossTypeList, variantList, expRatioDict):
    '''This function iterates through all of the genotypes
    and eliminates those that are not listed as possible genotypes for that
    cross type. For example, cross type D1.10 is ab x cc. The genotype 'cc' could 
    not occur here. Any instance of such a genotype will be changed to missing data
    for that sample for that variant'''
    print
    print 'Testing for impossible offspring genotypes...'
    impossibleCounter = 0
    goodCounter = 0
    impossibleCalls = []
    for i in range(len(crossTypeList)):
        crossType = crossTypeList[i]
        variant = variantList[i]
        genotypes = genotypeList[i]
        for x in range(len(genotypes)):
            sampleGeno = genotypes[x]
            if sampleGeno == '-':
                continue
            elif sampleGeno not in expRatioDict[crossType].keys():
                impossibleCounter += 1
                instance = '_'.join([crossType, variant, sampleGeno, 'sampleNumber: {}'.format(i)])
                impossibleCalls.append(instance)
                genotypeList[i][x] = '-' #replace the genotype with missing data
            else:
                goodCounter += 1
                continue
    print
    print '{} genotype calls were not possible based on crossType and replaced as missing data'.format(impossibleCounter)
    print
    percent = float(impossibleCounter) / (len(variantList) * len(genotypeList[0])) * 100
    print '{}% of genotype calls were not possible'.format(percent)
    return genotypeList
        

def read_header(HeaderFile):
    """Function to read in the header file to get the list of 
    sample sames that is coindexed with the genotype list for each 
    variant"""
    lineNumber = 0
    with open(HeaderFile, 'r') as infile:
        for line in infile:
            lineNumber +=1
            if lineNumber == 1:
                line = line.strip('\n').split('\t')[2].split(',') ##parse the list of sample names from the Header file
                sampleList = line
    return sampleList
            

def remove_poorly_typed_samples(sampleCount, markerCount, sampleGenotypedCutoff, genotypeList, sampleList):
    """Iterate through the samples' genotypes and count up how many are missing. Builds a list of sample numbers
    (as indices to genotypeList) that will be removed.
    Output is a new filtered list of genotypes, a list of cut sample indices, and a new filtered sample list.
    newGenotypeList is a list of lists, with the parent list coindexed with variantList, and the embedded lists
    coindexed with sampleList:
    newGenotypeList = [[list of genos for variant1 coindexed with sampleList], [list of genos for variant2]...[list of genos for variant(len(variantList))]]"""
    missingList = []
    cutList = []
    cutCount = 0
    newGenotypeList = []##empty list to store genotypes that pass filter
    newSampleList = []##empty list to store sample names to match the genotypes that pass filter. The two lists will be coindexed
    for i in range(sampleCount):
        missingList.append(0) #set up the list with zeros initially to add missing data counts to. Every entry here
    for i in range(markerCount):
        genotypes = genotypeList[i]
        for x in range(len(genotypes)):
            if genotypes[x] == '-':
                missingList[x] += 1
    for z in range(len(missingList)):
        if missingList[z] > markerCount - (sampleGenotypedCutoff * markerCount):
            cutCount += 1
            cutList.append(z)
    print
    print '{} samples had total genotypes below theshold and will be removed'.format(cutCount)
    for i in range(markerCount):
        newGenotypes = []
        newHeader = []
        genotypes = genotypeList[i]
        for x in range(len(genotypes)):##update the newGenotypeList with the filtered set of genotypes
            if x not in cutList:
                newGenotypes.append(genotypes[x])
        newGenotypeList.append(newGenotypes)
    for x in range(len(genotypes)): #use the same cutList to update the sampleList to stay coordinated with teh filtered set of genotypes
        if x not in cutList:
            newSampleList.append(sampleList[x])
    return newGenotypeList, cutList, newSampleList #
        
        
def remove_low_covered_variants(genotypeList, crossTypeList, variantList, variantGenotypedCutoff):
    nGenotypeList = []
    nCrossTypeList = []
    nVariantList = []
    cutCount = 0
    for i in range(len(variantList)):
        missingCount = 0
        for z in genotypeList[i]:
            if z == '-': missingCount += 1
        sampleCount = len(genotypeList[i])
        if missingCount < sampleCount - VariantGenotypedCutoff * sampleCount:
            nGenotypeList.append(genotypeList[i])
            nCrossTypeList.append(crossTypeList[i])
            nVariantList.append(variantList[i])
        else:
            cutCount += 1
            continue
    print
    print '{} variants were removed because too few samples were genotyped'.format(cutCount)
    return nGenotypeList, nCrossTypeList, nVariantList
                
def get_indices_for_mendelian_filter(SampleList, SelectSamplesFile):
    """This function gets a list of indices matching the sample
    names in the Header file. The indices can then be used to skip
    particular samples when doing the Mendelian filtering below."""
    print "\nReading the Select Set of Samples to use for Mendelian Filtering..."
    selectList = []
    selectIndexList = [] #set up a list to store the indices of the select samples to call them from SampleList
    selectCount = 0
    nonSelectCount = 0
    with open(SelectSamplesFile, 'r') as infile:
        for line in infile:
            line = line.strip('\n')
            selectList.append(line)
    for i in range(len(SampleList)):
        if SampleList[i] in selectList:
            selectCount += 1
            selectIndexList.append(i)
        else:
            nonSelectCount += 1
    print "\nFilter Passing Sample Count = {}".format(len(SampleList))
    print "\nFound {} Samples in the Select File".format(len(selectList))
    print "\n{} Filter passing Samples from sampleList were not found in the Select File".format(str(nonSelectCount))
    print "\nUsing {} Select Samples from the Sample List for Mendelian Filtering".format(len(selectIndexList))
    return selectIndexList
            
        
        
    

def remove_nonMendelian_variants(genotypeList, crossTypeList, variantList, countingDict, expRatiosDict, PvalueCutoff, SelectIndexList):
    """function that iterates through the variants and looks at the segregation ratios to determine if the variant has significant distortion.
    Variants with p values below the set threshold are removed. Calls on expRatiosDict, which gives the expected ratios of each genotype for
    each cross type.
    """
    nGenotypeList = []
    nCrossTypeList = []
    nVariantList = []
    pValueList = []
    cutCount = 0
    for i in range(len(variantList)):
        crossType = crossTypeList[i] #get the cross type for this variant
        genotypes = [] ##set up empty list to store the genotypes that match the selectedSamples from the selectFile
        missingCount = 0
        for INDEX in SelectIndexList: #iterate through the set of selected samples' indices
            genotypes.append(genotypeList[i][INDEX]) #add only the genotypes corresponding to selected samples to the list for doing the mendelian chi square test
        countDict = countingDict[crossType]
        for z in countDict.keys():
            countDict[z] = 0
        for z in genotypes:
            if z == '-':
                missingCount += 1
                continue
            countDict[z] += 1
        types = countDict.keys()
        sampleCount = len(genotypes) - missingCount
        observed = []
        expected = []
        for x in types:
            observed.append(countDict[x])
            expected.append(expRatiosDict[crossType][x] * sampleCount)
        observedArr = numpy.array(observed)
        expectedArr = numpy.array(expected)
        x2, P = scipy.stats.chisquare(observedArr, expectedArr)
        # #####---------- debugging
        # print '\n\n\n{}'.format(variantList[i])
        # # for G in range(len(SampleList)):
        # #     if G in SelectIndexList:
        # #         print SampleList[G] + 'IN'
        # #     else:
        # #         print SampleList[G] ######THIS SHOWS THAT YOUR PULLING THE RIGHT SAMPLES
        # #########
        # for G in range(len(SelectIndexList)):
        #     print genotypes[G]
        # print observedArr
        # print expectedArr
        # print P
        # print '----------------------'
        # ####################------- done
        pValueList.append(float(P))
    meanP = numpy.mean(pValueList)
    # p_adjust = stats.p_adjust(FloatVector(pValueList), method = 'BH') ##use rpy2 to call on R to do the FDR correction
    for i in range(len(pValueList)):
        P = pValueList[i]
        if P > PvalueCutoff:
            # print 'KEEP'
            nGenotypeList.append(genotypeList[i])
            nCrossTypeList.append(crossTypeList[i])
            nVariantList.append(variantList[i])
        else:
            # print 'CUT'
            cutCount += 1
            continue
    print
    print 'Mean P value across all variants:'
    print meanP
    print
    print '{} variants were removed because they did not conform to Mendelian Expectations'.format(cutCount)
    return nGenotypeList, nCrossTypeList, nVariantList
            
def output(outfileName, genotypeList, crossTypeList, variantList):
    with open(outfileName, 'w') as out:
        sampleCount = len(genotypeList[0])
        variantCount = len(variantList)
        print
        print 'New Filtered Sample Total = {}'.format(sampleCount)
        print 
        print 'New Filtered Variant Total = {}'.format(variantCount)
        uniqueScaffolds = []
        for i in variantList:
            scaffold = i.split('_')[0][1:]
            if scaffold not in uniqueScaffolds:
                uniqueScaffolds.append(scaffold)
        print
        print 'Total Number of Scaffolds to be Placed in Linkage Groups = {}'.format(len(uniqueScaffolds))
        out.write('{} {}\n'.format(sampleCount, variantCount))
        for i in range(len(variantList)):
            outString = variantList[i] + ' ' + crossTypeList[i] + '\t' + ','.join(genotypeList[i]) + '\n'
            out.write(outString)
                 
               
    
    



def print_stats():
    print "\n\n"
    print "-"*60
    print "\nThresholds Used for Filtering:"
    print "The P value Cutoff for ChiSqr for Mendelian Ratios: {}".format(PvalueCutoff)
    print "The Proportion of Variants that had to be Genotyped to Keep the Sample: {}".format(SampleGenotypedCutoff)
    print "The Proportion of Samples that had to be Genotyped to Keep a Variant: {}".format(VariantGenotypedCutoff)
    print "The Name of the Output File is: {}".format(OutfileName)
    

GenotypeList, VariantList, CrossTypeList = read_file(InfileName)  
ExpRatiosDict = build_expected_ratio_dict()
CountingDict= build_counting_dict() ##build a second copy of the dictionary to keep counts in
MarkerCount = len(GenotypeList) # assign global for the number of markers
SampleCount = len(GenotypeList[0]) # assign global for the number of samples
print '\nFound {} samples in file'.format(SampleCount)
if MarkerCount != len(CrossTypeList) or MarkerCount != len(GenotypeList): # double check that number of markers is consistent
    exit('Counts of data lists are not equal! Something is wrong :(')
GenotypeList = remove_bad_genotypes(GenotypeList, CrossTypeList, VariantList, ExpRatiosDict)
SampleList = read_header(HeaderFile)
# print
# print VariantList[0]
# for i in range(len(SampleList)):
#     print SampleList[i] + '\t' + GenotypeList[0][i]
# 
# exit()
GenotypeList, SampleCutList, SampleList = remove_poorly_typed_samples(SampleCount, MarkerCount, SampleGenotypedCutoff, GenotypeList, SampleList)
# print 
# print len(GenotypeList[0])
# print len(SampleList)
# exit()
SampleCount = len(GenotypeList[0])
GenotypeList, CrossTypeList, VariantList = remove_low_covered_variants(GenotypeList, CrossTypeList, VariantList, VariantGenotypedCutoff)
# print
# print
# print len(VariantList)
SelectIndexList = get_indices_for_mendelian_filter(SampleList, SelectSamplesFile)
GenotypeList, CrossTypeList, VariantList = remove_nonMendelian_variants(GenotypeList, CrossTypeList, VariantList, CountingDict, ExpRatiosDict, PvalueCutoff, SelectIndexList)
# print
# print
# print len(VariantList)
output(OutfileName, GenotypeList, CrossTypeList, VariantList)
print_stats()


EndTime = time.time()
Total = EndTime - StartTime
print
print 'COMPLETE.'
print
print 'Script took {} seconds to run'.format(Total)
print

            