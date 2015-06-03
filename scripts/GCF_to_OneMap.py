#!/usr/bin/env python
##GCF_to_OneMap.py
##written 11/6/13 by Groves Dixon
ProgramName = 'GCF_to_OneMap.py'
LastUpdated = '1/24/14'
By = 'Groves Dixon'
VersionNumber = '1.0'

##Assign Global User Help Variables

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script is intended to convert GCF files into input files for OneMap
Genotype data for parents from an input VCF file is used to assign 
Cross types. Then Genotype data for each sample in the VCF file is output
with the appropriate cross type and format for input into OneMap.
'''

AdditionalProgramInfo = '''
Additional Program Information:

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
parser.add_argument('-parentReps', '-preps', default = 'NoParentReps', metavar = 'list of parent replicates', dest = 'ParentReps', help = 'A list of all parent replicates with replicates \
for parent 1 in first column and those for parent 2 in second column')
parser.add_argument('-gcf', required = False, metavar = 'ouputName_forGCF', dest = 'gcfOut', help = 'Option to output the filtered dataset as a GCF also')
parser.add_argument('-pDataCut', required = False, default = 2, type = int, dest = 'pDataCut', help = 'The Number of Parent Replicates that Must be Genotyped to Keep a Given Variant')
parser.add_argument('-pMatchCut', required = False, default = 2, type = int, dest = 'pMatchCut', help = 'The Number of Parent Replicates that Must be Genotyped to Keep a Given Variant')
args = parser.parse_args()


#Assign Arguments
InfileName = args.In
OutfileName = args.Out
ParentReplicateList = args.ParentReps
if args.gcfOut > 0:
    GCFout = args.gcfOut
else:
    GCFout = InfileName[:-4] + '_GCF_to_OneMapOutput.gcf'
ParentDataThreshold = args.pDataCut
ParentMatchThreshold = args.pMatchCut
Verbose = 'N'






def read_GCF(infileName):
    print '\n\n'
    print 'Running {}...'.format(ProgramName)
    variantList = []
    genotypesList = []
    with open(infileName, 'r') as infile:
        for line in infile:
            if line[0] == '#':
                continue
            elif line[0] == '*':
                sampleList = line.strip('\n').split('\t')[1:]
                continue
            else:
                line = line.strip("\n").split("\t")
                variantList.append(line[0])
                genotypesList.append(line[1:])
    print '\nFound {} samples in the file'.format(len(sampleList))
    print '\nFound {} variants in the file'.format(len(variantList))
    return variantList, genotypesList, sampleList

def build_genotype_dictionary(variantList, genotypesList, sampleList):
    '''Use the data from the input file to build a dictionary keyed with the samples 
    and the list of genotypes as the values indexed the same as the variant list'''
    genotypeDict = {}
    for i in sampleList:
        genotypeDict[i] = []
    for i in range(len(variantList)):
        variant = variantList[i]
        genotypes  = genotypesList[i]
        for x in range(len(genotypes)):
            genotypeDict[sampleList[x]].append(genotypes[x])
    print '\nThere are genotypes for {} variants (should be same number as there are variants)'.format(len(genotypesList))
    return genotypeDict
        
def parent_rep_lists(parentRepList):
    '''Function to assign the two lists of parent replicates from the given parent replicate list
    '''
    p1Reps = []; p2Reps = []
    with open(parentRepList, 'r') as infile:
        for line in infile:
            line = line.strip('\n').split('\t')
            if line[0] != '':
                p1Reps.append(line[0])
            if line[1] != '':
                p2Reps.append(line[1])
    print
    print 'Parent1 replicates:'
    print p1Reps
    print 'Parent2 replicates:'
    print p2Reps
    print
    return p1Reps, p2Reps

def get_offspring_list(sampleList, p1Reps, p2Reps):
    '''pulls out the list of just the offspring'''
    offspringList = []
    for i in sampleList:
        if i not in p1Reps:
            if i not in p2Reps:
                offspringList.append(i)
    print '\nFound {} Offsrping in the File'.format(len(offspringList))
    return offspringList

def test_heterozygosity(genotype):
    '''Function to test if a genotype is heterozygous
    Returns either True or False
    Used to check parents for assigning cross types in 
    convert_genotypes()
    '''
    if genotype[0] == genotype[1]:
        answer = 'FALSE'
    else:
        answer = 'TRUE'
    return answer



def get_modal_geno(genotypeDict, samples, genoList):
    '''function to grab the modal genotype for parent replicates (ties are OK)
    genoList is a list of genotypes indexed the same as the list
    of samples
    '''
    matchScoreList = []
    for x in genoList:
        matchScoreList.append(genoList.count(x))
    maxIndex = matchScoreList.index(max(matchScoreList))
    mode = genoList[maxIndex]
    return mode

def get_anchor_geno(genotypeDict, pReps, variantNumber):
    '''function to grab the modal genotype for parent replicates
    For use in filter_variants_by_parent() below.
    Arguements:
    1. the genotype dictionary
    2. the list of replicates for 1 parent
    3. the variant number (the index for the variant for which we are getting an anchor genotype)
    '''
    genotypeList = []
    allMissing = 'TRUE'
    for i in pReps:
        pGeno = genotypeDict[i][variantNumber]
        genotypeList.append(pGeno)
    for i in genotypeList:
        if '-' not in i:
            allMissing = 'FALSE'
    modeGeno = get_modal_geno(genotypeDict, pReps, genotypeList)
    while '-' in modeGeno:
        if allMissing == 'TRUE':
            anchorGeno = modeGeno
            break
        else:
            newGenoList = []
            for z in genotypeList:
                if '-' not in z:
                    newGenoList.append(z)
            modeGeno = get_modal_geno(genotypeDict, pReps, newGenoList)
    anchorGenotype = modeGeno
    return anchorGenotype
    

def filter_variants_by_parent(pReps, genotypeDict, variantList, sampleList, dataThreshold, matchThreshold):
    '''This function filters variants based on the availability and agreement of 
    data from the parent replicates. It also  develops concensus parent data
    by pulling the 'best' genotype for each variant based on the parent replicates.
    
    '''
    keepList = []
    trashList = []
    startingTotal = len(variantList)
    trashedMissing = 0
    trashedDisagree = 0
    kept = 0
    concensusGenoList = [] #list to keep track of the best genotype for each variant for this set of parent replicates
    for i in range(len(variantList)):
        parentGenoList = []
        variant = variantList[i]
        matchScore = 0 #to keep track of how many matching genotypes there are
        dataScore = 0 #to keep track of how many replicates have data
        for p in pReps:
            parentGenoList.append(genotypeDict[p][i])
        anchorGeno = get_anchor_geno(genotypeDict, pReps, i)
        concensusGenoList.append(anchorGeno)
        for x in parentGenoList:
            if '-' not in x:
                dataScore += 1
            if x == anchorGeno:
                matchScore += 1
        if dataScore < dataThreshold:
            trashedMissing += 1
            trashList.append(i)
        elif matchScore < matchThreshold:
            trashedDisagree += 1
            trashList.append(i)
        else:
            kept += 1
            keepList.append(i)
    trashTotal = trashedMissing + trashedDisagree
    print "\nFiltering Variants Based on Parent Replicates with {}...".format(pReps[0])
    print "\n\tStarting with {} Variants".format(startingTotal)
    print '\n\t{} Variants Did Not Have Enough Replicates Genotyped'.format(trashedMissing)
    print '\n\t{} More Variants Had Too Much Disagreement Between Replicates'.format(trashedDisagree)
    print '\n\t{} Total Variants will be Removed'.format(trashTotal)
    print '\n\t{} Variants Passed Filter'.format(kept)
    return keepList, concensusGenoList
        
def remove_variants(variantList, genotypeDict, sampleList, keepList, verbose):
    '''functiont to iterate through the genotype data and only keep
    those in the keepList'''
    toRemove = len(variantList) - len(keepList)
    if verbose == 'Y':
        print '\n\tRemoving {} Variants from the Data...'.format(toRemove)
    newVariantList = []
    newGenotypeDict = {}
    removed = 0
    kept = 0
    for x in sampleList:
        newGenotypeDict[x] = []
    for i in range(len(VariantList)):
        if i not in keepList:
            removed += 1
            continue
        else:
            kept += 1
            variant = variantList[i]
            newVariantList.append(variant) ##updates the new variant list
            for z in sampleList:
                genotype = genotypeDict[z][i] ##get the genotype for each sample for this variant
                newGenotypeDict[z].append(genotype) ##update the new genotype dict
    newLength = len(newVariantList)
    newLength2 = len(newGenotypeDict[sampleList[1]])
    if newLength != newLength2:
        exit("PROBLEM WITH REMOVING VARIANTS!")
    if verbose == 'Y':
        print '\n\t{} variants have been removed for each sample'.format(len(variantList) - newLength)
        print '\n\t{} Variants are Left in the Dataset'.format(len(newVariantList))
    return newVariantList, newGenotypeDict
        
def filter_uninformative(variantList, genotypeDict, p1, p2):
    '''function to find the uninformative variants
    (those that are homozygous in both parents)'''
    print '\nFinding Uninformative Variants...'
    keepList = []
    uninformative = 0
    informative = 0
    for i in range(len(variantList)):
        p1Geno = genotypeDict[p1][i].split('/')
        p2Geno = genotypeDict[p2][i].split('/')
        if test_heterozygosity(p1Geno) == 'FALSE' and test_heterozygosity(p2Geno) == 'FALSE':
            uninformative += 1
            continue
        else:
            informative += 1
            keepList.append(i) ##put variant indices in the keep list only if at least one parent is a heterozygote
    print '\n\t{} Uninformative Variants were found'.format(uninformative)
    print '\n\t{} Were Informative'.format(informative)
    return keepList
        

def reassign_nonparental_alleles(genotypeDict, sampleList, variantList, offspringList, p1, p2):
    '''This function reassigns any alleles from the offspring 
    '''
    print '\nFinding Samples with Nonparental Alleles...'
    p1GenoList = genotypeDict[p1]
    p2GenoList = genotypeDict[p2]
    badCalls = 0
    goodCalls = 0
    missing = 0
    for x in range(len(variantList)):
        variant = variantList[x]
        p1Geno = p1GenoList[x].split('/')
        p2Geno = p2GenoList[x].split('/')
        pAlleles = p1Geno + p2Geno
        for i in offspringList:
            sampleGeno = genotypeDict[i][x].split('/')
            if '-' in sampleGeno:
                missing += 1
                continue
            allParental = 'TRUE'
            for allele in sampleGeno:
                if allele not in pAlleles:
                    allParental = 'FALSE'
            if allParental == "FALSE":
                badCalls += 1
                genotypeDict[i][x] = '-/-'
            else:
                goodCalls += 1
    print '\n\t{} genotype calls were changed to missing data because of disagreement with parents'.format(badCalls)
    print '\n\t{} calls were not changed'.format(goodCalls)
    print '\n\t{} calls were already missing'.format(missing)
    print '\n\tProduct of Number of Variants and Number of Offspring = {}'.format(len(variantList) * len(offspringList))
    return genotypeDict  


def output_GCF(GenotypeDict, SampleList, VariantList, GCFout, InfileName):
    print '\nOutputting Filtered Data as a GCF file...'
    with open(GCFout, 'w') as out:
        out.write('#data from {} as simple genotype table'.format(InfileName))
        out.write('\n*\t' + '\t'.join(SampleList))
        for i in range(len(VariantList)):
            genotypeList = []
            variant = VariantList[i]
            for x in SampleList:
                genotype = GenotypeDict[x][i]
                genotypeList.append(genotype)
                variantString = '\t'.join(genotypeList)
            out.write('\n' + variant + '\t' + variantString)
    print '\n\tGCF file saved as {}'.format(GCFout)


def do_conversions(variant, convertedVariantList, offspringList, converter, genotypeDict, convertedDict, index):
    '''function to perform conversions of sample alleles using a crossType specific converter
    This is ussed in the functions that assign the variatns to each crosstype and convert their 
    data to oneMap format'''
    convertedVariantList.append(variant)
    for sample in offspringList:
        genotype = genotypeDict[sample][index].split('/')
        a1 = genotype[0]
        a2 = genotype[1]
        if a1 == '-' and a2 == '-':
            convertedDict[sample].append('-')
            continue
        convertedGeno = [converter[a1]]
        if converter[a2] not in convertedGeno:
            convertedGeno.append(converter[a2]) ##this deals with possibility of a homozygous offspring, although technically this should not happen in some crosstypes
        convertedGeno.sort()
        convertedGeno2 = ''.join(convertedGeno)
        convertedDict[sample].append(convertedGeno2)
        

def gather_parent_genotype_data(genotypeDict, index):
    p1Geno = genotypeDict[p1][index].split('/')
    p2Geno = genotypeDict[p2][index].split('/')
    p1a1 = p1Geno[0]; p1a2 = p1Geno[1]; p2a1 = p2Geno[0]; p2a2 = p2Geno[1]
    p1Het = test_heterozygosity(p1Geno)
    p2Het = test_heterozygosity(p2Geno)
    parentAlleles = [p1a1, p1a2, p2a1, p2a2]
    return p1Geno, p2Geno, p1a1, p1a2, p2a1, p2a2, p1Het, p2Het, parentAlleles

def print_cross_type_results(count, condition, variant_list, crossType):
    print '\n\tFound {} Loci with Cross Type {}'.format(count, crossType)
    if condition == 'Y':
        print '\nVariants with {} CrossType:'.format(crossType)
        for var in variant_list:
            print '\t'+var


def get_A1_crosstypes(genotypeDict, sampleList, variantList, p1, p2, verbose):
    '''Rather than assign all cross types at once, we 
    will do it piece-wise. This function will deal only with the A.1 crosstypes.
    A.1 crosstype: ab X cd
    There are probably very few if any of these since most loci are biallelic.
    A.1 cross types have 4 parental alleles, so both parents are heterozygotes,
    and do not share any alleles.
    Two global variables are generated:
    ConvertDict and CrossTypeList
    which will be added to as variants identified as other cross types
    by the following functions'''
    convertedDict = {}
    convertedVariantList = []
    A1List = [] #keep the variants with A1 cross Type
    crossTypeList = []
    keepList = []
    assignmentCount = 0
    for i in sampleList:
        convertedDict[i] = []
    for i in range(len(variantList)):
        variant = variantList[i]
        p1Geno, p2Geno, p1a1, p1a2, p2a1, p2a2, p1Het, p2Het, parentAlleles = gather_parent_genotype_data(genotypeDict, i)
        if p1Het == 'FALSE':
            keepList.append(i) 
            continue
        if p2Het == 'FALSE': #skip any locus if a parent isn't heterozygous
            keepList.append(i)
            continue
        unique = []
        for x in parentAlleles:
            if x not in unique:
                unique.append(x)
        totalUnique = len(unique)
        if totalUnique == 4:
            assignmentCount += 1
            crossType = 'A.1'
            A1List.append(variantList[i])
            crossTypeList.append(crossType)
            converter = {p1a1: 'a', p1a2: 'b', p2a1: 'c', p2a2: 'd'}
            do_conversions(variant, convertedVariantList, OffspringList, converter, GenotypeDict, convertedDict, i)
        else:
            keepList.append(i)
    print_cross_type_results(assignmentCount, verbose, A1List, 'A.1')
    return convertedVariantList, convertedDict, crossTypeList, keepList, assignmentCount

def get_A2_crosstypes(convertedVariantList, genotypeDict, convertedDict, crossTypeList, sampleList, variantList, p1, p2, verbose):
    '''Function to deal with the A.2 cross types.
    A.2 crosstype: ab X ac
    Again expect very few if any of these.
    Both parents are heterozygotes and share a single allele.
    The shared allele is labeled as 'a', the unshared allele in
    parent1 is labeled 'b', the unshared in parent2 as 'c'.
    '''
    A2List = []
    keepList = []
    assignmentCount = 0
    for i in range(len(variantList)):
        variant = variantList[i]
        p1Geno, p2Geno, p1a1, p1a2, p2a1, p2a2, p1Het, p2Het, parentAlleles = gather_parent_genotype_data(genotypeDict, i)
        if p1Het == 'FALSE':
            keepList.append(i) 
            continue
        if p2Het == 'FALSE': #skip any locus if a parent isn't heterozygous
            keepList.append(i)
            continue
        unique = []
        for x in parentAlleles:
            if x not in unique:
                unique.append(x)
        totalUnique = len(unique)
        if totalUnique == 3:
            assignmentCount += 1
            crossType = 'A.2'
            A2List.append(variantList[i])
            crossTypeList.append(crossType)
            if p1a1 in p2Geno:
                p2Geno.remove(p1a1)
                converter = {p1a1: 'a', p1a2: 'b', p2Geno[0]: 'c'} #if p1a1 is shared, make it 'a' 
            else:
                if p1a2 not in p2Geno:
                    exit('ERROR!') #this should not happen unless there's a mistake
                p2Geno.remove(p1a2)
                converter = {p1a1: 'b', p1a2: 'a', p2Geno[0]: 'c'} #if p1a2 is shared, make it 'a'
            do_conversions(variant, convertedVariantList, OffspringList, converter, GenotypeDict, convertedDict, i)
        else:
            keepList.append(i) #the variant is not crossType A.2 so retain it for further iterations
    print_cross_type_results(assignmentCount, verbose, A2List, 'A.2')
    return keepList, assignmentCount
        
def get_B3_crosstypes(convertedVariantList, genotypeDict, convertedDict, crossTypeList, sampleList, variantList, p1, p2, verbose):
    '''function to get the B3.7 cross type variants. 
    B3.7 crosstype: ab X ab
    Here both parents are heterozygotes and have the same genotype'''
    B3List = []
    keepList = []
    assignmentCount = 0
    for i in range(len(variantList)):
        variant = variantList[i]
        p1Geno, p2Geno, p1a1, p1a2, p2a1, p2a2, p1Het, p2Het, parentAlleles = gather_parent_genotype_data(genotypeDict, i)
        if p1Het == 'FALSE':
            keepList.append(i) 
            continue
        if p2Het == 'FALSE': #skip any locus if a parent isn't heterozygous
            keepList.append(i)
            continue
        unique = []
        for x in parentAlleles:
            if x not in unique:
                unique.append(x)
        totalUnique = len(unique)
        if totalUnique == 2:
            assignmentCount += 1
            crossType = 'B3.7'
            B3List.append(variantList[i])
            crossTypeList.append(crossType)
            converter = {p1a1: 'a', p1a2: 'b'}
            if p1a1 not in p2Geno or p1a2 not in p2Geno:
                exit('ERROR!') #this should not happen
            
            # print
            # print p1Geno
            # print p2Geno
            
            do_conversions(variant, convertedVariantList, OffspringList, converter, GenotypeDict, convertedDict, i)
        else:
            keepList.append(i)
    print_cross_type_results(assignmentCount, verbose, B3List, 'B3.7')
    return keepList, assignmentCount
    




def get_D_crosstypes1(convertedVariantList, genotypeDict, convertedDict, crossTypeList, sampleList, variantList, p1, p2, verbose):
    '''function to gather the variants with D1.9 and D2.14 crosstypes
    D1.9 crossType: ab X cc
    D2.14 crossType: cc X ab
    These crosstypes are reciprocally equivalent and equally informative for
    linkage anaysis. So they are dealt with simultaneously in this funciton.
    1 parent is homozygous and the other is heterozygous and the parents do
    not share any alleles. The unshared homozygous allele is labeled "c", 
    the alleles of the heterozygous parent are labled "a" and "b".'''
    D9List = []
    D14List = [] #keep track of them individually
    DList = [] #keep track of both
    keepList = []
    D9assignmentCount = 0
    D14assignmentCount = 0 #keep count individually
    DassignmentCount = 0 #keep count of both
    for i in range(len(variantList)):
        variant = variantList[i]
        p1Geno, p2Geno, p1a1, p1a2, p2a1, p2a2, p1Het, p2Het, parentAlleles = gather_parent_genotype_data(genotypeDict, i)
        if p1Het == 'TRUE' and p2Het == 'FALSE':
            if p2a1 in p1Geno or p2a2 in p1Geno: #this designates D1.10 crosstype so skip it here
                keepList.append(i)
                continue
            else:
                if p2a1 != p2a2:
                    exit('ERROR!') #We should be dealing with D1.9, so this should not happen
                crossType = 'D1.9'
                D9assignmentCount += 1
                DassignmentCount += 1
                DList.append(variantList[i])
                D9List.append(variantList[i])
                crossTypeList.append(crossType)
                converter = {p1a1: 'a', p1a2: 'b', p2a1: 'c'}
                do_conversions(variant, convertedVariantList, OffspringList, converter, GenotypeDict, convertedDict, i)
        elif p1Het == 'FALSE' and p2Het == 'TRUE':
            if p1a1 in p2Geno or p1a2 in p2Geno: #this would be D2.15, so skip here
                keepList.append(i)
                continue
            else:
                if p1a1 != p1a2:
                    exit('ERROR!') #shouldn't happen
                crossType = 'D2.14'
                D14assignmentCount += 1
                D14List.append(variantList[i])
                DassignmentCount += 1
                DList.append(variantList[i])
                crossTypeList.append(crossType)
                converter = {p2a1: 'a', p2a2: 'b', p1a1: 'c'}
                do_conversions(variant, convertedVariantList, OffspringList, converter, GenotypeDict, convertedDict, i)
        else:
            keepList.append(i)
            print "this shouldn't happen..."
            print "look into variant {}".format(variantList[i])
            exit()
    print_cross_type_results(DassignmentCount, verbose, DList, 'D1.9/D2.14')
    return keepList, D9assignmentCount, D14assignmentCount
    
    
    
def get_D_crosstypes2(convertedVariantList, genotypeDict, convertedDict, crossTypeList, sampleList, variantList, p1, p2, verbose):
    D10List = []
    D15List = []
    DList = []
    keepList = []
    D10assignmentCount = 0
    D15assignmentCount = 0
    DassignmentCount = 0
    for i in range(len(variantList)):
        variant = variantList[i]
        p1Geno, p2Geno, p1a1, p1a2, p2a1, p2a2, p1Het, p2Het, parentAlleles = gather_parent_genotype_data(genotypeDict, i)
        if p1Het == 'TRUE' and p2Het == 'FALSE':
            if p2a1 in p1Geno and p2a2 in p1Geno: #this designates our D1.10 dudes
                D10List.append(variantList[i])
                DList.append(variantList[i])
                D10assignmentCount += 1
                DassignmentCount += 1
                crossType = 'D1.10'
                crossTypeList.append(crossType)
                if p1a1 in p2Geno:
                    converter = {p1a1: 'a', p1a2: 'b'}
                elif p1a2 in p2Geno:
                    converter = {p1a2: 'a', p1a1: 'b'}
                else:
                    exit('ERROR!') #shouldn't happen
                do_conversions(variant, convertedVariantList, OffspringList, converter, GenotypeDict, convertedDict, i)
        elif p1Het == 'FALSE' and p2Het == 'TRUE':
            if p1a1 in p2Geno and p1a2 in p2Geno:
                D15List.append(variantList[i])
                DList.append(variantList[i])
                D15assignmentCount += 1
                DassignmentCount += 1
                crossType = 'D2.15'
                crossTypeList.append(crossType)
                if p2a1 in p1Geno:
                    converter = {p2a1: 'a', p2a2: 'b'}
                elif p2a2 in p1Geno:
                    converter = {p2a2: 'a', p2a1: 'b'}
                else:
                    exit('ERROR!') #shouldn't happen
                do_conversions(variant, convertedVariantList, OffspringList, converter, GenotypeDict, convertedDict, i)
        else:
            keepList.append(i)
            print 'This shoud not happen'
            print "look into variant {}".format(variantList[i])
            exit()
    print_cross_type_results(DassignmentCount, verbose, DList, 'D1.10/D2.15')
    return keepList, D10assignmentCount, D15assignmentCount
    
                
        
    
def check_leftover(notDList2):
    print "\nFinished Assigning Cross Types"
    leftOver = len(notDList2)
    print "\n{} Genotypes Failed to be Assigned".format(leftOver)
    if leftOver != 0:
        print "\nHmmm... There should not be any left over"
        print "Please look into these Variants:"
        for i in notDList2:
            print i
    else:
        print "\nNice. All Variants Were Assigned"
def report_results(A1count, A2count, B3count, D9count, D10count, D14count, D15count):
    print '\n Filtering Parameters Used:'
    print 'Parent Replicates Genotyped Threshold: {}'.format(ParentDataThreshold)
    print 'Parent Replicates in Agreement Threshold: {}'.format(ParentMatchThreshold)
    print "\nFinal Results:"
    print "A.1 Crosstype Variants: {}".format(A1count)
    print "A.2 Crosstype Variants: {}".format(A2count)
    print "B3.7 Crosstype Variants: {}".format(B3count)
    print "D1.9 Crosstype Variants: {}".format(D9count)
    print "D1.10 Crosstype Variants: {}".format(D10count)
    print "D2.14 Crosstype Variants: {}".format(D14count)
    print "D2.15 Crosstype Variants: {}".format(D15count)
    total = A1count+A2count+B3count+D9count+D10count+D14count+D15count
    print "Total Variants Converted: {}".format(total)
    print "\nResults in OneMap Format Output as {}".format(OutfileName)
    return total



def output(outfileName, convertedVariantList, convertedDict, crossTypeList, offspringList, total):
    with open(outfileName, 'w') as out:
        header = "{} {}".format(len(offspringList), total) #the header of the OneMap input is #samples #loci
        out.write(header)
        for i in range(total):
            variant = convertedVariantList[i]
            crossType = crossTypeList[i]
            genotypeList = []
            for baby in offspringList:
                geno = convertedDict[baby][i]
                genotypeList.append(geno)
            genotypeString = ','.join(genotypeList)
            outString = "\n*{} {}\t{}".format(variant, crossType, genotypeString)
            out.write(outString)
            
            
    
    


VariantList, GenotypesList, SampleList = read_GCF(InfileName)
GenotypeDict = build_genotype_dictionary(VariantList, GenotypesList, SampleList)
P1reps, P2reps = parent_rep_lists(ParentReplicateList)
OffspringList = get_offspring_list(SampleList, P1reps, P2reps)
KeepList, P1concensus = filter_variants_by_parent(P1reps, GenotypeDict, VariantList, SampleList, ParentDataThreshold, ParentMatchThreshold)
p1 = 'PARENT1' #set global variable for concensus for parent 1
GenotypeDict[p1] = P1concensus ##update the genotypeDict with the best genotype data for the first set of parent replicates
SampleList.append(p1)
VariantList, GenotypeDict = remove_variants(VariantList, GenotypeDict, SampleList, KeepList, Verbose)
KeepList2, P2concensus = filter_variants_by_parent(P2reps, GenotypeDict, VariantList, SampleList, ParentDataThreshold, ParentMatchThreshold)
p2 = 'PARENT2'
GenotypeDict[p2] = P2concensus
SampleList.append(p2)
VariantList, GenotypeDict = remove_variants(VariantList, GenotypeDict, SampleList, KeepList2, Verbose)
KeepList3 = filter_uninformative(VariantList, GenotypeDict, p1, p2)
VariantList, GenotypeDict = remove_variants(VariantList, GenotypeDict, SampleList, KeepList3, Verbose)
reassign_nonparental_alleles(GenotypeDict, SampleList, VariantList, OffspringList, p1, p2)
print "\nFinished Filtering Variants"
print
print "-"*60
print "\nBegin Gathering the Variants that Match Each CrossType..."
print "\nStarting with {} Filter-Passing Variants".format(len(VariantList))
output_GCF(GenotypeDict, SampleList, VariantList, GCFout, InfileName)
ConvertedVariantList, ConvertedDict, CrossTypeList, notA1List, A1count = get_A1_crosstypes(GenotypeDict, SampleList, VariantList, p1, p2, Verbose)
VariantList, GenotypeDict = remove_variants(VariantList, GenotypeDict, SampleList, notA1List, Verbose)
notA2List, A2count = get_A2_crosstypes(ConvertedVariantList, GenotypeDict, ConvertedDict, CrossTypeList, SampleList, VariantList, p1, p2, Verbose)
VariantList, GenotypeDict = remove_variants(VariantList, GenotypeDict, SampleList, notA2List, Verbose)
notB3List, B3count = get_B3_crosstypes(ConvertedVariantList, GenotypeDict, ConvertedDict, CrossTypeList, SampleList, VariantList, p1, p2, Verbose)
VariantList, GenotypeDict = remove_variants(VariantList, GenotypeDict, SampleList, notB3List, Verbose)
notDList1, D9count, D14count = get_D_crosstypes1(ConvertedVariantList, GenotypeDict, ConvertedDict, CrossTypeList, SampleList, VariantList, p1, p2, Verbose)
VariantList, GenotypeDict = remove_variants(VariantList, GenotypeDict, SampleList, notDList1, Verbose)
notDList2, D10count, D15count = get_D_crosstypes2(ConvertedVariantList, GenotypeDict, ConvertedDict, CrossTypeList, SampleList, VariantList, p1, p2, Verbose)
check_leftover(notDList2)
Total = report_results(A1count, A2count, B3count, D9count, D10count, D14count, D15count)
# print len(CrossTypeList)
# print len(ConvertedDict[OffspringList[1]])
output(OutfileName, ConvertedVariantList, ConvertedDict, CrossTypeList, OffspringList, Total)

