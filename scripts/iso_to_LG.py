#!/usr/bin/env python
##iso_to_LGV2.py
##written 4/22/14 by Groves Dixon
ProgramName = 'iso_to_LG.py'
LastUpdated = '5/30/14'
Updates = """
1. Added the function read_iso2 from
read_iso2.py so that this script can 
automatically append gene names (or any annotation from
and iso2 type file) along with the map locations."""
By = 'Groves Dixon'
VersionNumber = '1.0'

##Assign Global User Help Variables

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This program takes a list of isogroups as input and uses transcriptome-to-genome blast output
to estimate their location on a linkage map relative to linkage mapped RAD markers.
'''

AdditionalProgramInfo = '''
Additional Program Information:

Arguments:
-select 
requires an input file where isogroups are given in first column and in this notation: isogroup=1234
'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
from random import choice
import numpy as np
import scipy.stats
StartTime = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-select', '-s', required = True, metavar = 'selected set of isogroups', dest = 'In', help = 'The name of the input file with the selected set of markers')
parser.add_argument('-align', '-a', required = True, metavar = 'transcriptome to genome blast output', dest = 'Blast', help = "The blastn results from blasting the transcriptome against the genome")
parser.add_argument('-output', '-o', required = True, metavar = 'outputName.out', dest = 'Out', help = 'The desired name for the output file')
parser.add_argument('-map', '-m', required = True, default = 't', metavar = 'The Linkage Map', dest = 'Map', help = "The linkage map")
parser.add_argument('-i', required = False, default = 'none', metavar = 'iso2 file', dest = 'Iso2', help = "Any iso2 file such as iso2gene or iso2go")
parser.add_argument('-tig2group', required = True, dest = 't2g', help = "A table connecting isotigs to isogroups, (probably amil_seq2iso.txt)")
args = parser.parse_args()


#Assign Arguments
InfileName = args.In #good place to get this is from R script selection_analysis2v2.R output 'round3'
OutfileName = args.Out
TranscriptomeCoords = args.Blast
LinkageMap = args.Map
Iso2 = args.Iso2
Tig2Group = args.t2g


def read_input(infileName):
    '''function reads in the set of selected isogorups.
    The input file should have a header and the isogroups
    should be listed in the first column with this style:
    isogroup=12345
    The table can have additional columns that are not read'''
    with open(infileName, 'r') as infile:
        lineNumber = 0
        isogroupList = []
        for line in infile:
            lineNumber += 1
            if lineNumber == 1: continue #skip header
            line = line.strip('\n').split('\t')
            iso = line[0]
            # print iso
            isogroupList.append(line[0])
    return isogroupList

def isotig_to_isogroup(Tig2Group):
    """Builds a dictionary connecting each isotig
    to its correct isogroup"""
    print "\nReading the isotig to isogroup table..."
    t2gDict = {}
    with open(Tig2Group, 'r') as infile:
        for line in infile:
            line = line.strip("\n").split('\t')
            t2gDict[line[0]] = line[1]
    return t2gDict
            
def read_blast_output(TranscriptomeCoords, Tig2groupDict):
    '''This function reads a blast output
    from blasting the transcriptome against a reference genome.
    This is to find out which genomic scaffolds the isogroups
    align to. The output is a dictionary with isogroups as
    keys and scaffolds as values. The blast output is assumed to have
    this format:
    
    isotig	contig	strt	end
    22931	2514691	36379	35410
    8113	1218474	4841	3847
    
    The position of the isogroup on the scaffold 
    is recorded as the middle between the start and end.
    This is so the RAD tag closest to the isogroup center can
    be used to estimate its map position.
    The output is a dictionary connecting isogroups that are in 
    the selected isogroup list with associated scaffolds and positions
    dict = {iso1 : [scaffold, position], iso2 ...}
    '''
    print "Reading the Blast Output to get Gene Coordinates..."
    iso2scaffDict = {}
    previousIsotig = 'none'
    with open(TranscriptomeCoords, 'r') as infile:
        lineNumber = 0
        for line in infile:
            lineNumber += 1
            if lineNumber == 1: continue
            line = line.strip('\n').split('\t')
            isotig = line[0]
            if isotig==previousIsotig: ##blast gives the best hit first, so only record the first incidence of each isotig
                previousIsotig = isotig
                continue
            else:
                previousIsotig = isotig
            isogroup = Tig2groupDict[isotig]
            scaffold = line[1]
            strt = int(line[2])
            end = int(line[3])
            points = [strt, end]
            pos = int(np.mean(points))
            iso2scaffDict[isogroup] = [scaffold, strt, end]##note this will result in replacing any repeats with the lowest entry down
    return iso2scaffDict

def read_linkage_map(LinkageMap):
    """This function reads the linkage map file.
    The file format should be a linearized linkage map (cumulative cM distances for LGs) like this:
    lg\tvariant\tcM
    1\tc1234_5678\t0
    1\tc1234_8903\t5.2135
    
    The output is a dictionary of scaffold map locations with this format:
    dict = {scaffold :[[list of LGs],[list of map positions],[list of scaffold positions]]} 
    the three lists for each scaffold are coindexed so that you can call on them later
    """
    print "Reading Linkage Map..."
    mapDict = {}
    lineNumber = 0
    with open(LinkageMap, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split('\t')  
            if lineNumber == 1:
                header = line
                continue
            lg = line[0]
            radTagList = line[1].split("_")
            scaff = radTagList[0] #grab the scaffold for the RAD locus from the RADtag name
            pos = radTagList[1] #grab the position of the radTag in the scaffold from name
            cM = line[2]
            try:
                mapDict[scaff][0].append(lg) #add LG to its list
                mapDict[scaff][1].append(cM) #add the map position to its list
                mapDict[scaff][2].append(int(pos)) #add the scaffold bp position 
            except KeyError:
                mapDict[scaff] = [[lg],[cM], [int(pos)]] #if the scaffold has not been entered in the dictionary yet, make a new key   
    return mapDict     
        
def read_iso2(Iso2):
    """outputs dictionary pairing isogroups with a
    gene or go term or whatever is in the iso2 file"""
    iso2Dict = {}
    with open(Iso2, 'r') as infile:
        for line in infile:
            line = line.strip('\n').split('\t')
            iso = line[0]
            annot = line[1]
            iso2Dict[iso] = annot
    return iso2Dict

    
def find_scaffs_on_map(IsogroupList, IsoScaffDict, MapDict, Iso2Dict, OutfileName):
    """outputs the start and end locations for each gene along with gene name and scaffold
    """
    print "Locating the Genes on the Map..."
    with open(OutfileName, 'w') as out:
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('isogroup', 'scaffold', 'start', 'end', "annotation", 'lg', 'cM'))
        for i in IsogroupList:
            resultsList = [] # set up list to store the isogroups map data in
            scaffold = IsoScaffDict[i][0] #grab scaffold for this isogroups
            scaffStart = int(IsoScaffDict[i][1]) #grab the scaffold start position
            scaffEnd = int(IsoScaffDict[i][2]) #grab the scaffold start position
            scaffPos = (scaffStart + scaffEnd)/2
            # print '------'
            # print scaffStart
            # print scaffEnd
            # print scaffPos
            # print '\nIsogroup: {}'.format(i)
            try:
                MapDict[scaffold]
            except KeyError:
                resultsList = ['NA','NA']
                continue
            dataLists = MapDict[scaffold]
            ##get the gene if there is one from Iso2Dict
            if Iso2 != "none":
                try:
                    annotation = Iso2Dict[i]
                except KeyError:
                    annotation = 'NA'
            # print 'isogroup position: {}'.format(scaffPos)
            # print MapDict[scaffold]
            scaffPosList = dataLists[2]
            # print scaffPosList
            closestPos = min(scaffPosList, key=lambda x:abs(x-scaffPos)) #find the closest RAD locus position to the isogroup scaffold position
            index = scaffPosList.index(closestPos)
            linkageGroup = dataLists[0][index]
            mapPosition = dataLists[1][index]
            out.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(i, scaffold, scaffStart, scaffEnd, annotation, linkageGroup, mapPosition))
            # print "LG: {}, cM: {}".format(linkageGroup, mapPosition)

Tig2groupDict = isotig_to_isogroup(Tig2Group)
Iso2Dict = read_iso2(Iso2)
MapDict = read_linkage_map(LinkageMap)
IsogroupList = read_input(InfileName)
IsoScaffDict = read_blast_output(TranscriptomeCoords, Tig2groupDict)
IsogroupList = IsoScaffDict.keys() #recode the list




find_scaffs_on_map(IsogroupList, IsoScaffDict, MapDict, Iso2Dict, OutfileName)








