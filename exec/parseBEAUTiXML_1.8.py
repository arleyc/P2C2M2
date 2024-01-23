#!/usr/bin/env python2
# '''Parsing XML files that were generated by BEAUTi v.1.8'''
# __author__ = "Michael Gruenstaeudl, PhD"
# __copyright__ = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
# __email__ = "gruenstaeudl.1@osu.edu"
# __version__ = "2015.03.31.1400"
# __status__ = "Working"

# I/P:
# InFn = filename of the xml input file
#
# O/P:
# SAmtrx = species allele matrix
# LPmtrx = locus ploidy matrix
# GPmtrx = gene partition matrix
# ALNdict = alignment dictionary

#####################
# IMPORT OPERATIONS #
#####################

import xml.etree.ElementTree as etree
from itertools import izip

############
# RAW CODE #
############

# TFL corrects the filename to be in ascii format
# inFn = inFn.encode('ascii')

# "inFn" is the filename of the xml input file
# etree.parse() parses a file from the harddisk
rootElem = etree.parse(inFn).getroot()

#########################################
# 1. Extract species-allele association #                               # Extract the association of "species" to "alleles" and save to SAmtrx         # 
#########################################
SAmtrx = []
if bool(rootElem.findall(".//sp")):                                     # Test if list generated by findall empty
    sp = [sp for sp in rootElem.iter("sp") if sp.attrib.keys()[0]=="id"]
    for s in sp:
        taxa = [taxon for taxon in s.iter("taxon")]
        for t in taxa:
            SAmtrx.append([s.get("id"), t.get("idref")])
else:
    raise ValueError('XML file malformed')

#######################################
# 3. Extract locus-ploidy association #                                 # Extract the association of "gene names" to "ploidy" and save to LPmtrx
#######################################
LPmtrx = []
if bool(rootElem.findall(".//gtree")):
    for entry in rootElem.iter("gtree"):
        ploidy = entry.get("ploidy")
        gene = [gene.get("idref") for gene in entry.iter("treeModel")]
        # Is the gene[0] correct here?
        gene = gene[0].rstrip("treeModel").rstrip(".")                  # .rstrip("treeModel").rstrip(".") preferential over .rstrip(".treeModel"), because R interprets the latter as "any treeModel"
        LPmtrx.append([gene, ploidy])
else:
    if bool(rootElem.findall(".//geneTrees")):
        for entry in rootElem.iter("geneTrees"):                        # There is only one entry with tag "geneTrees"
            genes = [gene.get("idref") for gene in entry.iter("treeModel")]
            gene_names = [gene.rstrip("treeModel").rstrip(".") for gene in genes]
            ploidy = 2                                                  # Generate diploid default ploidy levels
            for name in gene_names:
                LPmtrx.append([name, ploidy])
    else:
        raise ValueError('XML file malformed')

##################################
# 4. Extract gene-partition info #                                      # Extract gene-partition info from allDict and save to GPmtrx
##################################
GPmtrx = []
if bool(rootElem.findall(".//treeLikelihood")):
    tL = [tL for tL in rootElem.iter("treeLikelihood") if tL.attrib.keys()[0]=="id"]
    for l in tL:
        gene = [gene.get("idref") for gene in l.iter("treeModel")]
        gene = gene[0].rstrip("treeModel").rstrip(".")
        part = [part.get("idref") for part in l.iter("patterns")]
        part = part[0].rstrip("patterns").rstrip(".")
        GPmtrx.append([gene, part])
else:
    raise ValueError('XML file malformed')

#########################
# 5. Extract alignments #
#########################
# Extract alignments and save to ALNdict
# Note:  simple (or nested) Python-lists are being converted to simple
#        (or nested) unnamed R-lists by the rPython-package;
# Note:  Python-dictionaries are being converted to named R-lists by the
#        rPython-package, whereby the new R-list name is the former
#        Python-dictionary key
ALNdict = {}
# Some alignment elements only have a single attrib.; they can be cleaned out by this fact
alignms = [alignm for alignm in rootElem.iter("alignment") if len(alignm.attrib.keys()) > 1]
patterns = [pat for pat in rootElem.iter("patterns") if len(pat.attrib.keys()) >= 3]
for (a, p) in izip(alignms, patterns):
    taxDict = {}
    for taxon in a.iter("taxon"):
        taxDict[taxon.get("idref")] = list(taxon.tail.replace("\n","").replace("\t",""))
    patternDict = {}
    for alignm in p.iter("alignment"):
        patternDict[alignm.get("idref")] = p.get("id").rstrip("patterns").rstrip(".")
    ALNdict[patternDict[a.get("id")]] = taxDict

#for alignm in allDict["alignment"]:
#    tmpDict = {}
#    for seq in alignm.findall("sequence"):
#        tmp = list(seq[0].tail.replace("\t", "").replace("\n", ""))
#        tmpDict[seq[0].attrib["idref"]] = tmp
#
#    geneNum = int(alignm.attrib["id"].lstrip("alignment"))-1
#
#    geneName = allDict["treeLikelihood"][geneNum]
#    geneName = geneName.findall("treeModel")[0].attrib["idref"]
#    geneName = geneName.rstrip("treeModel").rstrip(".")
#
#    ALNdict[geneName] = tmpDict
