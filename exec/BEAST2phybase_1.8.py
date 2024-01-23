#!/usr/bin/env python2
'''Converting a species tree string from the format of *BEAST v.1.8 to the
   format required by PHYBASE'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.09.29.1900"
__status__ = "Working"


#####################
# IMPORT OPERATIONS #
#####################

import sys

# IMPORTANT: This script presently requires DendroPy (= 3.12.0) and does not 
# work correctly with DendroPy (= 3.12.2), because the initial parsing 
# "dendropy.DataSet.get_from_path" is no longer accustomed to the tree 
# annotations of BEAST v.1.8 and, hence, does not recognize all correctly.

import dendropy
import numpy


###############
# DEFINITIONS #
###############

def rmext(inStr, delim="."):
    return(inStr[:inStr.rfind(delim)])


def changing_location_of_metadata(pathAndInfilename):
    ''' Loading a set of trees in *BEAST v.1.8 format and saving them to
        phybase format. Simultaneously converting dmv values to theta value '''
    # Loading *BEAST file into DendroPy
    inData = dendropy.DataSet.get_from_path(pathAndInfilename,
                                            "beast-summary-tree",
                                            extract_comment_metadata=True)
    # Looping over all trees
    for tree in inData.tree_lists[0]:
        nodes = tree.nodes()
        # Looping over every node of a given tree
        for node in nodes:

            # Converting branch lengths in e-notation to decimal format
            if node.edge_length:
                node.edge_length = "%.12f" % (float(node.edge_length))

            # Dealing with dmv-values
            dmvList = node.annotations.get_value("dmv")
            if dmvList:
                # In BEAST v.1.8, dmv = ['0.1','0.2']
                theta = converting_dmv_to_theta(dmvList)
                # Converting theta values in e-notation to decimal format
                theta = "%.12f" % (float(theta))
                # Dropping existing metadata
                node.annotations.drop()
                # Saving new metadata as replacement for dropped metadata
                node.annotations.add_new("theta", theta)
    # Returning tree as string
    outData = inData.as_string("nexus")
    return outData


def converting_dmv_to_theta(inList):
    ''' Converting the dmv values (always 2 per node) to a theta value
        inList: dmv = ['0.1','0.2'] '''
    # Converting strings to floats in list
    tmpList = [float(num) for num in inList]
    # Calculating artithmetic mean
    tmpValue = numpy.mean(tmpList)
    # Multiplying by four because theta = dmv*4*ploidy
    tmpValue = float(tmpValue)*4*2
    # Restrict number of significands to 4 and then convert to string
    outValue = str(round(tmpValue, 4))
    return outValue


def formatting_tree_string(inString):
    ''' Formatting a tree string '''
    # Convert entire string into upper case
    # My Rule: "species tree stuff is always in upper case"
    handle = inString.upper()
    # Remove "taxa" section of nexus file
    handle = handle.replace(handle[handle.find("BEGIN TAXA;"):
                                   handle.find("END;",
                                   handle.find("BEGIN TAXA;"))+4],
                            "")
    handle = handle.replace("[&THETA=", "#")
    handle = handle.replace("[&R]", "")
    # Note: Remove "[&R]" prior to "]"; otherwise "[&R]" is no longer found.
    handle = handle.replace("]", "")
    return handle


########
# MAIN #
########

def main(pathAndInfilename):
    # Conducting main functions
    treeList = changing_location_of_metadata(pathAndInfilename)
    outD = formatting_tree_string(treeList)

    # Writing outF
    outFName = rmext(pathAndInfilename)+".P2C2M.phyb"
    outF = open(outFName, "w")
    outF.write(outD)
    outF.close()


###########
# EXECUTE #
###########

main(sys.argv[1])
