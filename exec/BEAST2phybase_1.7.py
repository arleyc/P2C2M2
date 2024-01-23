#!/usr/bin/env python2
'''Converting a species tree string from the format of *BEAST v.1.7 to the
   format required by PHYBASE'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.10.04.1900"
__status__ = "Working"


#####################
# IMPORT OPERATIONS #
#####################

import sys
import dendropy
import numpy


###############
# DEFINITIONS #
###############

def rmext(inStr, delim="."):
    return(inStr[:inStr.rfind(delim)])


def changing_location_of_metadata(pathAndInfilename):
    ''' Loading a set of trees in *BEAST v.1.7 format and saving them to
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
                node.edge_length = "%.12f" % (node.edge_length)

            # In BEAST v.1.7, dmv is parsed as a string (e.g. dmv='0.15')
            dmv = node.annotations.get_value("dmv")
            # Converting theta values in e-notation to decimal format
            dmv = "%.12f" % (float(dmv))
            theta = float(dmv)*4*2
            # Restrict number of significands to 4 and then convert back 
            # to string
            theta = str(round(float(theta), 4))
            # Dropping existing metadata
            node.annotations.drop()
            # Saving new metadata as replacement for dropped metadata
            node.annotations.add_new("theta", theta)
    # Returning tree as string
    outData = inData.as_string("nexus")
    return outData


def formatting_tree_string(inString):
    ''' Formatting a tree string '''
    # Convert entire string into upper case
    # My Rule: "species tree specs are always in upper case"
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
