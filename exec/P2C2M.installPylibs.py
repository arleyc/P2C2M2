#!/usr/bin/env python2
'''Automatic installation script for Python libraries of P2C2M'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.10.13.2000"
__status__ = "Working"

###########
# IMPORTS #
###########

import sys
from pkg_resources import WorkingSet, DistributionNotFound

###############
# DEFINITIONS #
###############

def main(pkgList):
    working_set = WorkingSet()
    for pkgName, pkgVersion in pkgList.iteritems():
        try:
            depends = working_set.require(pkgName)
        except DistributionNotFound:
            import os
            import urllib2
            print "\n -- Library " + pkgName + " needs to be installed --\n"
            # Prompt for user decision if a package requires installation
            allow = raw_input(" May I install the above library? ([y]/n): ")
            if allow.upper() == "N"  or allow.upper() == "NO":
                sys.exit("\n -- Please install package " + pkgName +
                         " manually. Aborting. --\n")
            else:
                try:
                    response = urllib2.urlopen('http://www.google.com', timeout=20)
                    os.system("easy_install-2.7 --user " + pkgName + "==" + pkgVersion)
                    # Important: After installation via easy_install, I must 
                    # restart the script for certain new modules (i.e. those 
                    # installed via egg,like dendropy) to be properly loaded
                    os.execv(__file__, sys.argv)
                    # TFLs, while nicer, fails because easy_install must be
                    # run as superuser
                    #   from setuptools.command.easy_install import main as install
                    #   install([pkg])
                except urllib2.URLError as err:
                    sys.exit("\n -- Please connect to the internet. Aborting. --\n")
        # Eventually, import the module
        exec("import " + pkgName)
        ## Alternative: modules = map(__import__, pkgList)

    ## Completion notice
    print "\n -- P2C2M: Installation of Python libraries complete --\n"

###########
# EXECUTE #
###########

main({"dendropy": "3.12.0", "numpy": "1.9.0"})
