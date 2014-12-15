#!/usr/bin/python

import commands
import fileinput
import sys
import os

TESTING        = 0
QUEUE          = '1nd'
EVENTS_PER_JOB = 10000
TOT_EVENTS     = 1000000
CMSSW_FOLDER   = '/afs/cern.ch/user/r/rgerosa/work/TP_ANALYSIS/DELPHES_ANALYSIS/CMSSW_6_2_0_SLHC20_patch1/src/'
DELPHES_FOLDER = CMSSW_FOLDER + '/Delphes'
GEN_MINBIAS    = DELPHES_FOLDER + '/genMinBias_14TeV'
TRANS_MINBIAS  = DELPHES_FOLDER + '/hepmc2pileup'
TUNE_MINIBIAS  = 5
CMSStagefolder = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/TP_ANALYSIS/PU_FILES_14TEV/'
EOSfolder      = '/eos/cms' + CMSStagefolder

rootfolder = os.getcwd ()
    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def runCommand (command, printIt = 0, doIt = 1) :
    if printIt : print ('> ' + command)
    if doIt : 
        commandOutput = commands.getstatusoutput (command)
        if printIt : print commandOutput[1]
        return commandOutput[0]
    else :    print ('    jobs not submitted')
    return 1
    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def replaceAll (file,searchExp,replaceExp) :
    for line in fileinput.input (file, inplace = 1) :
        if searchExp in line:
            line = line.replace (searchExp, replaceExp)
        sys.stdout.write (line)


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def prepareJob (tag) :
    filename = 'run_' + tag + '.job'
    f = open (filename, 'w')
    f.write ('cd /afs/cern.ch/user/r/rgerosa/work/TP_ANALYSIS/DELPHES_ANALYSIS/CMSSW_6_2_0_SLHC20_patch1/src/\n')
    f.write ('eval `scram run -sh`\n')
    f.write ('cd -\n')
    outfilename = tag + '.tempo'
    finalfilename = tag + '.mb'
    f.write (GEN_MINBIAS + ' ' + str (EVENTS_PER_JOB) + ' '  + outfilename + ' ' + str(TUNE_MINIBIAS) + '\n')
    f.write (TRANS_MINBIAS + ' ' + finalfilename +  ' '  + outfilename + '\n')
    f.write ('cmsStage ' + finalfilename + ' ' + CMSStagefolder + '\n')
    f.close ()
    return filename


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == "__main__":

    eoscmd = '/afs/cern.ch/project/eos/installation/cms/bin/eos.select' ;

    if (TESTING == 1) :     
        print '  --- TESTNG, NO submissions will happen ---  '
        print

    njobs = int (TOT_EVENTS) / EVENTS_PER_JOB
    print 'run : submitting ' + str (njobs) + ' jobs'

    folderName = 'MBgeneration'
    runCommand ('mkdir ' + folderName)
    jobtag = 'MB'
    for i in range (1, njobs + 1) :
        tag = jobtag + '_' + str (i)
        jobname = prepareJob (tag)
        runCommand ('bsub -J ' + tag + ' -u pippopluto -q ' + QUEUE + ' < ' + jobname, 1, TESTING == 0)
    runCommand ('mv *.job ' + folderName)


