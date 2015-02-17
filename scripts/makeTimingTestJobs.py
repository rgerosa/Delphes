#!/usr/bin/python

import commands
import fileinput
import argparse
import sys
import os

TESTING = 0

rootfolder = os.getcwd ()

def runCommand (command, printIt = 0, doIt = 1) :
    if printIt : print ('> ' + command)
    if doIt : 
        commandOutput = commands.getstatusoutput (command)
        if printIt : print commandOutput[1]
        return commandOutput[0]
    else :    print ('    jobs not submitted')
    return 1
    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def replaceAll (file, searchExp, replaceExp) :
    for line in fileinput.input (file, inplace = 1) :
        if searchExp in line:
            line = line.replace (searchExp, replaceExp)
        sys.stdout.write (line)


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def prepareJob (folder, delphesCommand) :
    filename = folder + '/run.job'
    f = open (filename, 'w')

    f.write ('cd ' + folder + '    \n')
    f.write ('scram setup fastjet  \n')
    f.write ('scram setup pythia8  \n')
    f.write ('eval `scram run -sh` \n')

    f.write (delphesCommand + '\n')

    f.close ()
    return filename


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == "__main__":

    parser = argparse.ArgumentParser (description = 'run delphes on AFS for timing studies')
    parser.add_argument ('-p', '--punum'      , default = '1',     help = 'number of PU events [1]')
    parser.add_argument ('-z', '--puzpos'     , default = '0',     help = 'PU z position in mm [0]')
    parser.add_argument ('-s', '--puzsmear'   , default = '0',     help = 'PU z smearing in mm [0]')
    parser.add_argument ('-t', '--putsmear'   , default = '0',     help = 'PU t smearing in ns [0]')
    parser.add_argument ('-c', '--delphescard', default = 'NONE',  help = 'delphes card [NONE]')
    parser.add_argument ('-i', '--lhefile'    , default = 'NONE',  help = 'input LHE file [NONE]')
    parser.add_argument ('-T', '--totEvents'  , default = '10000', help = 'total number of events to be generated [10000]')
    parser.add_argument ('-q', '--lsfQueue'   , default = '2nw',   help = 'LSF queue [2nw]')

    args = parser.parse_args ()
    
    folder = 'job_' + args.punum + '_' + args.puzpos + '_' + args.puzsmear + '_' + args.putsmear
    folder = folder.replace ('.', 'p')
    folder = rootfolder + '/' + folder
    
    # if the working directory is already there remove it
    if os.path.isdir (folder) :
        runCommand ("rm -r " + folder, 1)
    runCommand ("mkdir -p " + folder, 1)
    
    card = folder + '/card.tcl'
    runCommand ('cp ' + args.delphescard + ' ' + card)

    replaceAll (card, 'TEMP_PUMEAN',    args.punum)
    replaceAll (card, 'TEMP_PUZPOS',    args.puzpos)
    replaceAll (card, 'TEMP_PUZSPREAD', args.puzsmear)
    replaceAll (card, 'TEMP_PUTSPREAD', args.putsmear)
    replaceAll (card, 'TEMP_OUTFOLDER', folder)

    outfile = folder + '/output.root' ;
    delphesCommand = rootfolder  + '/DelphesPythia8 ' + card \
        + ' ' + args.lhefile + ' ' + outfile + ' 0 0 1 ' + args.totEvents

    jobname = prepareJob (folder, delphesCommand)
    runCommand ('bsub -u pippopluto -q ' + args.lsfQueue + ' < ' + jobname, 1, TESTING == 0)


    