#!/usr/bin/python

import sys
import os
import commands
from commands import getstatusoutput
import argparse
import datetime
import math
import random
import string

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def replaceParameterInFile (inputFile, outputFile, substitute): 
    f = open (inputFile)
    s = f.read ()
    f.close ()
    for k,v in substitute.items () :
        s = s.replace (k, v)
    f = open (outputFile, 'w')
    f.write (s)
    f.close ()


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def runCommand (command, printIt = 0, doIt = 1) :
    if printIt : print ('> ' + command)
    if doIt : 
        commandOutput = commands.getstatusoutput (command)
        if printIt : print commandOutput[1]
        return commandOutput[0]
    else :    print ('    jobs not submitted')
    return 1


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'

if __name__ == '__main__':

    debug = True

    parser = argparse.ArgumentParser (description = 'run phantom productions on lxplus')
    parser.add_argument('-d', '--dir'       , default= '', help = 'eos folder with input files'     )
# - delphes zip
# - delphes card
# - LHE file    
    args = parser.parse_args ()

    # get the proxy filename
    retCode = getstatusoutput('voms-proxy-info --path')
    proxypath=retCode[1]
    retcode = getstatusoutput('grep config.General.workArea '+args.dir+"/crabConfig.py")
    crabarea = string.strip(retcode[1].split("=")[1]," '")
    #print crabarea
    retcode = getstatusoutput('grep config.General.requestName '+args.dir+"/crabConfig.py")
    requestname = string.strip(retcode[1].split("=")[1]," '")
    #print requestname

    retcode = getstatusoutput('grep NJOBS '+args.dir+"/crabConfig.py")
    njobs = int(string.strip(retcode[1].split("=")[1].split("#")[0]," '"))
    #print njobs

    retcode = getstatusoutput('grep config.Data.outLFN '+args.dir+"/crabConfig.py")  
    storagebase = string.strip(retcode[1].split("=")[1]," '")
    #print storagebase

    retcode = getstatusoutput('grep config.Data.primaryDataset '+args.dir+"/crabConfig.py")  
    pd = string.strip(retcode[1].split("=")[1]," '")
    #print storagebase
    
    #generate detailed crab report
    runCommand("crab status --long --proxy="+proxypath+" --d "+crabarea+"/crab_"+requestname+" > crabresult.txt")
    
    retcode = getstatusoutput('grep "Task name:" crabresult.txt')
    submission = retcode[1].split()[2].split("_")[0] + "_" +  retcode[1].split()[2].split("_")[1]
    #print submission
    
    #readin extended carb log
    crabfile = open("crabresult.txt")
    foundlist = False
    nproc = 0
    running=[]
    toresubmit=[]
    todelete=[]
    finished=[]
    srmbase="srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=/pnfs/desy.de/cms/tier2"
    for line in crabfile:
        if (not foundlist):
            if len(line.split()) >= 2 and  line.split()[0] == "Job":
                foundlist = True
            #print "dismiss: " + line
            continue        
        else:
            if len(line.split()) < 3:
                break
            #print  len(line.split()) 
            print "processing: " + string.strip(line,"\n")
            nproc = nproc+1
            secs = line.split()
            jobnr=secs[0]
            status = secs[1]
            if status == "running":
                print " still running => nothing to do"
                running.append(jobnr)
            if status == "failed":
                #print "XXX",
                retcode = getstatusoutput('lcg-ls -D srmv2 -b '+srmbase+storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+jobnr+".root")
                #print retcode
                print " failed => add to resubmission list"
                toresubmit.append(jobnr)
                if not(retcode[1].find("No such file or directory") > 0):
                    todelete.append(jobnr)
                    print " Job Nr. "+jobnr+" failed, but output file: "+srmbase+storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+jobnr+".root exists"  
            if status =="finished":
                retcode = getstatusoutput('grep "\-#######  read events:" '+args.dir+'/pass/cmsRun-stdout-'+jobnr+'.log')
                if(retcode[1].find("-#######  read events:") >= 0):
                    print " logfile Ok => checking presence of tree"
                    retcode = getstatusoutput('lcg-ls -D srmv2 -b '+srmbase+storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+jobnr+".root")
                    if (retcode[1].find("No such file or directory") > 0):
                        print " file "+srmbase+storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+jobnr+".root missing"
                        toresubmit.append(jobnr)
                    else:
                        finished.append(jobnr)
                else:
                    print " logfile NOT OK  => mark for resuibmission"
                    toresubmit.append(jobnr)
                    retcode = getstatusoutput('lcg-ls -D srmv2 -b '+srmbase+storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+jobnr+".root")

                    if not (retcode[1].find("No such file or directory") > 0):
                        todelete.append(jobnr)
                        print " Job Nr. "+jobnr+" failed, but output file: "+srmbase+storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+jobnr+".root exists"  

    print "+++++++++++++++++++++++++++++"
    print "processsed "+ str(nproc)+ " jobs"
    if (nproc !=njobs):
        print "WARNING: number of processed jobs != number of requested jobs"
    if (len(toresubmit)+len(finished)+len(running))!=njobs:
        print "WARNING no action for some jobs????"
    print "DONE: "+str(len(finished))+" jobs: ",
    print finished
    print "RESUBMIT",
    for num in toresubmit:
        sys.stdout.write(num+",")
    
    print
    
    delfile = open("delete.txt","w")
    for job in todelete:
        delfile.write(srmbase+storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+job+".root\n")
    delfile.close()
    print "corrupt root files written to delete.txt"

    delfile = open("finished.txt","w")
    for job in finished:
        delfile.write(storagebase+"/"+pd+"/crab_"+requestname+"/"+submission+"/0000/delphesTree_"+job+".root\n")
    delfile.close()
    print "finished root files written to finished.txt"


