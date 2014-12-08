import os
import sys
import numpy
import glob
import random
import ROOT
import commands
from ROOT import TFile, TTree

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inputdir"   , dest="inputdir"   , type="string", help="Path where ntuples or LHE are located. Example: /store/user/govoni/LHE/phantom/TP14TeV/")
parser.add_option("-k","--fileKey"    , dest="fileKey"    , type="string", help="LHE keys to find them in the path")
parser.add_option("-w","--workdir"    , dest="workdir"    , type="string", default="mydir",   help = "Name of the directory for creating job folders")
parser.add_option("-o","--outputname" , dest="outputname" , type="string", default="outtree", help = "Name of the output file. Default is: outtree")
parser.add_option(""  ,"--eosdir"     , dest="eosdir"     , type="string", default="",help="Name of the eos output directory for jobs")
parser.add_option("-e","--executable" , dest="executable" , type="string", default="DelphesPythia8",help="Name of the executable. Default is: DelphesPythia8")

## njobs decide the number of the jobs as a function fo the number of the file in the inputdir
parser.add_option("-a","--njobs"    , dest="njobs"    , type="int"   , default = 0, help = "Number of jobs")

parser.add_option(""  ,"--checkJobs"  , dest="checkJobs"  , action="store_true", default=False,help="Checks job status")
parser.add_option(""  ,"--resubmit"   , dest="resubmit"   , action="store_true", default=False,help="Resubmit job ")
parser.add_option(""  ,"--submit"     , dest="submit"     , action="store_true", default=False,help="submit jobs after create them ")
parser.add_option("-j","--jobmin"     , dest="jobmin"     , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option("-x","--jobmax"     , dest="jobmax"     , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option("-q","--queue"      , dest="queue"      , type="string", default="1nh",help="Name of the queue on lxbatch")


(options,args)=parser.parse_args()

eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'

def makeFilesList(inputdir,workingdir,additionalString = "",key=".root"):
    
    list = [];
    command = ('%s find -f %s | grep %s > %s/list%s.txt' % (eos,inputdir,key,workingdir,additionalString));
    print command;    
    os.system(command);
    file = open('%s/list%s.txt'%(workingdir,additionalString), 'r');     
    for line in file:
       list.append(line.replace('/eos/cms/','root://eoscms.cern.ch//').replace('\n',''));
    print 'Found %d files for list%s.txt' %(len(list),additionalString);    
    return list

def writeJobs(workingdir,executable,inputdir,outputname,eosoutdir,njobs):

    print "==> Star Job Creation ==>";
    #----------------------------------------------------------------
    # --- prepare the list of files to be analyzed --> read from eos
    #-----------------------------------------------------------------
    listoffiles = [];
    listoffiles = makeFilesList(inputdir,workingdir,"",options.fileKey); ## make the file list for the input directory on eos
    
    #---------------------------------------------
    # --- now split the jobs
    #---------------------------------------------
    for job in range(njobs):

     jobdir = '%s/JOB_%d'%(workingdir,job)
     os.system("mkdir -p "+jobdir)

     #--- prepare the list of files for each job
     f = open('%s/input_%d.txt'%(jobdir,job), 'w')
     sublist = [file for i,file in enumerate(listoffiles) if (i%njobs==job)]
     for fname in sublist:
        f.write('%s \n'%fname)


     ### prepare the job scripts                                                                                                                                                          
     jobscript = open('%s/subJob_%d.sh'%(jobdir,job),'w')
     jobscript.write('cd %s \n'%jobdir)
     jobscript.write('eval ` scramv1 runtime -sh ` \n')
     jobscript.write('cd - \n')
     jobscript.write('scp '+os.getcwd()+"/"+executable+" ./ \n");     
     jobscript.write('if ( \n')
     jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,job))      
     jobscript.write('\t ./%s %s/input_%d.txt %s_%d.root'%(executable,jobdir,job,outputname,job))
     jobscript.write(') then \n')
     if (eosoutdir == ''):
       jobscript.write('\t cp ./%s_%d.root %s \n'%(outputname,job,jobdir))
     else:
       jobscript.write('\t cmsStage -f ./%s_%d.root %s/ \n'%(outputname,job,eosoutdir))
     jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,job))
     jobscript.write('else \n')
     jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,job))
     jobscript.write('fi \n')
     os.system('chmod a+x %s/subJob_%d.sh'%(jobdir,job))

    print "==> End of Job creation ==>";
    return njobs ;

def submitJobs(workingdir, njobs, queue):
    for job in range(njobs):
        print 'job %d' %job
        jobdir  = '%s/JOB_%d'%(workingdir,job)
        jobname = '%s/subJob_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/subJob_%d.log -e %s/subJob_%d.err %s'%(queue,jobdir,job,jobdir,job,jobname)
        os.system('bsub -q %s -o %s/subJob_%d.log -e %s/subJob_%d.err %s'%(queue,jobdir,job,jobdir,job,jobname))
        

def checkJobs(workingdir):
    ## find total number of jobs 
    jobs = glob.glob( '%s/JOB_*/subJob*.sh'% (workingdir) )
    print 'Total number of jobs: %d' %len(jobs)
        
    ## list of job that are done    
    listdone = [];
    for j in range(len(jobs)):
       if os.path.isfile('%s/JOB_%d/subJob_%d.done' % (workingdir,j,j)):
          listdone.append(j);

    print 'Total number of DONE jobs: %s ' % len(listdone)
        
    ## eliminate run jobs for jobs that are done    
    for j in listdone:
        f = '%s/JOB_%d/subJob_%d.run'%(workingdir,j,j)
        if (os.path.isfile(f)):
            os.system('rm %s'%f)
            
    ## print running jobs    
    listrun = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/subJob_%d.run' % (workingdir,j,j))]
    print 'Total number of RUNNING jobs: %d ' %len(listrun)
         
    ## print list failed
    listfailed = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/subJob_%d.fail' % (workingdir,j,j))]
    print 'Failed jobs: %s ' % listfailed
    print '   %s' %listfailed
                        
    for j in listfailed:
     print 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/subJob_%d.sh -o -e %s/JOB_%d/subJob_%d.err '%(queue,workingdir,j,j,workingdir,j,j,workingdir,j,j)
     os.system('bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/subJob_%d.sh -o -e %s/JOB_%d/subJob_%d.err '%(queue,workingdir,j,j,workingdir,j,j,workingdir,j,j));                    
    
#-----------------------------------------
#--- MAIN CODE
#-----------------------------------------


path       = os.getcwd()
workingdir = path+'/'+options.workdir
eosoutdir  = ''
njobs      = 0

if not options.checkJobs and not options.resubmit: ## if you want to create jobs

  print "==> Star the script ==>";

  # if the working director is already there remove it
  if os.path.isdir(workingdir):
   print "rm -r "+workingdir;
   os.system("rm -r "+workingdir);

  
  # create the new working dir
  print "mkdir -p "+workingdir;
  os.system("mkdir -p "+workingdir);

  # wrote the output on eos
  if (options.eosdir !=''):
      eosoutdir = options.eosdir+'/'+options.workdir;
      mkdir = 'cmsMkdir ';
      command = '%s %s'%(mkdir,eosoutdir);
      print command;
      os.system(command);
  
  # create jobs
  njobs = writeJobs(workingdir,options.executable,options.inputdir,options.outputname,eosoutdir,options.njobs);
   
  # -- submit jobs
  print "==> Start job submission ==>";
  if not options.resubmit:
   if options.submit: 
    submitJobs(workingdir, njobs, options.queue);
   # resubmit option
   else:
    print "submit options is not enables ";
  elif options.resubmit and options.jobmin >-1 and options.jobmax >-1 and not options.checkJobs:
    for ijob in range(options.jobmax-options.jobmin):
     print 'Resubmitting job %d ' %(ijob+options.jobmin)
     resubcmd = 'bsub -q %s -o %s/JOB_%d/subJob_%d.log %s/JOB_%d/subJob_%d.sh'%(options.queue,workingdir,options.jobmin+ijob,options.jobmin+ijob,workingdir,options.jobmin+ijob,options.jobmin+ijob)
    #print resubcmd
    os.system(resubcmd)

elif options.checkJobs:
    checkJobs(workingdir);
    
