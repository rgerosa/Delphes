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
parser.add_option("-k","--lheKey"     , dest="lheKey"     , type="string", help="LHE keys to find them in the path")
parser.add_option("-w","--workdir"    , dest="workdir"    , type="string", default="mydir",   help = "Name of the directory for creating job folders")
parser.add_option("-c","--configCard" , dest="configCard" , type="string", default="Cards/CMS_Phase_I_50PileUp_Tracker2p5.tcl", help = "Delphes card to be used")
parser.add_option("-o","--outputname" , dest="outputname" , type="string", default="outtree", help = "Name of the output file. Default is: outtree")
parser.add_option("-p","--inputPUdir" , dest="inputPUdir" , type="string", help="Path where PU files are located. Example: /afs/cern.ch/user/s/spigazzi/work/public/PU_14TeV/")
parser.add_option(""  ,"--eosdir"     , dest="eosdir"     , type="string", default="",help="Name of the eos output directory for jobs")
parser.add_option("-e","--executable" , dest="executable" , type="string", default="DelphesPythia8",help="Name of the executable. Default is: DelphesPythia8")
parser.add_option("-d","--executableDumper" , dest="executableDumper" , type="string", default="",help="Name of the Dumper executable. Default is: ''")

## njobs decide the number of the jobs as a function fo the number of the file in the inputdir
parser.add_option("-a","--njobmax"    , dest="njobmax"    , type="int"   , default = 0, help = "Number of jobs max")
## eventsToRun decide the number of events for each job
parser.add_option(""  ,"--eventsPerJob", dest="eventsPerJob", default = 25000, type="int", help="Number of events in each job .. don't look at the number of files")

## other options
parser.add_option("-m" ,"--mjjcut" , dest="mjjcut", default=0, type="int",help="cut value in GeV for the mjj at LHE")
parser.add_option("-f" ,"--filter" , dest="filter", default=1, type="int",help="filter all hadronic events at LHE level")

parser.add_option(""  ,"--checkJobs"  , dest="checkJobs"  , action="store_true", default=False,help="Checks job status")
parser.add_option(""  ,"--resubmit"   , dest="resubmit"   , action="store_true", default=False,help="Resubmit job ")
parser.add_option(""  ,"--submit"     , dest="submit"     , action="store_true", default=False,help="submit jobs after create them ")
parser.add_option("-j","--jobmin"     , dest="jobmin"     , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option("-x","--jobmax"     , dest="jobmax"     , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option("-q","--queue"      , dest="queue"      , type="string", default="1nh",help="Name of the queue on lxbatch")


(options,args)=parser.parse_args()

eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'

def makeFilesList(inputdir,workingdir,additionalString = "",key=".dat"):
    
    list = [];
    command = ('%s find -f %s | grep %s > %s/list%s.txt' % (eos,inputdir,key,workingdir,additionalString));
    print command;    
    os.system(command);
    file = open('%s/list%s.txt'%(workingdir,additionalString), 'r');     
    for line in file:
       list.append(line.replace('/eos/cms/','eos/cms/').replace('\n',''));
    print 'Found %d files for list%s.txt' %(len(list),additionalString);    
    return list

def writeJobs(workingdir,executable,executableDumper,configCard,inputdir,inputPUdir,outputname,eosoutdir,njobs):

    print "==> Star Job Creation ==>";
    #----------------------------------------------------------------
    # --- prepare the list of files to be analyzed --> read from eos
    #-----------------------------------------------------------------
    listoffiles = [];
    listoffiles = makeFilesList(inputdir,workingdir,"",options.lheKey); ## make the file list for the input directory on eos
    
    listofPUfiles = [];    
    listofPUfiles = makeFilesList(inputPUdir,workingdir,"_PU",".mb");
    
    #---------------------------------------------
    # --- now split the jobs
    #---------------------------------------------

    jobid = 0 ;
    nentries = 0;
    ## loop on the list
    for ifile in listoffiles:
     ## discover how many LHE events are there     
     res = commands.getstatusoutput("(./countEvents "+ifile+")");

     if res[0] == 0 : 
        print "not able to count lhe events --> continue "; 
        continue ;
     nentries = int(res[1]);

     ## get the number of jobs for this file       
     if nentries/options.eventsPerJob - int(nentries/options.eventsPerJob) > 0.5 :  
         njobs    = int(nentries/options.eventsPerJob)+1;  
     else:
         njobs    = int(nentries/options.eventsPerJob);  

     residualEvents = nentries-options.eventsPerJob*njobs ;
     if residualEvents > 0 : njobs = njobs+1;
     
     ## random number for pileup file      
     random.seed();
     pileupEntry = random.randint(0,len(listofPUfiles)-1);

     ## create this jobs
     print "create job for file ",jobid/njobs," total file ",len(listoffiles)," name ",ifile;

     for i in range(njobs):
      ## create job directory         
      ## avoid to many jobs
      if jobid >= options.njobmax and options.njobmax != 0: break; 
      jobdir = '%s/JOB_%d'%(workingdir,jobid)
      os.system("mkdir -p "+jobdir)        

      firstEvent = 0;
      if i != 0: 
          firstEvent = options.eventsPerJob*i;

      ## copy the config file in the job dir and sobstitute the pileup file
      os.system("cp %s %s/temp.tcl"%(configCard,jobdir));                       
      pileupName = listofPUfiles[pileupEntry].split("/"); 
      configName = configCard.split("/"); 
      os.system("cat "+jobdir+"/temp.tcl | sed -e s%/afs/cern.ch/user/s/spigazzi/work/public/PU_14TeV/MinBias_14TeV_100k_TunePP15.pileup%"+str(pileupName[len(pileupName)-1])+"%g > "+jobdir+"/"+configName[len(configName)-1]);     
      ## remove the temp card
      os.system("rm %s/temp.tcl"%jobdir);

      pileUpName = listofPUfiles[pileupEntry].replace('eos/cms','').replace('\n','');
      fileSplit  = ifile.replace('eos/cms','').replace('\n','');
      fileName   = fileSplit.split("/");  
      ### prepare the job scripts
      jobscript = open('%s/subJob_%d.sh'%(jobdir,jobid),'w')
      jobscript.write('cd %s \n'%jobdir)
      jobscript.write('eval ` scramv1 runtime -sh ` \n')
      jobscript.write('cd - \n')
      jobscript.write('cmsStage -f %s ./ \n'%(pileUpName)) 
      jobscript.write('cmsStage -f %s ./ \n'%(fileSplit))
      jobscript.write('scp '+os.getcwd()+"/"+executable+" ./ \n"); 
      if( executableDumper != '' ) : 
          jobscript.write('scp '+os.getcwd()+"/"+executableDumper+" ./ \n"); 
      jobscript.write('if ( \n')
      jobscript.write('\t touch %s/subJob_%d.run \n'%(jobdir,jobid))
      if( executableDumper != '' ) :
          jobscript.write('\t ./%s %s %s %s %d %d %d %d \n'%(executable,jobdir+"/"+str(configName[len(configName)-1]),str(fileName[len(fileName)-1]),outputname+"_"+str(jobid)+"_Delphes.root",options.mjjcut,options.filter,firstEvent,options.eventsPerJob));
          dumperSplit = executableDumper.split("/")
          jobscript.write('\t ./%s %s %s'%(dumperSplit[len(dumperSplit)-1],outputname+"_"+str(jobid)+"_Delphes.root", outputname+"_"+str(jobid)+".root")); 
      else :
          jobscript.write('\t ./%s %s %s %s %d %d %d %d'%(executable,jobdir+"/"+str(configName[len(configName)-1]),str(fileName[len(fileName)-1]),outputname+"_"+str(jobid)+".root",options.mjjcut,options.filter,firstEvent,options.eventsPerJob));    
      jobscript.write(') then \n')
      if (eosoutdir == ''):
            jobscript.write('\t cp ./%s_%s.root %s \n'%(outputname,jobid,jobdir))
      else:
            jobscript.write('\t cmsStage -f ./%s_%d.root %s/ \n'%(outputname,jobid,eosoutdir))
      jobscript.write('\t touch %s/subJob_%d.done \n'%(jobdir,jobid))
      jobscript.write('else \n')
      jobscript.write('\t touch %s/subJob_%d.fail \n'%(jobdir,jobid))
      jobscript.write('fi \n')
      os.system('chmod a+x %s/subJob_%d.sh'%(jobdir,jobid))

      jobid += 1; 

     if jobid >= options.njobmax and options.njobmax != 0: break; 

    njobs = jobid;
    
    print "==> End of Job creation ==>";

    return njobs ;

def submitJobs(workingdir, njobs, queue):
    for job in range(njobs):
        print 'job %d' %job
        jobdir  = '%s/JOB_%d'%(workingdir,job)
        jobname = '%s/subJob_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/subJob_%d.log -e %s/subJob_%d.err %s'%(queue,jobdir,job,jobdir,job,jobname)
        os.system('bsub -q %s -o %s/subJob_%d.log -e %s/subJob_%d.err %s'%(queue,jobdir,job,jobdir,job,jobname))
        


def checkJobs(wdir):

    ## find total number of jobs 
    jobs = glob.glob( '%s/JOB_*/subJob*.sh'% (wdir) )
    print 'Total number of jobs: %d' %len(jobs)
        
    ## list of job that are done    
    listdone = [];
    for j in range(len(jobs)):
       if os.path.isfile('%s/JOB_%d/subJob_%d.done' % (wdir,j,j)):
          listdone.append(j);

    print 'Total number of DONE jobs: %s ' % len(listdone)
        
    ## eliminate run jobs for jobs that are done
    
    for j in listdone:
        f = '%s/JOB_%d/subJob_%d.run'%(wdir,j,j)
        if (os.path.isfile(f)):
            os.system('rm %s'%f)
            
    ## print running jobs    
    listrun = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/subJob_%d.run' % (wdir,j,j))]
    print 'Total number of RUNNING jobs: %d ' %len(listrun)
         
    ## print list failed
    listfailed = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/subJob_%d.fail' % (wdir,j,j))]
    print 'Failed jobs: %s ' % listfailed
    print '   %s' %listfailed
                        
    for j in listfailed:
     print 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/subJob_%d.sh -o -e %s/JOB_%d/subJob_%d.err '%(queue,wdir,j,j,wdir,j,j,wdir,j,j)
     os.system('bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/subJob_%d.sh -o -e %s/JOB_%d/subJob_%d.err '%(queue,wdir,j,j,wdir,j,j,wdir,j,j));                    
    
#-----------------------------------------
#--- MAIN CODE
#-----------------------------------------


path       = os.getcwd()
configCard = path+'/'+options.configCard
workingdir = path+'/'+options.workdir
eosoutdir = ''
njobs     = 0

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
  njobs = writeJobs(workingdir,options.executable,options.executableDumper,options.configCard,options.inputdir,options.inputPUdir,options.outputname,eosoutdir,njobs);
   
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
    
