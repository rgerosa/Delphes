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
parser.add_option("-s","--dataset"    , dest="dataset"    , type="string", help="dataset to run on")
parser.add_option("-p","--proxy"      , dest="proxy"    , type="string", help="Grid proxy for xrdcp")
parser.add_option("-w","--workdir"    , dest="workdir"    , type="string", default="mydir",   help = "Name of the directory for creating job folders")
parser.add_option("-c","--configCard" , dest="configCard" , type="string", default="Cards/CMS_Phase_I_50PileUp_Tracker2p5.tcl", help = "Delphes card to be used")
parser.add_option("-o","--outputname" , dest="outputname" , type="string", default="outtree", help = "Name of the output file. Default is: outtree")
parser.add_option("-i","--inputPUdir" , dest="inputPUdir" , type="string", help="Path where PU files are located. Example: /afs/cern.ch/user/s/spigazzi/work/public/PU_14TeV/")
parser.add_option(""  ,"--eosdir"     , dest="eosdir"     , type="string", default="",help="Name of the eos output directory for jobs")
parser.add_option("-e","--executable" , dest="executable" , type="string", default="DelphesCMSFWLite",help="Name of the executable. Default is: DelphesPythia8")
parser.add_option("-d","--executableDumper" , dest="executableDumper" , type="string", default="",help="Name of the Dumper executable. Default is: ''")

## njobs decide the number of the jobs as a function fo the number of the file in the inputdir
parser.add_option("-a","--njobmax"    , dest="njobmax"    , type="int"   , default = 0, help = "Number of jobs max")
## eventsToRun decide the number of events for each job
#parser.add_option(""  ,"--eventsPerJob", dest="eventsPerJob", default = 25000, type="int", help="Number of events in each job .. don't look at the number of files")

## other options
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

def makeFilesListDAS(dataset,workingdir):
    
    list = [];
    command = ('python scripts/das.py --limit=0 --query="file dataset=%s" > %s/list.txt' % (dataset,workingdir))
    #print command;    
    os.system(command)
    file = open('%s/list.txt'%(workingdir), 'r')     
    for line in file:
        list.append(line.replace('\n',''));
    print 'Found %d files for list.txt' %(len(list));
    return list

def writeJobs(workingdir,executable,executableDumper,configCard,outputname,eosoutdir,njobs):

    print "==> Start Job Creation ==>";
    #----------------------------------------------------------------
    # --- prepare the list of files to be analyzed
    #-----------------------------------------------------------------
    listoffiles = [];
    listoffiles = makeFilesListDAS(options.dataset,workingdir); ## make the file list for the input directory using das client script

    listofPUfiles = [];
    listofPUfiles = makeFilesList(options.inputPUdir,workingdir,"_PU",".mb");

    jobid = 0
    nentries = 0 # possible to get this? Could submit multiple jobs with ::MaxEvents / ::SkipEvents in cfg
    random.seed()

    ## loop on the list
    for ifile in listoffiles:

     ## get the number of jobs for this file       
     njobs = 1

     ## create this jobs
     print "create job for file ",jobid/njobs," total file ",len(listoffiles)," name ",ifile;

     for i in range(njobs):
      ## create job directory         
      ## avoid too many jobs
      pileupEntry = random.randint(0,len(listofPUfiles)-1)

      if jobid >= options.njobmax and options.njobmax != 0: break; 
      jobdir = '%s/JOB_%d'%(workingdir,jobid)
      os.system("mkdir -p "+jobdir)        

      ## copy the config file in the job dir
      os.system("cp %s %s/temp.tcl"%(configCard,jobdir))    
      pileupName = listofPUfiles[pileupEntry].split("/")
      configName = configCard.split("/");
      # set correct pile-up file
      os.system("cat "+jobdir+"/temp.tcl | sed -e s%MB_1.mb%"+str(pileupName[len(pileupName)-1])+"%g "+"> %s/temp2.tcl"%(jobdir))
      # remove lhe particle branch from card
      os.system("cat "+jobdir+'/temp2.tcl | sed -e s%add\ Branch\ Delphes/LHEParticles%#add\ Branch\ Delphes/LHEParticles%g > '+jobdir+"/"+configName[len(configName)-1])
      ## remove the temp card
      os.system("rm %s/temp.tcl"%jobdir);
      os.system("rm %s/temp2.tcl"%jobdir);
      
      pileUpName = listofPUfiles[pileupEntry].replace('eos/cms','').replace('\n','');
      fileName   = ifile.split("/");  
      ### prepare the job scripts
      jobscript = open('%s/subJob_%d.sh'%(jobdir,jobid),'w')
      jobscript.write('cd %s \n'%jobdir)
      jobscript.write('eval ` scramv1 runtime -sh ` \n')
      jobscript.write('cd - \n \n')
      jobscript.write('cmsStage -f %s ./ \n'%(pileUpName))
      jobscript.write('export X509_USER_PROXY=%s \n'%(options.proxy) )
      jobscript.write('xrdcp root://xrootd.unl.edu/%s ./ \n'%(ifile))
      jobscript.write('scp '+os.getcwd()+"/"+executable+" ./ \n \n"); 
      if( executableDumper != '' ) : 
          jobscript.write('scp '+os.getcwd()+"/"+executableDumper+" ./ \n"); 
      jobscript.write('if ( \n')
      jobscript.write('\t touch %s/subJob_%d.run \n'%(jobdir,jobid))
      if( executableDumper != '' ) :
          jobscript.write('\t ./%s %s %s %s \n'%(executable,jobdir+"/"+str(configName[len(configName)-1]),outputname+"_"+str(jobid)+"_Delphes.root",str(fileName[len(fileName)-1])));
          dumperSplit = executableDumper.split("/")
          jobscript.write('\t ./%s %s %s'%(dumperSplit[len(dumperSplit)-1],outputname+"_"+str(jobid)+"_Delphes.root", outputname+"_"+str(jobid)+".root")); 
      else :
          jobscript.write('\t ./%s %s %s %s'%(executable,jobdir+"/"+str(configName[len(configName)-1]),outputname+"_"+str(jobid)+".root",str(fileName[len(fileName)-1])));    
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

    return jobid ;

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
  print "==> Start the script ==>";

  # if the working director is already there remove it
  if os.path.isdir(workingdir):
   print 'Workding directory already exist. Exiting'
   sys.exit(0)
  
  # create the new working dir
  print "mkdir -p "+workingdir;
  os.system("mkdir -p "+workingdir);

  # write the output on eos
  if (options.eosdir !=''):
      eosoutdir = options.eosdir+'/'+options.workdir;
      mkdir = 'cmsMkdir ';
      command = '%s %s'%(mkdir,eosoutdir);
      print command;
      os.system(command);
  
  # create jobs
  njobs = writeJobs(workingdir,options.executable,options.executableDumper,options.configCard,options.outputname,eosoutdir,njobs);
   
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
    