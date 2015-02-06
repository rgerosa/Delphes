echo "Preparing to run a Delphes job on LHE files"

echo "input arguments:"
echo $@

echo $1

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


# CRAB3 job wrapper executes the following steps on the worker node:
# 1. It sets up the CMSSW environment (basically running cmsrel and cmsenv).
# 2. It tweaks the CMSSW parameter-set configuration file provided by the user 
#    (JobType.psetName parameter in the CRAB configuration file) 
#    by setting up the input files and the input lumis the job should analyze.
# 3. It executes cmsRun -j FrameworkJobReport.xml -p PSet.py, 
#    where PSet.py is the tweaked CMSSW parameter-set configuration file.
# 4. It parses the framework job report produced by cmsRun 
#    looking for OutputModule and TFileService outputs. 
#    These output files will be staged out locally on the temporary storage area 
#    of the site where the job runs, and also remotely
#    (unless the user disabled the transfers by setting General.
#    transferOutputs = False). 

# instructions from here: https://github.com/rgerosa/Delphes/blob/master/README

# fastjet intallation
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


# get the LHE file
# install fastjet
# modify the CMSSW files and re-run the sourcing
# get Delphes
# compile Delphes
# run the code

touch output.txt

INIT_FOLDER=`pwd`

# setup the proxy for using lcg-cp
# --- --- --- --- --- --- --- --- --- --- --- --- ---

#export X509_USER_PROXY="PROXY_FIXME"
#cp $X509_USER_PROXY /tmp
#voms-proxy-init --noregen 
voms-proxy-info --all

#exit 0

# do everything from the CMSSW release folder
# --- --- --- --- --- --- --- --- --- --- --- --- ---

cd $CMSSW_BASE/src

# INSTALL FASTJET
# --- --- --- --- --- --- --- --- --- --- --- --- ---

mkdir fastjet
cd fastjet
FASTJET_TGZ="fastjet-3.1.0.tar.gz" 
FASTJET_DIR=`echo $PWD/$FASTJET_TGZ | sed 's/.tar.gz//'`
cp ${INIT_FOLDER}/${FASTJET_TGZ} ./
tar fzx $FASTJET_TGZ
rm -rf $FASTJET_TGZ

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
pwd >> $INIT_FOLDER/output.txt
echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
ls >> $INIT_FOLDER/output.txt

export FASTJET_BASE=$PWD
cd $FASTJET_DIR

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
pwd >> $INIT_FOLDER/output.txt
echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
ls >> $INIT_FOLDER/output.txt

./configure --prefix=$FASTJET_BASE
make -j
make check -j
make install -j

cd ..

mv $INIT_FOLDER/the_fjcontrib.tgz ./
tar xzf the_fjcontrib.tgz
cd fjcontrib
./configure --fastjet-config=$FASTJET_BASE/bin/fastjet-config CXXFLAGS="-I$FASTJET_BASE/include -I$FASTJET_BASE/tools"
make -j
make check -j
make install -j
make fragile-shared -j
make fragile-shared-install -j
cd ..
cd ..

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
echo looking into $CMSSW_BASE/src >> $INIT_FOLDER/output.txt
ls $CMSSW_BASE/src/ >> $INIT_FOLDER/output.txt

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
cat $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml >> $INIT_FOLDER/output.txt

mv $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml ./fastjet.xml.safe
cat fastjet.xml.safe | \
 sed -e s%\<environment\ name=\"FASTJET_BASE\"\ default=\".*%\<environment\ name=\"FASTJET_BASE\"\ default=\"$FASTJET_BASE\"/\>% | \
 sed -e s%\<environment\ name=\"LIBDIR\"\ default=\".*%\<environment\ name=\"LIBDIR\"\ default=\"$FASTJET_BASE/lib\"/\>% | \
 sed -e s%\<environment\ name=\"INCLUDE\"\ default=\".*%\<environment\ name=\"INCLUDE\"\ default=\"$FASTJET_BASE/include\"/\>% > $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
cat $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml >> $INIT_FOLDER/output.txt

scram setup fastjet

# SET UP PYTHIA
# --- --- --- --- --- --- --- --- --- --- --- --- ---
echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
cat $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/pythia8.xml >> $INIT_FOLDER/output.txt

cd $CMSSW_BASE/src
cp ${INIT_FOLDER}/pythia8201.tgz ./
tar xvfz pythia8201.tgz
rm pythia8201.tgz
cd pythia8201/    
./configure --enable-shared  --with-hepmc2=/afs/cern.ch/cms/slc6_amd64_gcc472/external/hepmc/2.06.07-cms3 --with-hepmc2-lib=/afs/cern.ch/cms/slc6_amd64_gcc472/external/hepmc/2.06.07-cms3/lib/ --with-hepmc2-include=/afs/cern.ch/cms/slc6_amd64_gcc472/external/hepmc/2.06.07-cms3/include --with-lhapdf5=$LHAPATH --with-lhapdf5-lib=$LHAPATH/../../../lib --with-lhapdf5-include=$LHAPATH/../../../include --with-lhapdf5-bin=$LHAPATH/../../../bin
make
export PYTHIA_BASE=$PWD 
cd ..
mv  $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/pythia8.xml ./pythia8.xml.safe
cat  ./pythia8.xml.safe | \
    sed -e s%\<environment\ name=\"PYTHIA8_BASE\"\ default=\".*%\<environment\ name=\"PYTHIA8_BASE\"\ default=\"$PYTHIA_BASE\"/\>% | \
    sed -e "s%<lib name=\"pythia8tohepmc\"/>%<lib name=\"pythia8lhapdf5\"/>%g" | \
    sed -e "s%$PYTHIA8_BASE/xmldoc%$PYTHIA8_BASE/share/Pythia8/xmldoc/%g" >   $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/pythia8.xml

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
cat $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/pythia8.xml # >> $INIT_FOLDER/output.txt

scram setup pythia8
eval `scramv1 runtime -sh`
export SRT_PYTHIA8DATA_SCRAMRTDEL=$PYTHIA8DATA


echo "--------+--------+--------+--------+--------+--------"
echo DONE pythia installation
echo "--------+--------+--------+--------+--------+--------"


# INSTALL DELPHES
# --- --- --- --- --- --- --- --- --- --- --- --- ---

cd $CMSSW_BASE/src
mv $INIT_FOLDER/Delphes.tgz ./
tar xzf Delphes.tgz
cd Delphes

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt 
pwd >> $INIT_FOLDER/output.txt

#./configure
#mv Makefile Makefile.safe
#cat Makefile.safe | \
# sed -e s%FASTJET_BASE.*%FASTJET_BASE\ =$FASTJET_BASE% | \
# sed -e s%FASTJET_LIB.*%FASTJET_LIB\ =fastjetplugins\ fastjettools\ siscone\ siscone_spherical\ fastjet% | \
# sed -e s%FASTJET_LIBDIR.*%FASTJET_LIBDIR\ =$FASTJET_BASE/lib% | \
# sed -e s%FASTJET_INCLUDE.*%FASTJET_INCLUDE\ =$FASTJET_BASE/include% > Makefile

echo "--------+--------+--------+--------+--------+--------"
echo "--------+--------+--------+--------+--------+--------"
printenv
echo "--------+--------+--------+--------+--------+--------"
echo "--------+--------+--------+--------+--------+--------"
#echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
#cat Makefile >> $INIT_FOLDER/output.txt

#echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
# make -j >>& $INIT_FOLDER/output.txt
#make -j 

cd ..

# END OF INSTALLATION, still in $CMSSW_BASE/src
# --- --- --- --- --- --- --- --- --- --- --- --- ---

#cp output.txt delphesTree.root

#mv $INIT_FOLDER/CMS_Phase_I_140PileUp_Tracker2p5.tcl ./
#mv $INIT_FOLDER/phamom.dat.gz ./
#gunzip phamom.dat.gz
#cp Delphes/DelphesPythia8 ./
#./DelphesPythia8 CARD LHEFILE OUTPUTFILE MJJCUT FILTRAFULLYHAD STARTINGEVENT NUMEVENTS

#
# copy input file
#

INFILE=$(tail -n+$1 $INIT_FOLDER/inputlist.txt | head -n+1 )
echo $INFILE
lcg-cp  -b -D srmv2 srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=${INFILE} ./infile.lhe

#
# copy pileup file
#

PUFILE=$(tail -n+$1 $INIT_FOLDER/pulist.txt | head -n+1 )
echo $PUFILE
lcg-cp  -b -D srmv2 srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=${PUFILE} ./MB_1.mb


$CMSSW_BASE/src/Delphes/DelphesPythia8 $CMSSW_BASE/src/Delphes/Cards/TP_CARDS/CMS_Phase_I_140PileUp_Tracker4p0.tcl infile.lhe delphesTree.root 0 0 0 

cp delphesTree.root $INIT_FOLDER

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
pwd >> $INIT_FOLDER/output.txt
echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
ls >> $INIT_FOLDER/output.txt

cd $INIT_FOLDER

echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
pwd >> $INIT_FOLDER/output.txt
echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt
ls >> $INIT_FOLDER/output.txt
echo "--------+--------+--------+--------+--------+--------" >> $INIT_FOLDER/output.txt


echo 'DONE' >> $INIT_FOLDER/output.txt
