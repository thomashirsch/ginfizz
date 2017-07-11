#!/bin/bash

# env var for dev machine
# MCRROOT=/homes_unix/hirsch/historique_fli_iam/essai_spm_stand_alone/mcr2016/v91
# SPMSAROOT=/homes_unix/hirsch/historique_fli_iam/essai_spm_stand_alone/spm12/run_spm12.sh
# CODEROOT=/homes_unix/hirsch/ginfizz/src
# DATAROOT=/scratch/user/hirsch/datadirp

# env var for docker machines, for VIP
MCRROOT=/opt/mcr/v91/
SPMSAROOT=/opt/spm12/run_spm12.sh
CODEROOT=/rstp_code
DATAROOT=/gasw-execution-dir/rstp_data



function die {
    local D=`date`
    echo "[ $D ] ERROR: $*"
    exit 1
}

function info {
    local D=`date`
    echo "[ $D ] INFO: $*"
}

if [ $# != 2 ]
then
    die "usage: $0 <input_tgz>  <output_dir>"
fi

info "Start running ginfizz wrapper"

INPUTFILE=$1
OUTPUTDIR=$2
info "parameters are ${INPUTFILE} and ${OUTPUTDIR} "
# it is lacking the path of the tar command

info "PATH is ${PATH}"
info "INPUTFILE is ${INPUTFILE}"

#DATAROOT= ${OUTPUTDIR}

cd ${DATAROOT}

pwd

# untar of the inputfile 
TOP=`tar --exclude '*/*' --overwrite -tzf ${INPUTFILE}` || die "Cannot get top level directory from ${INPUTFILE}!"
info "TOP is ${TOP}"
tar zxf ${INPUTFILE} || die "Cannot untargz ${INPUTFILE}!"

info "TOP is ${TOP}"
export TOP

# List top directory in subject directory. There must be only 1. 
NSECONDTOP=`ls ${DATAROOT}/${TOP} | wc -l| awk '{print $1}'` || die "Cannot count directories in ${DATAROOT}/${TOP}!"
if [ ${NSECONDTOP} -ne 1 ] 
then
    die "Found 0 or more than 1 directory in ${PWD}!"
fi

SECONDTOP=`ls ${DATAROOT}/${TOP} || die "Cannot find directory in ${DATAROOT}/${TOP}"`
info "SECONDTOP is ${SECONDTOP}"

TROISTOP=`ls ${DATAROOT}/${TOP}/${SECONDTOP} || die "Cannot find directory in ${DATAROOT}/${TOP}/${SECONDTOP}"`
info "TROISTOP is ${TROISTOP}"


# Export FLI base directory expected by the matlab script. 
export FLIBASEDIR=${DATAROOT}/${TOP}/${SECONDTOP}/${TROISTOP}
info "FLIBASEDIR is ${FLIBASEDIR}"

# Export ATLAS base directory expected by the matlab script. 
export ATLASBASEDIR=${DATAROOT}/${TOP}/${SECONDTOP}/${TROISTOP}/Atlases
info "ATLASBASEDIR is ${ATLASBASEDIR}"

# modification for roi Atlas
# Export ROI ATLAS base directory expected by the matlab script. 
export ROIATLASBASEDIR=${DATAROOT}/${TOP}/${SECONDTOP}/${TROISTOP}/ROIAtlases
info "ROIATLASBASEDIR is ${ROIATLASBASEDIR}"


# find the XML file
# Search for the XML file of the subject
NXMLFILES=`ls ${FLIBASEDIR}/*.xml | wc -l | awk '{print $1}'` || die "Cannot count xml files in ${FLIBASEDIR}!"
if [ ${NXMLFILES} -ne 1 ] 
then
    die "Found 0 or more than 1 xml file in ${FLIBASEDIR}!"
fi
XMLFILE=`ls ${FLIBASEDIR}/*.xml` || die "Cannot find xml file in ${FLIBASEDIR}!"
info "XMLFILE is ${XMLFILE}"

# find the Atlas file
# Search for the .nii file of the subject
ATLASFILE=`ls ${ATLASBASEDIR}/*.nii | wc -l | awk '{print $1}'` || die "Cannot count nii files in ${ATLASBASEDIR}!"
if [ ${ATLASFILE} -ne 1 ] 
then
    die "Found 0 or more than 1 nii file in ${ATLASBASEDIR}!"
fi
ATLASFILE=`ls ${ATLASBASEDIR}/*.nii` || die "Cannot find nii file in ${ATLASBASEDIR}!"
info "ATLASFILE is ${ATLASFILE}"

# modification for roi Atlas
# find the Atlas file
# Search for the .nii file of the subject
ROIATLASFILE=`ls ${ROIATLASBASEDIR}/*.nii | wc -l | awk '{print $1}'` || die "Cannot count nii files in ${ROIATLASBASEDIR}!"
if [ ${ROIATLASFILE} -ne 1 ] 
then
    die "Found 0 or more than 1 nii file in ${ROIATLASBASEDIR}!"
fi
ROIATLASFILE=`ls ${ROIATLASBASEDIR}/*.nii` || die "Cannot find nii file in ${ROIATLASBASEDIR}!"
info "ROIATLASFILE is ${ROIATLASFILE}"





# preparation of the results
cd ${DATAROOT};
mkdir results;
RESULTSDIR=${DATAROOT}/results;
info "RESULTSDIR is ${RESULTSDIR}"
export RESULTSDIR;

    #spm_standalone =  sys.argv[1] 
    #mcr =  sys.argv[2] 
    #flibasedir  =  sys.argv[3]
    #atlasfile  =  sys.argv[4]
    #roiatlasfile =  sys.argv[5]
    #resultdir =  sys.argv[6]

# 1 - make the python preprocess run
(cd ${CODEROOT};
pwd;
exec   python ./ginfizz_main.py  ${SPMSAROOT}   ${MCRROOT}  ${FLIBASEDIR}   ${ATLASFILE}   ${ROIATLASFILE} ${RESULTSDIR} ;
info "1 - python ginfizz_main sent") &&



# 4 - create a tarball from the results then give it to VIP by the outputdir argument
( info "RESULTSDIR is ${RESULTSDIR}";
  cd  ${RESULTSDIR}; 
  pwd;
  ls ;
  
  info  "we are sending the tar fct for the following repository"
  
  tar czf ${OUTPUTDIR}/data_results.tar.gz   ${RESULTSDIR}/functionnal ${RESULTSDIR}/structural ${RESULTSDIR}/logs  ${DATAROOT}/*.log;
  info "4 eval has been sent we get the results in a tarball we give the tarball to VIP";)

info "End running ginfizz wrapper"