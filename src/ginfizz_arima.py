#!/usr/bin/env python 
# coding: utf-8

# de V2 à V3 travail sur les logs , les exceptions et produire une image .nii.gz au lieu d'un fichier .npy

# # Modèle ARIMA en R sur tous les voxels de la matière grise 

import nibabel as nib
import pandas as pd
import numpy as np

import nibabel as nib

import matplotlib.pylab as plt

import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.matlab as mlab      # how to run matlab
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine

from nipype.interfaces.utility import Function, IdentityInterface


# ## Workflow start


# creation of a subworflow to calculate the arima residus
arimaResidus = pe.Workflow(name='arimaResidus')


# ## logging management 


from nipype import config
cfg = dict(logging=dict(workflow_level = 'DEBUG'),
           execution={'stop_on_first_crash': False,
                       'hash_method': 'content'})
config.update_config(cfg)
arimaResidus.config['execution'] = {'stop_on_first_rerun': 'True',
                                   'hash_method': 'timestamp'}
import os
from nipype import config, logging
config.update_config({'logging': {'log_directory': os.getcwd(),
                      'log_to_file': True }})
logging.update_logging(config)

from nipype import logging
iflogger = logging.getLogger('interface')
message = "Start of arima "
iflogger.info(message)

# todo
sourcedir = '/scratch/user/hirsch/datadir4/data_results_py'


from nipype import SelectFiles, Node
templates = dict(gmMask=sourcedir+ "/structural/normalized_files/" + "wc1*.nii",
                 bpfile=sourcedir+ "/functionnal/bandpassedFile/" + "wcra*_merged_bp.nii.gz")

filesource = Node(SelectFiles(templates), "filesource")
filesource.inputs.subject_id = "subj1"
filesource.outputs.get()



def computeGm(gmMask):
    
    import nibabel as nib
    import numpy as np
    import os
    
    from nipype import logging
    
    # on regarde le grey matter
    i1=nib.load(gmMask)         
    i1array=np.asarray(i1.dataobj).copy() 
    i1array[(i1array )< 0.2] = 0
    # binary mask the resulting image
    i1array[i1array > 0] = 1
    gm_coord = np.transpose(np.nonzero(i1array))
    #print "nb of gm voxels" + str(len(gm_coord))
    iflogger = logging.getLogger('interface')
    iflogger.info("nb of gm voxels" + str(len(gm_coord)))
    
    out_file = os.getcwd() + '/' + 'gm_coord_file.npy'
    np.save(out_file, gm_coord)
    return out_file
    

# identify all the gm voxels
computeGmVoxels = Node(Function(input_names=['gmMask'],
                                output_names=['out_file'],
                                function=computeGm),
                                name='computeGmVoxels')

arimaResidus.connect(filesource, "gmMask", computeGmVoxels, "gmMask")


def computeArimaResidu(gmPointsFile, signalFile):
    
    import nibabel as nib
    import numpy as np
    import pandas as pd
    import os
    
    from nipype import logging
    iflogger = logging.getLogger('interface')
    
    import readline
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    ro.r('library(stats)')
   
    gmPoints = np.load(gmPointsFile)
    print gmPoints.size
    print gmPoints[0]
    
    bp =nib.load(signalFile) 
    bparray = np.asarray(bp.dataobj).copy()
    
    result = np.zeros(bparray.shape,  dtype=np.float)
    
    nb_errors_arima=0
    
    for i in range(len(gmPoints)):
        
        #print gmPoints[i][0], gmPoints[i][1],gmPoints[i][2]
    
        # time serie associated to point i
        bp_ts = bparray[gmPoints[i][0], gmPoints[i][1],gmPoints[i][2], : ]
        #print bp_ts
    
        # covert numpy array to panda
        l = len(bp_ts)
        index = ['Row'+str(j) for j in range(1, l+1)]
        df = pd.DataFrame(bp_ts, index=index)
        #print df
        rdf = pandas2ri.py2ri(df)
        ro.globalenv['r_timeseries'] = rdf
        
        try:
            # model the time serie with ARIMA
            ro.r('fit <- arima(r_timeseries, order = c(1,1,1))')
            # get the residu
            residu = ro.r('fit$residuals')
            # result update
            result[gmPoints[i][0], gmPoints[i][1],gmPoints[i][2], : ] = residu
        except:
            nb_errors_arima += 1
            iflogger.info("i, and exception arima / nb errors arima :" + str(i) + " ,  " + str(nb_errors_arima))
            #print "exception arima"
    
    iflogger.info( "End iteration on i" )
    iflogger.info( "nb errors arima : " + str(nb_errors_arima))
      
    # create the resulting image
    residu_image = nib.Nifti1Image(result, bp.affine, bp.header)
    # save the residu array
    out_file = os.getcwd() + '/' + 'arima_residu.nii.gz'
    nib.save(residu_image, out_file)
    
    return out_file
    
 
# identify all the gm voxels
computeArimaResiduNode = Node(Function(input_names=['gmPointsFile', 'signalFile'],
                                output_names=['out_file'],
                                function=computeArimaResidu),
                                name='computeArimaResiduNode')

computeArimaResiduNode.inputs.ignore_exception = True

arimaResidus.connect(computeGmVoxels,  'out_file', computeArimaResiduNode, "gmPointsFile")
arimaResidus.connect(filesource,  'bpfile', computeArimaResiduNode, "signalFile")

# output node 
outputNode = Node(IdentityInterface(fields=['gmVoxels', 'arimaResidus']), name="outputNode")
arimaResidus.connect(computeGmVoxels,  'out_file', outputNode, 'functionnal.arima.gmVoxels')
arimaResidus.connect(computeArimaResiduNode,  'out_file', outputNode, 'functionnal.arima.arimaResidus')


## data sink
datasink = pe.Node(nio.DataSink(), name='datasink')
datasink.inputs.base_directory = '/scratch/user/hirsch/datadir4/data_results_py'

## for dumping result files
arimaResidus.connect(computeGmVoxels,  'out_file', datasink, 'functionnal.arima.gmVoxels')
#arimaResidus.connect(computeArimaResiduNode,  'out_file', datasink, 'functionnal.arima.arimaResidus')

# The RUN
arimaResidus.run()

