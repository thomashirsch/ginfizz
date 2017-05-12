#!/usr/bin/env python 
# coding: utf-8

# de V2 à V3 travail sur les logs , les exceptions et produire une image .nii.gz au lieu d'un fichier .npy

# # Modèle ARIMA en R sur tous les voxels de la matière grise 
# maintenant sur tous les voxels

#----------------------------------------------------------------------

def arima(sourcedir, resultdir):
    """compute ARIMA model and residu on all brain time courses with previous bandpass filtering"""
    
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
    import os


    # ## Workflow start
    
    
    # creation of a subworflow to calculate the arima residus
    arimaResidus = pe.Workflow(name='arimaResidus')
    
    # 0 - we fix the working dir to get the commands later
    
    arimaResidus.base_dir = os.path.join(resultdir, "_report")    


    # ## logging management 
    
    from nipype import config
    cfg = dict(logging=dict(workflow_level = 'DEBUG'),
               execution={'stop_on_first_crash': False,
                           'hash_method': 'content'})
    config.update_config(cfg)
    arimaResidus.config['execution'] = {'stop_on_first_rerun': 'True',
                                       'hash_method': 'timestamp'}

    from nipype import config, logging
    config.update_config({'logging': {'log_directory': resultdir,
                          'log_to_file': True }})
    logging.update_logging(config)
    
    from nipype import logging
    iflogger = logging.getLogger('interface')
    message = "Start of arima "
    iflogger.info(message)


    from nipype import SelectFiles, Node
    templates = dict(gmMask=sourcedir+ "/structural/normalized_files/" + "wc1*.nii",
                     bpfile=sourcedir+ "/functionnal/bandpassedFile/" + "wcra*_merged_bp.nii.gz")
    
    filesource = Node(SelectFiles(templates), "filesource")
    filesource.inputs.subject_id = "subj1"
    filesource.outputs.get()


    def computeArimaResidu(signalFile):
        
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
         
        bp =nib.load(signalFile) 
        bparray = np.asarray(bp.dataobj).copy()
        result = np.zeros(bparray.shape,  dtype=np.float)
        
        nb_errors_arima=0
        
        print bparray.shape[0]
        print bparray.shape[1]
        print bparray.shape[2]
        # lets iterate on all the points of the bp image in 3D
        for i in range(bparray.shape[0]):
            
            for j in range(bparray.shape[1]):
                
                for k in range(bparray.shape[2]):
                    
                    # time serie associated to point i
                    bp_ts = bparray[i, j, k, : ]
                    #print bp_ts
        
                    # convert numpy array to panda
                    l = len(bp_ts)
                    index = ['Row'+str(m) for m in range(1, l+1)]
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
                        result[i, j, k, : ] = residu
                    except:
                        nb_errors_arima += 1
                                                      
        
        iflogger.info( "End iteration on " + str(i) + " ,"  + str(j) + " ,"  + str(k) + " ," )
        iflogger.info( "nb errors arima : " + str(nb_errors_arima))                     
        # create the resulting image
        residu_image = nib.Nifti1Image(result, bp.affine, bp.header)
        # save the residu array
        out_file = os.getcwd() + '/' + 'arima_residu.nii.gz'
        nib.save(residu_image, out_file)
        
        return out_file
    
 

    computeArimaResiduNode = Node(Function(input_names=['signalFile'],
                                    output_names=['out_file'],
                                    function=computeArimaResidu),
                                    name='computeArimaResiduNode')
    
    computeArimaResiduNode.inputs.ignore_exception = True
    
    arimaResidus.connect(filesource,  'bpfile', computeArimaResiduNode, "signalFile")

    # output node 
    outputNode = Node(IdentityInterface(fields=[ 'arimaResidus']), name="outputNode")
    
    arimaResidus.connect(computeArimaResiduNode,  'out_file', outputNode, 'functionnal.arima.arimaResidus')
    
    
    ## data sink
    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = resultdir
    
    ## for dumping result files
    
    arimaResidus.connect(computeArimaResiduNode,  'out_file', datasink, 'functionnal.arima.arimaResidus')
    
    # The RUN
    arimaResidus.run()
    
    return


if __name__ == "__main__":
    import sys 
    print sys.argv
    arima(sys.argv[1], sys.argv[2])

