#!/usr/bin/env python 
# -*- coding: utf-8 -*-

from ginfizz_preprocess import preprocess
from ginfizz_bandpass import bandpass
import ginfizz_config

#    preprocess(spm_standalone, mcr, flibasedir, atlasfile, resultdir
#    main(spm_standalone, mcr, flibasedir, atlasfile, resultdir, hpass, lpass, cthr)
#    return




if __name__ == "__main__":
    import sys 
    import os
    print sys.argv
    
    spm_standalone =  sys.argv[1] 
    mcr =  sys.argv[2] 
    flibasedir  =  sys.argv[3]
    atlasfile  =  sys.argv[4]
    roiatlasfile =  sys.argv[5]
    resultdir =  sys.argv[6]
    
    import nipype.pipeline.engine as pe          # pypeline engine
    from nipype import Node
    
    restingState = pe.Workflow(name='restingState')
    
    # working directory
    working_dir =  os.path.join(resultdir, "_report")
    restingState.base_dir = working_dir    

    # ## logging management 

    from nipype import config
    cfg = dict(logging=dict(workflow_level = 'DEBUG'),
               execution={'stop_on_first_crash': False,
                          'hash_method': 'content'})
    config.update_config(cfg)
    restingState.config['execution'] = {'stop_on_first_rerun': 'False',
                                   'hash_method': 'timestamp'}

    # create logging dir under resultdir os.path.join('dir','other-dir') os.makedirs(newpath)
    import os
    logsdir = os.path.join(resultdir, 'logs')
    try:
        if not os.path.exists(logsdir):
            os.makedirs(logsdir)
    except:
        print "exception while creating logs directory"

    from nipype import config, logging
    config.update_config({'logging': {'log_directory': logsdir,
                                      'log_to_file': True }})
    logging.update_logging(config)

    from nipype import logging
    iflogger = logging.getLogger('interface')
    message = "Start of resting state workflow"
    iflogger.info(message)            
    
    
    # 1 - preprocess images
    result_preproc = preprocess(spm_standalone, mcr, flibasedir,atlasfile, resultdir)
 
    
    # 2 - bandpass
    result_prebandpass = bandpass(resultdir)

    
    # connection
    
    restingState.connect(result_preproc, 'norm12bis.normalized_files'  ,result_prebandpass, 'inputNode.normalized_masks')
    restingState.connect(result_preproc,  'norm12.normalized_files'  ,result_prebandpass, 'inputNode.filesToMerge')
    restingState.connect(result_preproc,  'realign.realignment_parameters'  ,result_prebandpass, 'inputNode.mocoVariables')
    
    # the run
    restingState.run()
    
    restingState.write_graph(graph2use='colored', format='svg', simple_form=True)