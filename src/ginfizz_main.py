#!/usr/bin/env python 
# -*- coding: utf-8 -*-

from ginfizz_preprocess import preprocess
from ginfizz_bandpass import bandpass

#    preprocess(spm_standalone, mcr, flibasedir, atlasfile, resultdir
#    main(spm_standalone, mcr, flibasedir, atlasfile, resultdir, hpass, lpass, cthr)
#    return




if __name__ == "__main__":
    import sys 
    import os
    print sys.argv
    
    resultdir =  sys.argv[5]
    
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
    result_preproc = preprocess(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
 
    
    # 2 - bandpass
    hpass =  0.01
    lpass = 0.1 
    acqNb = 240
    result_prebandpass = bandpass(resultdir, hpass, lpass, acqNb)

    
    # connection
    
    restingState.connect(result_preproc, 'norm12bis.normalized_files'  ,result_prebandpass, 'inputNode.normalized_masks')
    restingState.connect(result_preproc,  'norm12.normalized_files'  ,result_prebandpass, 'inputNode.filesToMerge')
    restingState.connect(result_preproc,  'realign.realignment_parameters'  ,result_prebandpass, 'inputNode.mocoVariables')
    
    # the run
    restingState.run()
    
    restingState.write_graph(graph2use='colored', format='svg', simple_form=True)