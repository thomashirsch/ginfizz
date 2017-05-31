#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ginfizz_config


# # Start of the new subworkflow - PREBANDPASS


# -------------------------------------------
#
#     Start of bandpass pipeline

def bandpass(resultdir):
    '''
    Description:
        the  functional bandpass pipeline
            - plots mocoVariables
            - computes regressors
            - FSL merges normalized_files
            - AFNI 3D bandpass merged_file
    
    Inputs:
        - resultdir: the results directory
        - hight pass
        - low pass
        - functionnal images acquisition number
        
    Outputs:
        - functionnal bandpassed signal time courses

    '''    

    import nipype.interfaces.io as nio           # Data i/o
    import nipype.interfaces.spm as spm          # spm
    import nipype.interfaces.matlab as mlab      # how to run matlab
    import nipype.interfaces.utility as util     # utility
    import nipype.pipeline.engine as pe          # pypeline engine
    
    from nipype.interfaces.utility import Function, IdentityInterface
    
    from nipype.interfaces.fsl.maths import MathsCommand
    from nipype.interfaces.fsl.utils import PlotMotionParams   # to plot moco variables
    
    import ginfizz_config
    
    print ginfizz_config.hpass, ginfizz_config.lpass
    
    # creation of a subworflow to calculate the bandpass parameters
    prebandpass = pe.Workflow(name='prebandpass')
    
    # 0 - we fix the working dir to get the commands later
    import os
    prebandpass.base_dir = os.path.join(resultdir, "_report")    
    
    # ## logging management 

    from nipype import config
    cfg = dict(logging=dict(workflow_level = 'DEBUG'),
               execution={'stop_on_first_crash': False,
                          'hash_method': 'content'})
    config.update_config(cfg)
    prebandpass.config['execution'] = {'stop_on_first_rerun': 'True',
                                        'hash_method': 'timestamp'}

    # create logging dir under resultdir os.path.join('dir','other-dir') os.makedirs(newpath)
    import os
    logsdir = os.path.join(resultdir, 'logs')
    
    from nipype import config, logging
    config.update_config({'logging': {'log_directory': logsdir,
                                      'log_to_file': True }})
    logging.update_logging(config)

    from nipype import logging
    iflogger = logging.getLogger('interface')
    message = "Start of bandpass workflow"
    iflogger.info(message)   


    # todo remplace this node by an identity node that get input from preprocess pipeline / node 7 eg. fsl merge
    # input node get the good files first the tissues normalised files
    # the files to merge are in /scratch/user/hirsch/datadir/data_results/functionnal
    #sourcedir = '/scratch/user/hirsch/datadir4/data_results_py/structural/normalized_files'
    # second the merged functionnal file
    #sourcemergeddir = '//scratch/user/hirsch/datadir4/data_results_py/functionnal/normalized_files'
    #funcdir = '/scratch/user/hirsch/datadir4/data_results_py/functionnal/realignment_parameters'

    from nipype import Node

    # intput node
    field_list=['normalized_masks', 
                'filesToMerge',
                'mocoVariables']


    inputNode = Node(IdentityInterface(fields=field_list), name="inputNode")


    # to erase this just for memory while transforming code
    #templates = dict(wmMask=sourcedir+ "/" + "wc2*.nii",
                     #csfMask=sourcedir+ "/" + "wc3*.nii",
                     #filesToMerge=sourcemergeddir + "/" + "wcra*.nii",
                     #mocoVariables=funcdir+ "/" + "rp*.txt")
    
    #filesource = Node(SelectFiles(templates), "filesource")
    #filesource.inputs.subject_id = "subj1"
    #filesource.outputs.get()



    # ## 1 - compute moco file to feed with ortho.txt file the bandpass node of preprocess workflow 
    # derivatives
    #dx(t)= x(t) - x(t-1)
    #x(t0) = x(t1)
       


    # ## Node 1 - compute moco
 

    
    def computeMoco(mocoFile):
        import pandas as pd
        import numpy as np
        import os
        import ginfizz_config
        from ginfizz_tools  import plusDerivative, plusSquare
    
        # read the moco file to put it in a panda dataframe 
        mocodf = pd.read_csv(mocoFile, header=None, sep='  ',engine='python') 
        
        # first, we derivate the 6 colunms of dataframe of moco file, and append the 6 new colums to df
        dfderivate = plusDerivative(mocodf)
        
        # then, we compute the square of the 6 first colums, and append them to df. it makes 18 colums that are going to 
        # participate in the ortho file to make 18 regressors bandpassed
        dfsquare = plusSquare(dfderivate)
        
        #     out_file = os.path.join(os.getcwd(), "merged_array.txt")
        # np.savetxt(out_file, cvectors, fmt="%.10f")

        h = os.path.join(os.getcwd(), "ortho.txt")
        
        #g = dfsquare.to_csv('ortho.txt', sep=' ',index=False, header=False)
        np.savetxt(h, dfsquare, fmt=str("%.10f"))

        return h
    

    computeMoco = Node(Function(input_names=['mocoFile'],
                                    output_names=['out_file'],
                                    function=computeMoco),
                                    name='computeMoco')
    
    
    
    prebandpass.connect(inputNode, "mocoVariables", computeMoco, "mocoFile")


    # ##  2 - Moco plots
    
    
    # plot moco variables rotations
    MotionCorrectionPlot1 = Node(PlotMotionParams(), name="MotionCorrectionPlot1")
    MotionCorrectionPlot1.inputs.ignore_exception = False     
    MotionCorrectionPlot1.inputs.in_source = 'spm'     
    MotionCorrectionPlot1.inputs.output_type = 'NIFTI_GZ'     
    MotionCorrectionPlot1.inputs.plot_size = (500, 1000)     
    MotionCorrectionPlot1.inputs.plot_type = 'rotations'     
    MotionCorrectionPlot1.inputs.terminal_output = 'stream'     
    prebandpass.connect(inputNode, "mocoVariables", MotionCorrectionPlot1, "in_file")

    
    # plot moco variables translations
    MotionCorrectionPlot2 = Node(PlotMotionParams(), name="MotionCorrectionPlot2")
    MotionCorrectionPlot2.inputs.ignore_exception = False     
    MotionCorrectionPlot2.inputs.in_source = 'spm'     
    MotionCorrectionPlot2.inputs.output_type = 'NIFTI_GZ'     
    MotionCorrectionPlot2.inputs.plot_size = (500, 1000)     
    MotionCorrectionPlot2.inputs.plot_type = 'translations'     
    MotionCorrectionPlot2.inputs.terminal_output = 'stream'     
    prebandpass.connect(inputNode, "mocoVariables", MotionCorrectionPlot2, "in_file")


    # plot moco variables displacement
    MotionCorrectionPlot3 = Node(PlotMotionParams(), name="MotionCorrectionPlot3")
    MotionCorrectionPlot3.inputs.ignore_exception = False     
    MotionCorrectionPlot3.inputs.in_source = 'spm'     
    MotionCorrectionPlot3.inputs.output_type = 'NIFTI_GZ'     
    MotionCorrectionPlot3.inputs.plot_size = (500, 1000)     
    MotionCorrectionPlot3.inputs.plot_type = 'displacement'     
    MotionCorrectionPlot3.inputs.terminal_output = 'stream'     
    prebandpass.connect(inputNode, "mocoVariables", MotionCorrectionPlot3, "in_file")


    # ## 3 - get wm and csf mask (normalised), erode them and calculate signal mean on both masks
    # (input from segment + normalyse wmMask = '/homes_unix/hirsch/_pypipe/datadir/data_results/structural/norm_files/wc2t0009_t1_s03.nii')
    # remark: erosion is done once; in connectomics, it is done 3 times
    
    def regexfilter(files_list,patern):
        import re
        
        for f in files_list:
            if re.search(patern, str(f)):
                    res = f         
        return res            
    
    # calculate eroded binary mask for wm
    from nipype.interfaces.fsl.maths import MathsCommand

    erosion = pe.Node(interface=MathsCommand(), name='erosion')
        
    erosion.inputs.args = '-thr 0 -uthr 111 -bin -ero  '     
    erosion.inputs.ignore_exception = False     
    erosion.inputs.output_type = 'NIFTI'     
    erosion.inputs.terminal_output = 'stream'     
    prebandpass.connect(inputNode,('normalized_masks',regexfilter,r'wc2.*nii'), erosion, "in_file")
    

    # calculate eroded binary mask for csf
    
    def regexfilter(files_list,patern):
        import re
        
        for f in files_list:
            if re.search(patern, str(f)):
                    res = f         
        return res             
    
    erosioncsf = pe.Node(interface=MathsCommand(), name='erosioncsf')
           
        
    erosioncsf.inputs.args = '-thr 0 -uthr 111 -bin -ero '     
    erosioncsf.inputs.ignore_exception = False     
    erosioncsf.inputs.output_type = 'NIFTI'     
    erosioncsf.inputs.terminal_output = 'stream'     
    prebandpass.connect(inputNode,('normalized_masks', regexfilter,r'wc3.*nii') , erosioncsf, "in_file")
    
    # lets merge the functionnal normalysed files with fsl merge
    
    from nipype.interfaces.fsl import Merge
    fsl_merge = pe.Node(interface=Merge(), name="fsl_merge")     
     
    fsl_merge.inputs.dimension = 't'
    fsl_merge.inputs.output_type = 'NIFTI_GZ'
    # todo recuperer le TR deja exploité du xml
    fsl_merge.inputs.tr = ginfizz_config.TR
    #fsl_merge.inputs.force_even = False     
    fsl_merge.inputs.ignore_exception = False     
    #fsl_merge.inputs.sort = False
    
    prebandpass.connect(inputNode,  'filesToMerge', fsl_merge, 'in_files')

    # lets calculate the mean signal on these eroded masks first on wm
    from nipype.interfaces.fsl.utils import ImageMeants
    
    wmMeants = Node(ImageMeants(), name="wmMeants")     
    wmMeants.inputs.ignore_exception = False     
    wmMeants.inputs.order = 1     
    wmMeants.inputs.output_type = 'NIFTI_GZ'     
    wmMeants.inputs.terminal_output = 'stream'     
    prebandpass.connect(fsl_merge, "merged_file" , wmMeants, "in_file")   
    prebandpass.connect(erosion, "out_file", wmMeants, "mask")


    # lets calculate the mean signal on these eroded masks econd on csf
    
    csfMeants = Node(ImageMeants(), name="csfMeants")     
    csfMeants.inputs.ignore_exception = False     
    csfMeants.inputs.order = 1     
    csfMeants.inputs.output_type = 'NIFTI_GZ'     
    wmMeants.inputs.terminal_output = 'stream'     
    prebandpass.connect(fsl_merge, "merged_file" , csfMeants, "in_file")   
    prebandpass.connect(erosioncsf, "out_file", csfMeants, "mask")

    from computations import  mergeTables
    from flow import createList3
    
    # lets make a list from the 3 file moco+derivative+square, wmmeants, csfmeants
    # Node: func2.createListFrom3
    createListFrom3 = Node(Function(input_names=['item1', 'item2', 'item3'],                                     
                                    output_names=['out_list'],                                     
                                    function=createList3),                            
                                    name="createListFrom3")     
    createListFrom3.inputs.ignore_exception = False     
    prebandpass.connect(computeMoco, "out_file", createListFrom3, "item1")     
    prebandpass.connect(wmMeants, "out_file", createListFrom3, "item2")     
    prebandpass.connect(csfMeants, "out_file", createListFrom3, "item3")         

    # Node: func2.buildNuisanceTable     
    buildNuisanceTable = Node(Function(input_names= ['in_files'],                                        
                                       output_names=['out_file'],                                        
                                       function=mergeTables),                               
                                       name="buildNuisanceTable")     
    buildNuisanceTable.inputs.ignore_exception = False     
    prebandpass.connect(createListFrom3, "out_list", buildNuisanceTable, "in_files")

  
    
    # lets calculate the band pass with 3 ortho files moco+derivative+square, wm_meants and csf_meants
    from nipype.interfaces import afni

    afniBandpass = pe.Node(interface=afni.Bandpass(), name='afniBandpass')
      
    afniBandpass.inputs.automask = False     
    afniBandpass.inputs.environ = {}     
    afniBandpass.inputs.highpass = ginfizz_config.hpass     
    afniBandpass.inputs.ignore_exception = False     
    afniBandpass.inputs.lowpass = ginfizz_config.lpass     
    afniBandpass.inputs.outputtype = 'NIFTI_GZ'     
    afniBandpass.inputs.terminal_output = 'stream'   
    
    prebandpass.connect(fsl_merge, "merged_file", afniBandpass, "in_file")
    prebandpass.connect(buildNuisanceTable, "out_file", afniBandpass, "orthogonalize_file")        

    # lets make an output node with the 3 files to be pass to afni bandpass: 
    # file 1 with 1 colums moco, derivative moco, moco and derivative square
    # file 2 with wm mean signal on wm
    # file 3 with lcf mean signal on lcf

    field_list=['rotations_plot', 
                'translations_plot',
                'displacement_plot',
                'moco_gradient_square',
                'wm_normalized_eroded_mask',
                'wm_meants',
                'csf_normalized_eroded_mask',
                'csf_meants',
                "merged_file",
                'bandpassedFile']
                
            
    outputNode = Node(IdentityInterface(fields=field_list), name="outputNode")
    
    # for plot files
    prebandpass.connect(MotionCorrectionPlot1,  'out_file', outputNode, 'rotations_plot')
    prebandpass.connect(MotionCorrectionPlot2,  'out_file', outputNode, 'translations_plot')
    prebandpass.connect(MotionCorrectionPlot3,  'out_file', outputNode, 'displacement_plot')
    
    # for moco , gradient and square
    prebandpass.connect(computeMoco,  'out_file', outputNode, 'moco_gradient_square')

    # for segmented normalised eroded wm and lcf mask
    prebandpass.connect(erosion,  'out_file', outputNode, 'wm_normalized_eroded_mask')
    prebandpass.connect(erosioncsf,  'out_file', outputNode, 'csf_normalized_eroded_mask')
    
    # for wm and csf mean signal to text files in functionnal repository
    prebandpass.connect(wmMeants,  'out_file', outputNode, 'wm_meants')
    prebandpass.connect(csfMeants,  'out_file', outputNode, 'csf_meants')
    
    prebandpass.connect(fsl_merge, "merged_file", outputNode, "merged_file")
    
    prebandpass.connect(afniBandpass,  'out_file', outputNode, 'bandpassedFile')
    
    
    # data sink
    datasink = pe.Node(nio.DataSink(), name='datasink')
    
    
    datasink.inputs.base_directory = resultdir
    
    # for motion correction  plot files
    prebandpass.connect(MotionCorrectionPlot1,  'out_file', datasink, 'functionnal.motion_plot_files.rotation')
    prebandpass.connect(MotionCorrectionPlot2,  'out_file', datasink, 'functionnal.motion_plot_files.translation')
    prebandpass.connect(MotionCorrectionPlot3,  'out_file', datasink, 'functionnal.motion_plot_files.displacement')

    # for segmented normalised eroded wm and lcf mask
    prebandpass.connect(erosion,  'out_file', datasink, 'functionnal.bandpass_wm_mask')
    prebandpass.connect(erosioncsf,  'out_file', datasink, 'functionnal.bandpass_csf_mask')
    
    # for wm and lcf mean signal to text files in functionnal repository
    prebandpass.connect(wmMeants,  'out_file', datasink, 'functionnal.bandpass_wm_meants')
    prebandpass.connect(csfMeants,  'out_file', datasink, 'functionnal.bandpass_csf_meants')
    
    # fsl merged file 
    prebandpass.connect(fsl_merge, "merged_file", datasink, "functionnal.merged_file")
    
    # orthogonalize_file file 
    prebandpass.connect(buildNuisanceTable, "out_file", datasink, "functionnal.orthogonalize_file")    
    
    # bandpassed merged normalised  file sinks
    prebandpass.connect(afniBandpass,  'out_file', datasink, 'functionnal.bandpassedFile')

    # the run
    #prebandpass.run()
    
    prebandpass.write_graph(graph2use='colored', format='svg', simple_form=True)
    
    return prebandpass

