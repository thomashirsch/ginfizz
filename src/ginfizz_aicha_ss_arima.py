#!/usr/bin/env python 
# coding: utf-8

# cette version est faite pour traiter le signal bp sans arima

# nouvelles specs on part d'un nouvel atlas de 10 régions et on voit quels sont les voxels qui sont en correlations avec ces régions  <br> 
# on va d'abord determiner combien il y de regions, puis on va les calculer à partir de l'atlas utilisateur
def identifyRegionAtlases(atlas_file):
        """compute Regions of Interest number, their integer values from atlas file given by the user,
        then compute the atlas image for each of n values, then
        output is the array of regions atlases, each one representing a different ROI """

        import numpy as np     
        import nibabel as nib     
        import os 

        #------------------------------
        # fcts used by current function
        #------------------------------

        #    
        def selectRegion(atlas_file, n):
                ''' lets select region i input is atlasFile, output is maskFile with region n  '''

                atlas_img=nib.load(atlas_file)         
                atlas_array=np.asarray(atlas_img.dataobj).copy() # Avoid caching the proxy image
                region_array = atlas_array
                # binary mask the resulting image
                region_array[atlas_array <> n] = 0
                region_array[atlas_array == n] = 1
                region_image = nib.Nifti1Image(region_array, atlas_img.affine, atlas_img.header)
                result = os.getcwd() + '/' + 'region_image' + str(n) + '.nii'
                nib.save(region_image, result)
                return result   

        def identifyRegionNb(atlas_file):
                """compute Regions of Interest number, their integer values from atlas file given by the user, 
                output is the array of int number, each one representing a different ROI """

                atlas_img=nib.load(atlas_file)         
                atlas_array=np.asarray(atlas_img.dataobj).copy() # Avoid caching the proxy image
                atlas_array = atlas_array[atlas_array>0]
                # np.unique(a)
                regions = np.unique(atlas_array)
                print regions
                regionsNb = len(regions) 
                return regions


        #------------------------------
        # end fcts used by current function
        #------------------------------    

        roi_img_list =[]
        roi_nbs = identifyRegionNb(atlas_file)
        for i in roi_nbs:

                roi_img_list.append(selectRegion(atlas_file, i))

        return roi_img_list

# --------------------------------------------------------------------

atlas_file = '/scratch/user/hirsch/datadir/data_set/t0009/repos01/Atlases/atlas_2reg.nii'

regions_images_list = identifyRegionAtlases(atlas_file)

# --------------------------------------------------------------------

import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.matlab as mlab      # how to run matlab
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine

from nipype.interfaces.utility import Function, IdentityInterface

from nipype.interfaces.fsl.maths import MathsCommand
from nipype.interfaces.fsl.utils import PlotMotionParams   # to plot moco variabl
from nipype import SelectFiles, Node, MapNode



# creation of a subworflow to process ROI Correlations
connectivity = pe.Workflow(name='connectivity')

# comprendre bp dir 
arimadir = '/scratch/user/hirsch/datadir4/data_results_py/functionnal/bandpassedFile'

resultdir = '/scratch/user/hirsch/datadir4/data_results_py/'


from nipype import SelectFiles, Node
templates = dict(arimaFile=arimadir+ "/" + "*.nii.gz",
                 normalized_c1_file=resultdir+ "/" + "structural/normalized_files/wc1*.nii",
                 normalized_c2_file=resultdir+ "/" + "structural/normalized_files/wc2*.nii")

filesource = Node(SelectFiles(templates), "filesource")
filesource.inputs.subject_id = "subj1"
filesource.outputs.get()

# lets compute the brain mask c1 + c2 threshold at 0.2 to compute later the target time courses

from nipype.interfaces.fsl import MultiImageMaths

addFiles = pe.Node(interface=MultiImageMaths(), name='addFiles')

addFiles.inputs.op_string = "-add %s"   
addFiles.inputs.output_datatype = 'short'
addFiles.inputs.ignore_exception = False     
addFiles.inputs.output_type = 'NIFTI'     
addFiles.inputs.terminal_output = 'stream'     

connectivity.connect(filesource, 'normalized_c1_file', addFiles, "in_file")
connectivity.connect(filesource, 'normalized_c2_file' , addFiles, "operand_files")

from nipype.interfaces.fsl import Threshold

thrFile = pe.Node(interface=Threshold(), name='thrFile')

thrFile.inputs.thresh = 0.2   
thrFile.inputs.ignore_exception = False     
thrFile.inputs.output_type = 'NIFTI'     
thrFile.inputs.terminal_output = 'stream'     

connectivity.connect(addFiles,"out_file" , thrFile, "in_file")






# here comes the map node, we have to iterate the treatment for each region identified by previous step
# b = pe.MapNode(interface=B(), name="b", iterfield=['in_file']) 
# http://nipype.readthedocs.io/en/latest/users/mapnode_and_iterables.html

# lets calculate the mean residu signal eg. time courses in every region 
# lets use iterables -> startnode.iterables = ('subject_id', subjects)

from nipype.interfaces.fsl.utils import ImageMeants

regMeants = Node(ImageMeants(), name="regMeants")  

regMeants.iterables = ('mask', regions_images_list)

regMeants.inputs.ignore_exception = False     
regMeants.inputs.order = 1     
regMeants.inputs.output_type = 'NIFTI_GZ'     
regMeants.inputs.terminal_output = 'stream'     
connectivity.connect(filesource, "arimaFile" , regMeants, "in_file")   


# correlation computations
def computeCorrelations(residus, residusRegMean, brainMask):
    
        '''Function that takes 3 parameters
           residus: the residu of arima time courses
           residuRegMean: the residus mean on determined region i
           brainMask: the mask of gm and wm thresholded at 0,2
           Computes pearson correlations between the seed eg. residuRegMean
           and all the voxels of the brainMask
           Returns an Nifti images containing the correlations'''
        
        import os
        import numpy as np
        import matplotlib.pyplot as plt
        import nibabel as nib
        from scipy.stats.stats import pearsonr
        
        # first we get the seed mean signal
        seed_ts_array = np.loadtxt(residusRegMean)
        
        # from an other hand we get the residus 4D matrix
        fmri_data=nib.load(residus) 
        fmri_array=np.asarray(fmri_data.dataobj)
        
        # we get the coordinnates of voxels in all gm and wm normalized todo 
        reg_data=nib.load(brainMask) 
        regarray=np.asarray(reg_data.dataobj)
        # transpose(nonzero(a))
        reg_coords = np.transpose(np.nonzero(regarray))
        volume_shape = reg_coords.shape
        print volume_shape
        coords = list(np.ndindex(volume_shape))
        print len(coords)
        
        # the we iterate the correlation calculation on all voxels of the brain mask
    
        # the correlation matrix is initialized with all values to 0 
        corr_matrix = np.full(reg_data.shape, 0, dtype=float)
    
        for i in range(reg_coords.shape[0]):
            target_array = fmri_array[reg_coords[i, 0], reg_coords[i, 1],reg_coords[i,2], :]
            #print target_array
            non_zero_nb = np.count_nonzero(target_array)
            
            # if target time courses are all null, we do not compute correlation
            if non_zero_nb:
                try:
                    p = pearsonr(seed_ts_array,target_array) 
                    if p[0] > 0.5:
                        print p[0]
                    corr_matrix[reg_coords[i, 0], reg_coords[i, 1],reg_coords[i,2]] = p[0] 
                except:
                    print "exception"
        
        # save matrix in a file
        # create the resulting image
        corr_image = nib.Nifti1Image(corr_matrix,affine=reg_data.affine, header=reg_data.header)
        # save the correlation array
        out_file = os.getcwd() + '/' + 'corr_regv5.nii'
        nib.save(corr_image, out_file)
        
        return out_file
        
# node to compute the correlation matrix as each region i of user atlas mean signal serves as a seed 
# and brain mask signal residuserves as the targets
# def computeCorrelations(residus, residusRegMean, brainMask):
correlationsComputeNode = Node(Function(input_names=['residus', 'residusRegMean', 'brainMask'],
                                output_names=['out_file'],
                                function=computeCorrelations),
                                name='correlationsComputeNode')

connectivity.connect(filesource, "arimaFile", correlationsComputeNode, "residus")
connectivity.connect(regMeants, "out_file", correlationsComputeNode, "residusRegMean")
connectivity.connect(thrFile, "out_file", correlationsComputeNode, "brainMask")

# data sink
datasink = pe.Node(nio.DataSink(), name='datasink')
datasink.inputs.base_directory = '/scratch/user/hirsch/datadir4/data_results_py'

# data sink brain mask to compute target time courses
connectivity.connect(thrFile, "out_file", datasink, 'structural.normalized_c1c2_file')

# data sink of mean signal on ROIs that will be used as seed signal
connectivity.connect(regMeants,  'out_file', datasink, 'functionnal.regMeants')

# data sink of correlations matrix for ROI i
connectivity.connect(correlationsComputeNode,  'out_file', datasink, 'functionnal.correlationMatrix')

connectivity.run()


# # fin du pipe

