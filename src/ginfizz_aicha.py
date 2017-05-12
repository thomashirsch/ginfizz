#!/usr/bin/env python 
# coding: utf-8

# nouvelles specs on part d'un nouvel atlas de 10 régions et on voit quels sont les voxels qui sont en correlations avec ces régions  <br> 
# on va d'abord determiner combien il y de regions, puis on va les calculer à partir de l'atlas utilisateur



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


arimadir = '/scratch/user/hirsch/datadir4/data_results_py/functionnal/arima/arimaResidus'
atlasdir = '/scratch/user/hirsch/datadir/data_set/t0009/repos01/Atlases'


from nipype import SelectFiles, Node
templates = dict(arimaFile=arimadir+ "/" + "*.nii.gz",
                 atlasFile=atlasdir + "/" + "atlas*.nii")

filesource = Node(SelectFiles(templates), "filesource")
filesource.inputs.subject_id = "subj1"
filesource.outputs.get()


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


# node to compute a mask of region 1 of user atlas
IdentifyRegions = Node(Function(input_names=['atlas_file'],
                                output_names=['regions_images_list'],
                                function=identifyRegionAtlases),
                                name='IdentifyRegions')



connectivity.connect(filesource, "atlasFile", IdentifyRegions, "atlas_file")


# here comes the map node, we have to iterate the treatment for each region identified by previous step
# b = pe.MapNode(interface=B(), name="b", iterfield=['in_file']) 
# http://nipype.readthedocs.io/en/latest/users/mapnode_and_iterables.html

# lets calculate the mean residu signal eg. time courses in every region 


from nipype.interfaces.fsl.utils import ImageMeants

regMeants = MapNode(ImageMeants(), name="regMeants", iterfield=['mask'])     
regMeants.inputs.ignore_exception = False     
regMeants.inputs.order = 1     
regMeants.inputs.output_type = 'NIFTI_GZ'     
regMeants.inputs.terminal_output = 'stream'     
connectivity.connect(filesource, "arimaFile" , regMeants, "in_file")   
connectivity.connect(IdentifyRegions, "regions_images_list", regMeants, "mask")


# data sink
datasink = pe.Node(nio.DataSink(), name='datasink')
datasink.inputs.base_directory = '/scratch/user/hirsch/datadir4/data_results_py'
connectivity.connect(regMeants,  'out_file', datasink, 'functionnal.regMeants')

connectivity.run()


# # fin du pipe

