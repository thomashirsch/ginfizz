
# coding: utf-8
de v1 à v2 on introduit les masks wm et lcf et les meants
de v2 à v3 on plot les mocos, on wrappe le code python, on wrappe le code mathlab bramilla. on refait tout le pipeline
de v3 à v4 on change les premier et dernier noeud en identity interfacepour se brancher direct sur le proprocess pipeline
# # To calculate parameters for afni bandpass function and write them in ortho file.txt
# 

# # Start of the new subworkflow - PREBANDPASS
# rp_at0009_epi_s04_d0001.txt
# todo recuperer le seul fichier texte de results functionnal
# rp_at0009_epi_s+04_d0001.txt

mocoFile = '/homes_unix/hirsch/_pypipe/datadir/data_results/functionnal/rp_at0009_epi_s04_d0001.txt'
# In[1]:

import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.matlab as mlab      # how to run matlab
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine

from nipype.interfaces.utility import Function, IdentityInterface

from nipype.interfaces.fsl.maths import MathsCommand
from nipype.interfaces.fsl.utils import PlotMotionParams   # to plot moco variables


# In[2]:

# creation of a subworflow to calculate the bandpass parameters
prebandpass = pe.Workflow(name='prebandpass')


# In[3]:

# todo remplace this node by an identity node that get input from preprocess pipeline / node 7 eg. fsl merge
# input node get the good files first the tissues normalised files
sourcedir = '/homes_unix/hirsch/_pypipe/datadir/data_results/structural/norm_files'
# second the merged functionnal file
sourcemergeddir = '/homes_unix/hirsch/_pypipe/datadir/data_results/functionnal'


from nipype import SelectFiles, Node
templates = dict(wmMask=sourcedir+ "/" + "wc2*.nii",
                 lcfMask=sourcedir+ "/" + "wc3*.nii",
                 mergedFile=sourcemergeddir+ "/" + "*_merged.nii.gz",
                 mocoVariables=sourcemergeddir+ "/" + "rp*.txt")

filesource = Node(SelectFiles(templates), "filesource")
filesource.inputs.subject_id = "subj1"
filesource.outputs.get()


# ## 1 - compute moco file to feed with ortho.txt file the bandpass node of preprocess workflow 
derivatives
lets assume that the time of acquisition between 2 mesures is TR = 2 000 ms
dx(t)= x(t+1) - x(t-1) / 2 * TR
x(t0) = x(t1)
acqNb = 240


# ## Node 1 - compute moco

# In[4]:


def computeMoco(mocoFile):
    import pandas as pd
    import numpy as np
    import os
    

    # read the moco file to put it in a panda dataframe 
    mocodf = pd.read_csv(mocoFile, header=None, sep='  ',engine='python')
    print(mocodf.head())
    # todo recuperer ces infos de l'autre pipeline
    TR = 2000 
    acqNb = 240

    def vectorDerivative(v):
        dv = {}
        for i in range(acqNb):
            # print mocodf['x'][i]
            if i== 0:
                dv[i]= (v[i+1]-v[i]) / 2*TR
            elif i== acqNb-1:
                dv[i]= (v[i]-v[i-1]) / 2*TR 
            else:
                dv[i]= (v[i+1]-v[i-1]) / 2*TR
                #print 'derivative' + str(i)
                #print  v[i]
        return v

    def plusDerivative(df):
        lg = len(df.columns.values)
        dg = df
        for j in list(df.columns.values):
            vprime = vectorDerivative(df[j])
            dg[lg+j]=vprime
        return dg
    
    def plusSquare(df):
        lg = len(df.columns.values)
        ds = df
        for j in list(df.columns.values):
            vs = df[j]**2
            ds[lg+j]=vs
        return ds    

    
    # first, we derivate the 6 colunms of dataframe of moco file, and append the 6 new colums to df
    dfderivate = plusDerivative(mocodf)
    print(dfderivate.head())
    
    # then, we compute the square of the now 12 colums, and append them to df. it makes 24 colums that are going to 
    # participate in the ortho file to make 24 regressors bandpassed
    dfsquare = plusSquare(dfderivate)
    g = dfsquare.to_csv('ortho.txt', sep=' ', index=False,header=False)
    print g
    h = os.getcwd() + '/' + 'ortho.txt'
    return h
    

computeMoco = Node(Function(input_names=['mocoFile'],
                                output_names=['out_file'],
                                function=computeMoco),
                                name='computeMoco')



prebandpass.connect(filesource, "mocoVariables", computeMoco, "mocoFile")


# ##  2 - Moco plots

# In[5]:

# plot moco variables
MotionCorrectionPlot1 = Node(PlotMotionParams(), name="MotionCorrectionPlot1")
MotionCorrectionPlot1.inputs.ignore_exception = False     
MotionCorrectionPlot1.inputs.in_source = 'spm'     
MotionCorrectionPlot1.inputs.output_type = 'NIFTI_GZ'     
MotionCorrectionPlot1.inputs.plot_size = (500, 1000)     
MotionCorrectionPlot1.inputs.plot_type = 'rotations'     
MotionCorrectionPlot1.inputs.terminal_output = 'stream'     
prebandpass.connect(filesource, "mocoVariables", MotionCorrectionPlot1, "in_file")


# ## 3 - get wm and lcf mask (normalised), erode them and calculate signal mean on both masks
# (input from segment + normalyse wmMask = '/homes_unix/hirsch/_pypipe/datadir/data_results/structural/norm_files/wc2t0009_t1_s03.nii')
# remark: erosion is done once; in connectomics, it is done 3 times

# In[6]:

# calculate eroded binary mask for wm
from nipype.interfaces.fsl.maths import MathsCommand

erosion = pe.Node(interface=MathsCommand(), name='erosion')
    
erosion.inputs.args = '-thr 0 -uthr 111 -bin -ero  '     
erosion.inputs.ignore_exception = False     
erosion.inputs.output_type = 'NIFTI_GZ'     
erosion.inputs.terminal_output = 'stream'     
prebandpass.connect(filesource,"wmMask" , erosion, "in_file")


# In[7]:

# calculate eroded binary mask for lcf

erosionLcf = pe.Node(interface=MathsCommand(), name='erosionLcf')
    
erosionLcf.inputs.args = '-thr 0 -uthr 111 -bin -ero '     
erosionLcf.inputs.ignore_exception = False     
erosionLcf.inputs.output_type = 'NIFTI_GZ'     
erosionLcf.inputs.terminal_output = 'stream'     
prebandpass.connect(filesource,"lcfMask" , erosionLcf, "in_file")


# In[8]:

# lets calculate the mean signal on these eroded masks first on wm
from nipype.interfaces.fsl.utils import ImageMeants

wmMeants = Node(ImageMeants(), name="wmMeants")     
wmMeants.inputs.ignore_exception = False     
wmMeants.inputs.order = 1     
wmMeants.inputs.output_type = 'NIFTI_GZ'     
wmMeants.inputs.terminal_output = 'stream'     
prebandpass.connect(filesource,"mergedFile" , wmMeants, "in_file")   
prebandpass.connect(erosion, "out_file", wmMeants, "mask")


# In[9]:

# lets calculate the mean signal on these eroded masks econd on lcf

lcfMeants = Node(ImageMeants(), name="lcfMeants")     
lcfMeants.inputs.ignore_exception = False     
lcfMeants.inputs.order = 1     
lcfMeants.inputs.output_type = 'NIFTI_GZ'     
wmMeants.inputs.terminal_output = 'stream'     
prebandpass.connect(filesource,"mergedFile" , lcfMeants, "in_file")   
prebandpass.connect(erosionLcf, "out_file", lcfMeants, "mask")


# In[10]:

# lets make an output node with the 3 files to be pass to afni bandpass: 
# file 1 with 24 colums moco, derivative moco, moco and derivative square
# file 2 with wm mean signal on wm
# file 3 with lcf mean signal on lcf

field_list=['rotations_plot', 
            'moco_gradient_square',
            'wm_normalized_eroded_mask',
            'wm_meants',
            'lcf_normalized_eroded_mask',
            'lcf_meants']
            
            
outputNode = Node(IdentityInterface(fields=field_list), name="outputNode")

# for plot files
prebandpass.connect(MotionCorrectionPlot1,  'out_file', outputNode, 'rotations_plot')

# for moco , gradient and square
prebandpass.connect(computeMoco,  'out_file', outputNode, 'moco_gradient_square')

# for segmented normalised eroded wm and lcf mask
prebandpass.connect(erosion,  'out_file', outputNode, 'wm_normalized_eroded_mask')
prebandpass.connect(erosionLcf,  'out_file', outputNode, 'lcf_normalized_eroded_mask')

# for wm and lcf mean signal to text files in functionnal repository
prebandpass.connect(wmMeants,  'out_file', outputNode, 'wm_meants')
prebandpass.connect(lcfMeants,  'out_file', outputNode, 'lcf_meants')


# In[11]:

# the run
prebandpass.run()


# # End of precalculations dump in ortho.txt file

# # Brouillons

# In[ ]:

# data sink
datasink = pe.Node(nio.DataSink(), name='datasink')
datasink.inputs.base_directory = '/homes_unix/hirsch/_pypipe/datadir/data_results'

# for plot files
prebandpass.connect(MotionCorrectionPlot1,  'out_file', datasink, '.plotfiles')

# for segmented normalised eroded wm and lcf mask
prebandpass.connect(erosion,  'out_file', datasink, 'functionnal.bandpass_wm_mask')
prebandpass.connect(erosionLcf,  'out_file', datasink, 'functionnal.bandpass_lcf_mask')

# for wm and lcf mean signal to text files in functionnal repository
prebandpass.connect(wmMeants,  'out_file', datasink, 'functionnal.bandpass_wm_meants')
prebandpass.connect(lcfMeants,  'out_file', datasink, 'functionnal.bandpass_lcf_meants')


# In[ ]:

get_ipython().magic(u'pwd')


# In[ ]:

mocoFile = '/homes_unix/hirsch/_pypipe/datadir/data_results/functionnal/rp_at0009_epi_s04_d0001.txt'


# In[ ]:

computeMoco.inputs.mocoFile=mocoFile
res = computeMoco(mocoFile).run()


# In[ ]:

import os
def mocoOrtho(f):
    df = pd.read_csv(f, header=None, sep='  ',engine='python')
    dfd = plusDerivative(df)
    dfs = plusSquare(dfd)
    g = dfs.to_csv('ortho.txt', sep=' ', index=False,header=False)
    print g
    h = os.getcwd() + '/' + 'ortho.txt'
    return h


# In[ ]:

mocoOrtho('/homes_unix/hirsch/_pypipe/datadir/data_results/functionnal/rp_at0009_epi_s04_d0001.txt')


# In[ ]:

ImageMeants().help()


# In[ ]:

# Node: fsconv.editWmMask     
editWmMask = Node(MathsCommand(), name="editWmMask")     
editWmMask.inputs.args = '-thr 0 -uthr 111 -bin -ero -ero -ero '     
editWmMask.inputs.ignore_exception = False     
editWmMask.inputs.output_type = 'NIFTI_GZ'     
editWmMask.inputs.terminal_output = 'stream'     
fsconv.connect(fsWm2Nii, "out_file", editWmMask, "in_file")

from web
##############################################################################################
##erode a mask or image by zeroing non-zero voxels when zero voxels found in kernel
##############################################################################################
fslmaths 'mask.nii.gz' -kernel box 5x5x5 -ero 'output_image.nii.gz'


# In[ ]:

get_ipython().magic(u"cd '/homes_unix/hirsch/_pypipe/datadir/data_results/functionnal'")
get_ipython().magic(u'pwd')


# In[ ]:

cwd = os.getcwd()
print cwd


# In[ ]:

square = plusSquare(derivate)
print(square.head())


# In[ ]:

vectorSquare(mocodf[0])


# In[ ]:

list(mocodf.columns.values)


# In[ ]:

vectorderivative(mocodf[0])


# In[ ]:

TR = 2000 
acqNb = 240

dx = {}
for i in range(acqNb):
    # print mocodf['x'][i]
    if i== 0:
        dx[i]= (mocodf['x'][i+1]-mocodf['x'][i]) / 2*TR
    elif i== acqNb-1:
        dx[i]= (mocodf['x'][i]-mocodf['x'][i-1]) / 2*TR 
    else:
        dx[i]= (mocodf['x'][i+1]-mocodf['x'][i-1]) / 2*TR
    print 'derivative' + str(i)
    print  dx[i]


# In[ ]:

sdx = pd.Series(dx)


# In[ ]:

sdx


# In[ ]:



