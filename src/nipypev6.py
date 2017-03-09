
# coding: utf-8

# plus slice timing

# ## preliminaries

# In[1]:

import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.matlab as mlab      # how to run matlab
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.algorithms.rapidart as ra      # artifact detection
import nipype.algorithms.modelgen as model   # model specification
import os                                    # system functions


# In[2]:

spm_standalone='/homes_unix/hirsch/historique_fli_iam/essai_spm_stand_alone/spm12/run_spm12.sh'
mcr='/homes_unix/hirsch/historique_fli_iam/essai_spm_stand_alone/mcr2016/v91'


# In[3]:

from nipype.interfaces import spm

# essai du mercredi soir matlab_cmd = ' '.join([spm_standalone,mcr,'batch','script'])
matlab_cmd = ' '.join([spm_standalone,mcr,'batch'])
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)


# ## Nipype verbosity

# In[4]:

import nipype
# Optional: Use the following lines to increase verbosity of output
nipype.config.set('logging', 'workflow_level',  'INFO')
nipype.config.set('logging', 'interface_level', 'INFO')
nipype.logging.update_logging(nipype.config)


# ## 1 - we decode the xml file

# In[5]:

flibasedir = '/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01'


# In[6]:

# we get the xml file in this directory
def getXmlFile(flibasedir):
    import os
    dirs = os.listdir( flibasedir )
    #print dirs

    import re
    xmlmatch = r'(.*)xml'
    for dir in dirs:
        xmlfound=re.search(xmlmatch,dir,flags=0)
        if xmlfound:
            #print xmlfound.group()
            return flibasedir + '/' + xmlfound.group()


# In[7]:

xmlfile = getXmlFile(flibasedir)

print xmlfile


# In[8]:

# the we decode the xmlfile
import xml.etree.ElementTree as ET
tree = ET.parse(xmlfile)
root = tree.getroot()


# In[9]:

print root.tag

def getEpiInfos(xmlfile):
    result = {}
    try:
        import xml.etree.ElementTree as ET
        tree = ET.parse(xmlfile)
        root = tree.getroot()
        # epibold directory
        for epibold in root.iter('EPIBOLD'):
            for child in epibold:
                print child.text
                if child.tag == 'file':
                    epidir = child.text
                    result['epidir'] = epidir
                    print epidir
                    
                # parameters
                if child.tag == 'parameters':
                    print "parameters"
                    
                    for child2 in child:
                                                     
                        # TR
                        if child2.tag == 'TR':
                            for child3 in child2:
                                if child3.tag == 'value':
                                    tr = child3.text
                                    print tr
                                    result['TR'] = tr
                    
                        # dynamics
                        if child2.tag == 'dynamics':
                            for child3 in child2:
                                if child3.tag == 'value':
                                    dynamics = child3.text
                                    result['dynamics'] = dynamics
                                    
                        # sliceTimingVector
                        if child2.tag == 'sliceTimingVector':
                            for child3 in child2:
                                if child3.tag == 'value':
                                    sliceTimingVector = child3.text
                                    result['sliceTimingVector'] = sliceTimingVector
                                    
                        # nb_slices
                        if child2.tag == 'nb_slices':
                            for child3 in child2:
                                if child3.tag == 'value':
                                    nb_slices = child3.text
                                    result['nb_slices'] = nb_slices
                
        
    except:
        print 'exception'
    return result


# In[10]:

epiResults = getEpiInfos(xmlfile)
print epiResults


# In[11]:

def getT1file(xmlfile):
    try:
        import xml.etree.ElementTree as ET
        tree = ET.parse(xmlfile)
        root = tree.getroot()
        # T1 directory and file
        for t1 in root.iter('T1'):
            for child in t1:
                #print child.text
                if child.tag == 'file':
                    t1 = child.text
                    t1 = flibasedir + t1
                    # print  t1
    except:
        print 'exception'
    return t1
    
    


# In[12]:

t1 = getT1file(xmlfile)
print t1


# # Pileline declaration

# In[65]:

preproc = pe.Workflow(name='preproc')


# ## 1 - first node data grabbing by select files 

# In[66]:

epidir = getEpiInfos(xmlfile)['epidir']
epidir = flibasedir + epidir

from nipype import SelectFiles, Node
templates = dict(T1=t1,
                 epi= epidir + "/" + "*.nii")

filesource = Node(SelectFiles(templates), "filesource")
filesource.inputs.subject_id = "subj1"
filesource.outputs.get()


# ## 2 - segment T1

# In[67]:

segment = pe.Node(interface=spm.NewSegment(), name='segment')

#seg.inputs.channel_files = 't0009_t1_s03.nii'
segment.inputs.channel_info =(0.001,60.0,(True,True))          
segment.inputs.sampling_distance= 3.0
# todo donner un chemin de TPM.nii comme argument de la future fct
tissue1 = (('/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM.nii', 1), 1, (True,False), (False, False))
tissue2 = (('/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM.nii', 2), 1, (True,False), (False, False))
tissue3 = (('/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM.nii', 3), 2, (True,False), (False, False))
tissue4 = (('/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM.nii', 4), 3, (True,False), (False, False))
tissue5 = (('/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM.nii', 5), 4, (True,False), (False, False))
tissue6 = (('/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM.nii', 6), 2, (False,False), (False, False))
segment.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
segment.inputs.affine_regularization = 'mni'
#segment.inputs.warping_regularization = [0.0, 0.001, 0.5, 0.05, 0.2]
segment.inputs.write_deformation_fields = [False, True]


# ## 3 - slice timing of epibold images

# In[68]:

l = epiResults['sliceTimingVector'].split()
lint =[]
for i in l:
    lint.append(int(i))
print lint


# In[69]:

from nipype.interfaces.spm import SliceTiming

st = pe.Node(interface=spm.SliceTiming(), name='st')
#st.inputs.in_files = 'functional.nii'
st.inputs.num_slices = int(epiResults['nb_slices'])
st.inputs.time_repetition = float(epiResults['TR'])
st.inputs.time_acquisition = float(epiResults['TR']) - float(epiResults['TR'])/2.0
st.inputs.slice_order = lint
st.inputs.ref_slice = 1


# ## 4 - realign of epibold images

# In[70]:

realign = pe.Node(interface=spm.Realign(), name='realign')

realign.inputs.register_to_mean = False
realign.inputs.quality = 0.9
realign.inputs.separation = 4
realign.inputs.fwhm = 5
realign.inputs.interp = 2
realign.inputs.wrap = [0, 0, 0]
# essai au defaut 2, 1 realign.inputs.write_which = [0, 1]
realign.inputs.write_which = [2, 1]
realign.inputs.write_interp = 4
realign.inputs.wrap = [0, 0, 0]
realign.inputs.write_mask = True
realign.inputs.out_prefix = 'r'
realign.inputs.jobtype = 'write'


# ## 5 - Coregistration

# In[71]:

# spm.Coregister()
coregister = pe.Node(interface=spm.Coregister(), name='coregister' )

coregister.inputs.cost_function = 'nmi'
coregister.inputs.separation = [4.0, 2.0]
coregister.inputs.jobtype = 'write'
coregister.inputs.tolerance = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
coregister.inputs.fwhm = [7.0, 7.0]
coregister.inputs.write_interp = 4
coregister.inputs.write_wrap = [0, 0, 0]
coregister.inputs.write_mask = False 
coregister.inputs.out_prefix = 'c'


# ## 6 - Normalize

# In[72]:

# passage à normalize sans 12 non on repasse en 12 et on prend des param de pierre y
norm12 = pe.Node(interface=spm.Normalize12(), name='norm12')

norm12.inputs.bias_regularization = 0.0001
norm12.inputs.bias_fwhm = 60
# todo putthis parameter in args of a fct
# norm12.inputs.tpm = '/homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM.nii'
norm12.inputs.affine_regularization_type = 'mni'
norm12.inputs.sampling_distance = 3.0
norm12.inputs.write_bounding_box = [[-90.0, -126.0, -72.0],[ 90.0, 90.0, 108.0]]
norm12.inputs.write_voxel_sizes = [2.0, 2.0, 2.0]
norm12.inputs.write_interp = 4
norm12.inputs.jobtype = 'write'
norm12.inputs.use_v8struct = True
#norm12.inputs.out_prefix = 'w'


# # connection du pipeline

# In[73]:

preproc.connect(filesource,"T1" ,segment, 'channel_files')

preproc.connect(filesource,"epi" ,st, 'in_files')

preproc.connect(st,"timecorrected_files" ,realign, 'in_files')

preproc.connect(filesource,"T1" ,coregister, 'source')
preproc.connect(realign,"mean_image" ,coregister, 'target')
preproc.connect(realign,"realigned_files" ,coregister, 'apply_to_files')

preproc.connect(segment,"forward_deformation_field" ,norm12, "deformation_file")
preproc.connect(coregister,"coregistered_files" ,norm12, 'apply_to_files')


# ## Put the results in the datasink

# In[74]:

datasink = pe.Node(nio.DataSink(), name='datasink')
datasink.inputs.base_directory = '/homes_unix/hirsch/_pypipe/datadir/data_results'


# structural results

# In[75]:

# segment results
preproc.connect(segment,  'normalized_class_images', datasink, 'structural')
preproc.connect(segment,  'bias_corrected_images' , datasink, 'structural.@par')
preproc.connect(segment,  'bias_field_images' , datasink, 'structural.@par1')
preproc.connect(segment,  'forward_deformation_field' , datasink, 'structural.@par2')
preproc.connect(segment,  'inverse_deformation_field' , datasink, 'structural.@par3')
preproc.connect(segment,  'modulated_class_images' , datasink, 'structural.@par4')
preproc.connect(segment,  'native_class_images' , datasink, 'structural.@par5')
preproc.connect(segment,  'transformation_mat' , datasink, 'structural.@par6')
preproc.connect(segment,  'dartel_input_images' , datasink, 'structural.@par7')


# functional results

# In[76]:

# slice timing results
preproc.connect(st,  'timecorrected_files', datasink, 'functionnal')

# realign mean_image and realignment_parameters
preproc.connect(realign,  'mean_image', datasink, 'functionnal.@par')
preproc.connect(realign,  'realignment_parameters', datasink, 'functionnal.@par1')
preproc.connect(realign,  'realigned_files', datasink, 'functionnal.@par2')
preproc.connect(coregister,  'coregistered_files', datasink, 'functionnal.@par3')

preproc.connect(norm12,  'normalized_files', datasink, 'functionnal.@par4')
preproc.connect(norm12,  'normalized_image', datasink, 'functionnal.@par5')


# # THE RUN

# In[77]:

preproc.run()


# ## the end

# # Brouillons

# In[26]:

spm.Normalize12().help()


# In[ ]:

spm.Coregister.help()


# In[ ]:

spm.Realign().help()


# In[ ]:

spm.Realign().input_spec()


# In[ ]:

l = epiResults['sliceTimingVector'].split()
lint =[]
for i in l:
    lint.append(int(i))
print lint


# In[ ]:

float(epiResults['TR']) - float(epiResults['TR'])/2.0


# In[ ]:

spm.SliceTiming().help()


# In[ ]:

spm.SliceTiming().input_spec()


# In[ ]:

bias_corrected_images = <undefined>
bias_field_images = <undefined>
dartel_input_images = <undefined>
forward_deformation_field = <undefined>
inverse_deformation_field = <undefined>
modulated_class_images = <undefined>
native_class_images = <undefined>
normalized_class_images = <undefined>
transformation_mat = <undefined>


# In[ ]:

#advanced way to connect multiple nodes
workflowname.connect([(nodename1, nodename2, [('output_node1', 'input_node2')]),
                      (nodename1, nodename3, [('output_node1', 'input1_node3')]),
                      (nodename2, nodename3, [('output1_node2', 'input1_node3'),
                                              ('output2_node2', 'input2_node3')
                                              ])
                      ])


# In[ ]:




# In[ ]:

preproc.connect([(filesource, segment, [('T1', 'channel_files')]),
                      (segment, datasink, [('normalized_class_images', 'normalized'),
                                              ('bias_corrected_images', 'bias_corrected')
                                              ])
                      ])


# In[ ]:

datasink = pe.Node(nio.DataSink(), name='datasink')
datasink.inputs.base_directory = '/homes_unix/hirsch/_pypipe/datadir/data_results'
preproc.connect(segment,  'normalized_class_images', datasink, 'structural')
preproc.connect(segment,  'bias_corrected_images' , datasink, 'structural.@par')


# In[ ]:




# In[ ]:

preproc.run()


# In[ ]:

spm.NewSegment.help()


# In[ ]:

spm.NewSegment.input_spec()


# In[ ]:

spm.NewSegment.output_spec()


# In[ ]:

homes_unix(/hirsch/_pypipe/datadir/data_set/t0009/repos01/Atlases/TPM)


# In[ ]:

#Where should the output data be stored at?
sink = nipype.DataSink()
sink.inputs.base_directory = experiment_dir + '/output_folder'


# In[ ]:

get_ipython().magic(u'cd /homes_unix/hirsch/_pypipe/datadir/data_set/t0009/repos01/EPIBOLD')
get_ipython().magic(u'pwd')

