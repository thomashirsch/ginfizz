#!/usr/bin/env python 
# -*- coding: utf-8 -*-



##---------------------------------------------------------------------------------------
#
#        UTILS functions
#
##---------------------------------------------------------------------------------------


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
    
    



##---------------------------------------------------------------------------------------
#
#    end    UTILS functions
#
##---------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------
#
#    Function PREPROCESS 
#
##---------------------------------------------------------------------------------------

def preprocess(spm_standalone, mcr, flibasedir,atlasfile, resultdir):


# ## preliminaries



    import nipype.interfaces.io as nio           # Data i/o
    import nipype.interfaces.spm as spm          # spm
    import nipype.interfaces.matlab as mlab      # how to run matlab
    import nipype.interfaces.utility as util     # utility
    import nipype.pipeline.engine as pe          # pypeline engine
    import nipype.algorithms.rapidart as ra      # artifact detection
    import nipype.algorithms.modelgen as model   # model specification
    import os                                    # system functions


    

    from nipype.interfaces import spm

    # essai du mercredi soir matlab_cmd = ' '.join([spm_standalone,mcr,'batch','script'])
    matlab_cmd = ' '.join([spm_standalone,mcr,'batch'])
    spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)


    # ## Nipype verbosity


    import nipype
    # Optional: Use the following lines to increase verbosity of output
    nipype.config.set('logging', 'workflow_level',  'INFO')
    nipype.config.set('logging', 'interface_level', 'INFO')
    nipype.logging.update_logging(nipype.config)


# ## 1 - we decode the xml file


    xmlfile = getXmlFile(flibasedir)

    print xmlfile


    # In[8]:

    # the we decode the xmlfile
    import xml.etree.ElementTree as ET
    tree = ET.parse(xmlfile)
    root = tree.getroot()


    print root.tag


    epiResults = getEpiInfos(xmlfile)
    print epiResults



    t1 = getT1file(xmlfile)
    print t1


    # # Start Pileline -> declaration


    preproc = pe.Workflow(name='preproc')
    
    # 0 - we fix the working dir to get the commands later
    
    preproc.base_dir = os.path.join(resultdir, "_report")


    # ## 1 - first node data grabbing by select files 

    #os.path.join(dir_name, base_filename + suffix)

    epidir = epiResults['epidir']
    
    abs_epidir = os.path.join(flibasedir, epidir)
    abs_t1 = os.path.join(flibasedir, t1)

    from nipype import SelectFiles, Node
    templates = dict(T1=abs_t1,
                     epi= abs_epidir + "/" + "*.nii")

    filesource = Node(SelectFiles(templates), "filesource")
    filesource.inputs.subject_id = "subj1"
    filesource.outputs.get()


    # ## 2 - segment T1


    segment = pe.Node(interface=spm.NewSegment(), name='segment')

    #seg.inputs.channel_files = 't0009_t1_s03.nii'
    segment.inputs.channel_info =(0.001,60.0,(True,True))          
    segment.inputs.sampling_distance= 3.0
    # todo donner un chemin de TPM.nii comme argument de la future fct
    tissue1 = ((atlasfile , 1), 1, (True,False), (False, False))
    tissue2 = ((atlasfile, 2), 1, (True,False), (False, False))
    tissue3 = ((atlasfile, 3), 2, (True,False), (False, False))
    tissue4 = ((atlasfile, 4), 3, (True,False), (False, False))
    tissue5 = ((atlasfile, 5), 4, (True,False), (False, False))
    tissue6 = ((atlasfile, 6), 2, (False,False), (False, False))
    segment.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
    segment.inputs.affine_regularization = 'mni'
    #segment.inputs.warping_regularization = [0.0, 0.001, 0.5, 0.05, 0.2]
    segment.inputs.write_deformation_fields = [False, False]

    # 2 - segment
    preproc.connect(filesource,"T1" ,segment, 'channel_files')


# ## 2 bis - image calculator



    # def computeStructuralImage(gm_img, wm_img, mean_img):
    #     ''' 
    #     takes 3 images arguments: grey matter, white matter, and bias corrected T1 img.
    #     compute a binary mask with gm and wm thresholded at 0.2 combined with the mean image coming out of segment.
    #     used to realigne properly the EPI imgs. 
    #     '''

    #     import numpy as np     
    #     import nibabel as nib     
    #     import os 
        
    #     img1 = gm_img
    #     img2 = wm_img
    #     img3 = mean_img
        
    #     i1=nib.load(img1)         
    #     i1array=np.asarray(i1.dataobj).copy() # Avoid caching the proxy image
    
    #     i2=nib.load(img2)         
    #     i2array=np.asarray(i2.dataobj).copy()
    #     hdr2 = i2.header
    
    #     i3=nib.load(img3)         
    #     i3array=np.asarray(i3.dataobj).copy()
    #     hdr3 = i3.header
    
    #     print(hdr2.get_data_dtype())
    #     print(hdr3.get_data_dtype())
    
    #     gi = i1array + i2array
    #     # threshold image gm + wm at 0.2
    #     gi[(10 * gi )< 2] = 0
    #     # binary mask the resulting image
    #     gi[gi > 0] = 1
    
    #     struct_image = i3
    
    #     hdr_struct = struct_image.header
    #     print(hdr_struct.get_data_dtype())
    
    #     print"apres"
        
    #     print(hdr_struct.get_data_dtype())
    
    #     # apply the binary mask to the bias corrected T1
    #     struct_image_array = np.multiply(gi, i3array)
   
    #     hdr3.set_data_dtype(np.int16)
    
    #     # crate the resulting image
    #     struct_image = nib.Nifti1Image(struct_image_array.astype(np.int16), i3.affine, i3.header)
        
    #     # save the image to disk
    #     out_file = os.path.join(os.getcwd(), 'struct_image.nii')
    #     nib.save(struct_image, out_file)
            
    #     return out_file

    def regexfilter(files_list,patern):
        import re
        for f in files_list:
            for g in f:
                if re.search(patern, str(g)):
                    res = g
        return res

    
    

    # # In[17]:

    # from nipype.interfaces.utility import Function, IdentityInterface

    # computeStructuralImage = Node(Function(input_names=['gm_img', 'wm_img', 'mean_img'],
    #                                 output_names=['out_file'],
    #                                 function=computeStructuralImage),
    #                                 name='computeStructuralImage')

    # # 2 bis - image calculator -> to ger a nice structurak images to get good registration results
    # # input_names=['gm_img', 'wm_img', 'bias_corrected_img'],
    # preproc.connect(segment, ('native_class_images', regexfilter,r'.*c1.*.nii') ,computeStructuralImage, 'gm_img')
    # preproc.connect(segment, ('native_class_images', regexfilter,r'.*c2.*.nii') ,computeStructuralImage, 'wm_img')
    # # todo retrouver eventuellement la mean image je ne sais pas pourquoi elle n est pas produite par segment
    # preproc.connect(filesource,"T1" ,computeStructuralImage, 'mean_img')

    # new image calculator with fslmaths

    from nipype.interfaces.fsl import MultiImageMaths

    addFiles = pe.Node(interface=MultiImageMaths(), name='addFiles')
    
    addFiles.inputs.op_string = "-add %s"   
    addFiles.inputs.output_datatype = 'short'
    addFiles.inputs.ignore_exception = False     
    addFiles.inputs.output_type = 'NIFTI'     
    addFiles.inputs.terminal_output = 'stream'     

    preproc.connect(segment, ('native_class_images', regexfilter,r'.*c1.*.nii'), addFiles, "in_file")
    preproc.connect(segment, ('native_class_images', regexfilter,r'.*c2.*.nii') , addFiles, "operand_files")

    from nipype.interfaces.fsl import Threshold

    thrFile = pe.Node(interface=Threshold(), name='thrFile')
          
    thrFile.inputs.thresh = 0.2   
    thrFile.inputs.ignore_exception = False     
    thrFile.inputs.output_type = 'NIFTI'     
    thrFile.inputs.terminal_output = 'stream'     

    preproc.connect(addFiles,"out_file" , thrFile, "in_file")



    from nipype.interfaces.fsl import ApplyMask

    maskFiles = pe.Node(interface=ApplyMask(), name='maskFiles')
      
    maskFiles.inputs.ignore_exception = False     
    maskFiles.inputs.output_type = 'NIFTI'     
    maskFiles.inputs.terminal_output = 'stream' 
    maskFiles.inputs.output_datatype ='short'

    preproc.connect(segment,  'bias_corrected_images' , maskFiles, "in_file")
    preproc.connect(thrFile,"out_file" , maskFiles, "mask_file")





    # ## 3 - slice timing of epibold images

  

    l = epiResults['sliceTimingVector'].split()
    lint =[]
    for i in l:
        lint.append(int(i))
    print lint



    from nipype.interfaces.spm import SliceTiming

    st = pe.Node(interface=spm.SliceTiming(), name='st')
    #st.inputs.in_files = 'functional.nii'
    st.inputs.num_slices = int(epiResults['nb_slices'])
    st.inputs.time_repetition = float(epiResults['TR'])
    st.inputs.time_acquisition = float(epiResults['TR']) - float(epiResults['TR'])/float(epiResults['nb_slices'])
    # todo check previous line
    print st.inputs.time_acquisition
    st.inputs.slice_order = lint
    st.inputs.ref_slice = int(epiResults['nb_slices'])
    #st.inputs.jobtype = 'estwrite'

    preproc.connect(filesource,"epi" ,st, 'in_files')


    # ## 4 - realign of epibold images


    realign = pe.Node(interface=spm.Realign(), name='realign')

    realign.inputs.register_to_mean = True
    realign.inputs.quality = 0.9
    realign.inputs.separation = 4
    realign.inputs.fwhm = 5
    realign.inputs.interp = 2
    realign.inputs.wrap = [0, 0, 0]
    # essai au defaut 2, 1 realign.inputs.write_which = [0, 1] je reviens au 2 1 pour dumper
    realign.inputs.write_which = [2, 1]
    realign.inputs.write_mask = True
    realign.inputs.write_interp = 4
    realign.inputs.wrap = [0, 0, 0]
    realign.inputs.write_mask = True
    realign.inputs.out_prefix = 'r'
    realign.inputs.jobtype = 'estwrite'

    preproc.connect(st,"timecorrected_files" ,realign, 'in_files')


    # ## 5 - Coregistration

    coregister = pe.Node(interface=spm.Coregister(), name='coregister' )

    coregister.inputs.cost_function = 'nmi'
    coregister.inputs.separation = [4.0, 2.0]
    coregister.inputs.jobtype = 'estwrite'
    coregister.inputs.tolerance = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
    coregister.inputs.fwhm = [7.0, 7.0]
    coregister.inputs.write_interp = 4
    coregister.inputs.write_wrap = [0, 0, 0]
    coregister.inputs.write_mask = False 
    coregister.inputs.out_prefix = 'c'
    coregister.inputs.jobtype = 'estwrite'
    
    # try  to invert better
    # preproc.connect(computeStructuralImage,"out_file" ,coregister, 'target')
    
    preproc.connect(maskFiles,"out_file" ,coregister, 'target')  
    preproc.connect(realign,"mean_image" ,coregister, 'source')    
    preproc.connect(realign,"realigned_files" ,coregister, 'apply_to_files')


    # ## 6 - Normalize

    # In[22]:

    # passage à normalize sans 12 non on repasse en 12 et on prend des param de pierre y
    norm12 = pe.Node(interface=spm.Normalize12(), name='norm12')

    norm12.inputs.bias_regularization = 0.0001
    norm12.inputs.bias_fwhm = 60
    
    norm12.inputs.tpm = atlasfile
    norm12.inputs.affine_regularization_type = 'mni'
    norm12.inputs.sampling_distance = 3.0
    norm12.inputs.write_bounding_box = [[-90.0, -126.0, -72.0],[ 90.0, 90.0, 108.0]]
    norm12.inputs.write_voxel_sizes = [2.0, 2.0, 2.0]
    norm12.inputs.write_interp = 4 
    norm12.inputs.jobtype = 'estwrite'
    norm12.inputs.use_v8struct = True
    
    # reg = [0 0.001 0.5 0.05 0.2];
    norm12.inputs.warping_regularization = [0.0, 0.001, 0.5, 0.05, 0.2]
    # fwhm = 0;  i don t known how this parameter is called in nipype 
    
    #norm12.inputs.out_prefix = 'w'

    preproc.connect(segment,'bias_corrected_images' ,norm12, 'image_to_align')
    preproc.connect(coregister,"coregistered_files" ,norm12, 'apply_to_files')


    # deuxieme normalisation pour normaliser les masks wm et lcf
    norm12bis = pe.Node(interface=spm.Normalize12(), name='norm12bis')

    norm12bis.inputs.bias_regularization = 0.0001
    norm12bis.inputs.bias_fwhm = 60
    
    norm12bis.inputs.tpm = atlasfile
    norm12bis.inputs.affine_regularization_type = 'mni'
    norm12bis.inputs.sampling_distance = 3.0
    norm12bis.inputs.write_bounding_box = [[-90.0, -126.0, -72.0],[ 90.0, 90.0, 108.0]]
    norm12bis.inputs.write_voxel_sizes = [2.0, 2.0, 2.0]
    norm12bis.inputs.write_interp = 4 
    norm12bis.inputs.jobtype = 'estwrite'
    norm12bis.inputs.use_v8struct = True
    
    # reg = [0 0.001 0.5 0.05 0.2];
    norm12bis.inputs.warping_regularization = [0.0, 0.001, 0.5, 0.05, 0.2]
    # fwhm = 0;  i don t known how this parameter is called in nipype 

    # rajout de normalisation des tissus
    preproc.connect(segment,'bias_corrected_images' ,norm12bis, 'image_to_align')
    preproc.connect(segment,"native_class_images" ,norm12bis, 'apply_to_files')


    # In[24]:

    # data sink provisoire
    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = resultdir

    # segment results
    preproc.connect(segment,  'normalized_class_images', datasink, 'structural')
    preproc.connect(segment,  'bias_corrected_images' , datasink, 'structural.bias_corrected_images')
    preproc.connect(segment,  'bias_field_images' , datasink, 'structural.bias_field_images')
    preproc.connect(segment,  'forward_deformation_field' , datasink, 'structural.forward_deformation_field')
    preproc.connect(segment,  'inverse_deformation_field' , datasink, 'structural.inverse_deformation_field')
    preproc.connect(segment,  'modulated_class_images' , datasink, 'structural.mod_class')
    preproc.connect(segment,  'native_class_images' , datasink, 'structural.native_class')
    preproc.connect(segment,  'transformation_mat' , datasink, 'structural.transformation_mat')
    preproc.connect(segment,  'dartel_input_images' , datasink, 'structural.dartel_input_images')

    # compute image preproc.connect(computeStructuralImage,"out_file" ,coregister, 'source')
    preproc.connect(maskFiles, "out_file", datasink, 'structural.t1_masked_files')

    # rajout de normalisation des tissues
    preproc.connect(norm12bis,  'normalized_files', datasink, 'structural.normalized_files')
    preproc.connect(norm12bis,  'normalized_image', datasink, 'structural.normalized_image')

    # slice timing results
    preproc.connect(st,  'timecorrected_files', datasink, 'functionnal')

    # realign mean_image and realignment_parameters
    preproc.connect(realign,  'mean_image', datasink, 'functionnal.mean_image')
    preproc.connect(realign,  'realignment_parameters', datasink, 'functionnal.realignment_parameters')
    preproc.connect(realign,  'realigned_files', datasink, 'functionnal.realigned_files')
    
    # coregistration forgot coregistered_source
    preproc.connect(coregister,  'coregistered_files', datasink, 'functionnal.coregistered_files')
    preproc.connect(coregister,  'coregistered_source', datasink, 'functionnal.coregistered_source')

    preproc.connect(norm12,  'normalized_files', datasink, 'functionnal.normalized_files')
    preproc.connect(norm12,  'normalized_image', datasink, 'functionnal.normalized_image')


    # In[25]:

    preproc.run()

    # visualisation workflow
    preproc.write_graph(graph2use='colored', format='svg', simple_form=True)

    return
    # The end


#def main(name, spm_standalone, mcr, flibasedir, resultdir):
#    preprocess(name, spm_standalone, mcr, flibasedir, atlasfile, resultdir)
#    return




if __name__ == "__main__":
    import sys 
    print sys.argv
    preprocess(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])