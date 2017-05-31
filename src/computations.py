# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 17:30:27 2015
Functions for the REST pipeline
@author: herve
"""

# formatting functions

def bigMerge(in_files,tr,force_even,sort=False):
    """
    Equivalent to fslmerge -t
    This functions can sort the in_files list alphabetically if needed.
    BUT you will not want  to sort if your files are not zero-padded.
    Maybe wise to reserve memory in the qsub/sbatch args.
    """    
    import numpy as np
    import nibabel as nb
    import os 
    
    if type(in_files)==list:
        if sort==True:
            in_files.sort()
        first=nb.load(in_files[0])
        bignifti=np.asarray(first.dataobj).copy() # Avoid caching the proxy image, thanks G. V. !!
        print(in_files[0])
        del first
        if len(np.shape(bignifti))==3:
            bignifti=bignifti[...,np.newaxis]
        for in_file in in_files[1:]:
            print(in_file)
            nft=nb.load(in_file)
            ima=np.asarray(nft.dataobj).copy()
            if len(np.shape(ima))==3:
                ima=ima[...,np.newaxis]                # because otherwise, awful memory leaks occur
            bignifti=np.concatenate((bignifti,ima), 3) # concatenate 4D volumes
            del ima
            del nft  # otherwise massive memory leak.
        n_slices=np.shape(bignifti)[2]
        if ((np.mod(n_slices,2) != 0) & (force_even==True)):
            bignifti=bignifti[:,:,range(n_slices-1),:]  # discard last slice if odd (for topup subsampling schemes)
        nft=nb.load(in_files[0])
        niout= nb.Nifti1Image(bignifti, nft.get_affine())
        niout.header['pixdim'][4]=float(tr)
        keys=['scl_inter', 'slice_end', 'slice_code','qform_code', 'sform_code','intent_p1', 
                   'intent_p2', 'intent_p3', 'intent_code', 'slice_end', 'slice_code', 'xyzt_units', 
                   'cal_max', 'cal_min', 'slice_duration','toffset', 'glmax', 'glmin']
        for key in keys:
            niout.header[key]=nft.header[key]
        fname=os.path.abspath('./'+os.path.split(in_files[0])[1].split('.')[0]+'_merged.nii.gz')
        niout.to_filename(fname)
        del nft
        del niout
    else:
        fname=in_files
        print("ERROR: can't output a 4D file from a single 3D file. Please supply a list !! Now just forwarding the input file.")
    return fname


    
# get the tr from a nifti file
def get_tr_info(nifti_file):
    
    import nibabel as nb
    func4D = nb.load(nifti_file)
    tr=func4D.header['pixdim'][4]
    return tr
    

# DICOM parsing utilities
def get_func_info(dcm_file):
    '''
    Extracts slice timing info from SIEMENS proprietary DICOM field. 
    Writes down a text file with the slice times in ms. Also returns a few useful values by the same token.
    This is intended for EPI-BOLD time series.    
    
    Parameters
    ----------
        dcm_file: sample dicom file
    
    Outputs
    -------
        out_file: path to output text file with slice times in seconds
        tr: repetition-time in seconds
        dwelltime: the EFFECTIVE echo-spacing, regardless of PAT factor
        PEinfo: a string describing the axis and direction of phase-encoding
    
    '''
    import dicom
    import numpy as np
    import os 
    import subprocess
    import re
    
    # Siemens only for now
    # Will have to account for GE and Philips in the future, but how ?
    dcm = dicom.read_file(dcm_file)
    try: 
        slice_times=dcm[0x0019, 0x1029][:]
        st=np.array(slice_times)/1000.0  # milliseconds to seconds, 3dTShift wants same unit for TR and slice times.
        out_file=os.path.abspath('slice_times.txt')
        np.savetxt(out_file,np.atleast_2d(st),'%5.5f')
    except:
        out_file='not created, DICOM field 0019,1029 not present...not a mosaic dicom!!'
    # other useful information    
    tr = dcm[0x018, 0x0080].value/1000.0  # ms to seconds...
    dwelltime = 1.0/(dcm[0x0019, 0x1028].value*float(re.findall('(\d*)',dcm[0x0051,0x100b].value.split('*')[0])[0]))  # haaa...dicom....
    csa_dump=subprocess.Popen(['gdcmdump','--csa', dcm_file], stdout=subprocess.PIPE).stdout.read().split('\n') #haaa...siemens shadow header...
    for chunk in csa_dump:
        if re.match('.*(PhaseEncodingDirectionPositive).*',chunk):
            peDirPositive= int(re.findall(".*'(.*)'",chunk)[0])
            # now determine the axis, whether AP,LR,IS...
            acqPlane=dcm[0x0051, 0x100e].value[:3]
            phaseRC=dcm[0x0018,0x1312].value
            if (acqPlane=='Tra') & (phaseRC=='COL'):
                PEinfo='y'
            elif (acqPlane=='Tra') & (phaseRC=='ROW'):
                PEinfo='x'
            elif (acqPlane=='Cor') & (phaseRC=='COL'):
                PEinfo='z'
            elif (acqPlane=='Cor') & (phaseRC=='ROW'):
                PEinfo='x'
            elif (acqPlane=='Sag') & (phaseRC=='COL'):
                PEinfo='z' 
            elif (acqPlane=='Sag') & (phaseRC=='ROW'):
                PEinfo='y'
            if peDirPositive == 1:
                PEinfo='-'+PEinfo 
            del peDirPositive   

    return (out_file, tr, dwelltime, PEinfo)    
    
def identify_b0s(dcm_files,thr=6):
    '''
    From a set a of dicom files, this function will identify the series or 
    mosaic frames with a b-value equal to zero and output two ordered lists 
    of lists of dicom files for the b0 and actual dwi series.
    This function is SIEMENS based. 
    
    Arguments:
    
    dcm_files:    A list of dicom files
    thr      :    The threshold for b-value=0 files (the DWI savvy
                  folks don't specify b=0, rather b0=5 because the scanners
                  can't reach a b-value of 0 anyway)
                  
    Outputs:
    
    b0_images: b0 DICOM files, in serial order
    dwi_series: dwi DICOM files, in serial order
    '''
    
    import dicom    

    bvalues=dict()
    series=dict()
    b0_series=dict()
    print('scanning DICOM instances')
    for dcm_file in dcm_files:
        dcm = dicom.read_file(dcm_file)
        try:
            bval=int(dcm[(0x0019,0x100c)].value)
        except:
            # if field is absent, SE-EPI and not quite a dwi b0 file
            bval=0
            print('No b-value field!')
        series_number=dcm[(0x0020,0x0011)].value
        inst_number=dcm[(0x0020,0x0013)].value
        # build dictionaries for b=0 and b>0, per serie and instance number
        if bval<=thr:
            if series_number in b0_series.keys():
                b0_series[series_number][inst_number]=dcm_file   
            else:
                b0_series[series_number]=dict()
                b0_series[series_number][inst_number]=dcm_file
        else:
            # That's a real diffusion weighted image
            if bval in bvalues.keys():
                bvalues[bval].append(dcm_file)
            else:
                bvalues[bval]=[dcm_file]
            if series_number in series.keys():
                series[series_number][inst_number]=dcm_file  
            else:
                series[series_number]=dict()
                series[series_number][inst_number]=dcm_file       
    
    print('reinjecting b0 instances in  dwi series')
    # reinject mosaic b0 in mosaic series
    for key in b0_series.keys():
        if key in series.keys():
            for instance in b0_series[key].keys():
                series[key][instance]=b0_series[key][instance]
    
    # now sort per series and instance number
    # list of lists for dwi mosaic series including b0s
    # list of lists for all b0s (individual 3D b0s + mosaic b0s)
    print('Sorting dwi series')  
    dwi_series=[]
    keys=map(int,series.keys())
    keys.sort()
    for key in keys:
            keys2=map(int,series[key].keys())
            keys2.sort()
            dwi_series.append([series[key][key2] for key2 in keys2])
    
    print('Sorting b0 series')               
    b0_images=[]
    keys=map(int,b0_series.keys())
    keys.sort()
    for key in keys:
        keys2=map(int,b0_series[key].keys())
        keys2.sort()
        b0_images.append([b0_series[key][key2] for key2 in keys2])   
    print('FINISHED !')
         
    return (b0_images, dwi_series)
        
def get_series(dcm_files,run_name):
    '''
    Parses all DICOM instances and selects those
    whose series description match 'run_name'
    Outputs a list of paths to DICOM instances, 
    in the serial order.
    
    run_name can also be a tupple of the form ('name',number)
    where number is the nth occurrence of the series in the study.
    
    Intended primarily for functional neuroimaging runs
    '''    
    import dicom 
    
    # build a dictionary of the series and instances 
    # with the correct 'run_name' 
    
    fMRI_series=dict()
    idx = 0
    if type(run_name) is tuple:
        idx = run_name[1]
        run_name = run_name[0]
   
    for dcm_file in dcm_files:
        dcm = dicom.read_file(dcm_file)
        series_description=dcm[(0x0008,0x103E)].value
        if series_description==run_name:
            series_number=dcm[(0x0020,0x0011)].value
            inst_number=dcm[(0x0020,0x0013)].value
            if series_number in fMRI_series.keys():
                fMRI_series[series_number][inst_number]=dcm_file   
            else:
                fMRI_series[series_number]=dict()
                fMRI_series[series_number][inst_number]=dcm_file
    
    # making an ordered list of instances, sorted per series and instance numbers
    # if no index is supplied we take the first matching series

    fMRI_runs=[]
    keys=map(int,fMRI_series.keys())
    keys.sort()
    for key in keys:
            keys2=map(int,fMRI_series[key].keys())
            keys2.sort()
            fMRI_runs.append([fMRI_series[key][key2] for key2 in keys2])
    return fMRI_runs[idx]
    
def get_topup_parameters(dcm_files,nifti_files,output_file):
    """Extracts the parameters necessary for topup from the SIEMENS DICOM headers
    We have to rely on gdcmdump at some point (through the shell as I could not 
    get the python wrappers working yet).
    This is intended for pairs of b-value=0 images from DWI series.
    
    Parameters
    ----------
    
    dcm_files   : a list of dcm filenames, from all the series
    nifti_files : the corresponding NIfTI files as obtained with dcm2nii 
                  with the series number in the filename (this is really important!)
    output_filename : file name for the output
    
    Outputs:
    --------
    parameter_file
    b0names     : name of the series in the same order as b0
    
    """
    import dicom as dcm
    import numpy as np
    import os
    import re    
    import subprocess  
    import nibabel as nb
    
    # match nifti file with dcm series
    series=list()
    b0names=list()
    encdir_vecs=list()
    for nifti_file in nifti_files:
        nifti_vol=nb.load(nifti_file)
        dims= nifti_vol.header.get_data_shape()
        if len(dims)==3:
            series.append([int(re.findall('s(.*)a',nifti_file.split('.')[0][-9:])[0])])
        else:
            series.append([int(re.findall('s(.*)a',nifti_file.split('.')[0][-9:])[0]) for i in range(dims[3])])
    for (sidx,serie) in enumerate(series):
        for (fidx,frame) in enumerate(serie):
            dcm_file=dcm_files[sidx][fidx]
            bv0 = dcm.read_file(dcm_file)
            b0names.append(str(bv0[(0x0020,0x0011)].value))
            ETlength=float(re.findall('(\d*)',bv0[0x0051,0x100b].value.split('*')[0])[0]) # take the first field of the matrix size, remove any attached letter.
            total_readout=1.0/(bv0[0x0019, 0x1028].value*ETlength)*(ETlength-1.0)
            #phase encoding direction (positive or negative)
            #Well, what matters here is that we get the same convention for all files entered in topup or epi_reg
            # we want to parse the Siemens ASCII Shadow Header
            # This information is apparently absent with GE scanners
            # Don't know about Philips/bruker
            csa_dump=subprocess.Popen(['gdcmdump','--csa', dcm_file], stdout=subprocess.PIPE).stdout.read().split('\n')
            for chunk in csa_dump:
                if re.match('.*(PhaseEncodingDirectionPositive).*',chunk):
                    peDirPositive= int(re.findall(".*'(.*)'",chunk)[0])
                    # now determine the axis, whether AP,LR,IS...
                    acqPlane=bv0[0x0051, 0x100e].value[:3]
                    phaseRC=bv0[0x0018,0x1312].value
                    if (acqPlane=='Tra') & (phaseRC=='COL'):
                        orVec=np.array([[0, 1, 0]])
                    elif (acqPlane=='Tra') & (phaseRC=='ROW'):
                        orVec=np.array([[1, 0, 0]])
                    elif (acqPlane=='Cor') & (phaseRC=='COL'):
                        orVec=np.array([[0, 0, 1]])   
                    elif (acqPlane=='Cor') & (phaseRC=='ROW'):
                        orVec=np.array([[1, 0, 0]])
                    elif (acqPlane=='Sag') & (phaseRC=='COL'):
                        orVec=np.array([[0, 0, 1]])   
                    elif (acqPlane=='Sag') & (phaseRC=='ROW'):
                        orVec=np.array([[0, 1, 0]]) 
                    if (peDirPositive == 1):
                        orVec=orVec * -1  
                    # stack-up the values for this image...
                    encdir_vecs.append(np.concatenate((orVec,np.array([[total_readout]])),axis=1))
                    del orVec
                    del peDirPositive
    # build and write parameter table
    parameters = np.concatenate(encdir_vecs,0)
    np.savetxt(output_file,parameters,fmt='%1.8f')
    parameter_file=os.path.abspath(output_file)
                       
    return (parameter_file,b0names)
    
def check_fix_nifti(ima_file):
    """Corrects the sform/qform matrices if translations set to 0 
    because SPM really hates it when the volumes are not properly centered.
    Correction consists in setting the origin to the center of the volume, 
    using the built-in nibabel functionality.

    Parameters
    ----------

    ima_file: a NIfTI file name (str)

    Returns
    -------
    status: [sform, qform] False if nothing done, True if matrix fixed.
        
    """
    import nibabel as nb
    import numpy as np
    
    ima=nb.load(ima_file)
    status=[False]*2
    if np.array(ima.get_sform()[0:3,3] == [ 0.,  0.,  0.]).all():
        # needs a fix
        status[0] = True
        sform = np.array(ima.get_sform())
        sform = ima.header.get_base_affine()
        ima.set_sform(sform)
    else:
        status[0] = False
    if np.array(ima.get_qform()[0:3,3] == [ 0.,  0.,  0.]).all():
        # image needs a fix
        status[1] = True
        qform = np.array(ima.get_qform())
        qform = ima.header.get_base_affine()
        ima.set_qform(qform)
    else:
        status[1] = False
    if np.array(status).any():
        # we overwrite the image !
        nb.save(ima,ima_file)
    return status    


def fix_nifti_origin(in_file,out_filename):
    '''
    Applies nifti convention regarding the origin of the referential
    
    in_file:   input file
    out_filename:   output file
    
    Returns the absolute path to the output file
    '''
    
    import nibabel as nb
    import numpy as np    
    import os
    
    ima=nb.load(in_file)
    
    # Get sform and qform
    sform = ima.header.get_sform()
    qform = ima.header.get_qform()
    
    # Compute coordinates of middle voxel (keep subvoxel precision)
    midvox = np.concatenate(((np.array(ima.get_shape())-1)/2.0,[1.0]))
    midvoxcoords= np.dot(sform,midvox)

    # Recenter sform referential    
    nsform  = sform
    nsform[:3,3]=sform[:3,3]-midvoxcoords[:3]
    ima.set_sform(nsform)
    
    # Exact same thing with qform (with the same algebra)
    nqform  = qform
    nqform[:3,3]=qform[:3,3]-midvoxcoords[:3]
    ima.set_qform(qform)  
    
    # Write out    
    ima.to_filename(out_filename)
    out_file = os.path.abspath(out_filename)
    
    return out_file
    
def fix_crop_origin(in_file,coords,out_filename):
    '''
    Useful for fixing scilpy cropped volumes.
    
    in_file:   input file
    out_filename:   output file
    
    Returns the absolute path to the output file
    DEPRECATED AS OF MAY 2016
    '''
    
    import nibabel as nb
    import numpy as np    
    import os
    import pickle
    
    ima=nb.load(in_file)
    coords_file=open(coords,'r')
    crop=pickle.load(coords_file)
    coords_file.close()
    
    # Get sform and qform
    sform = ima.header.get_sform()
    qform = ima.header.get_qform()
    
    # compute the translation of the referential
    zero_vox= np.dot(np.linalg.inv(ima.get_sform()),[0.0,0.0,0.0,1.0])
    
    new_zerovox=zero_vox[0:3] - crop[0]
    new_zerovox=np.concatenate((new_zerovox,[1.0]))
    new_zerocoords= np.dot(sform,new_zerovox)

    # Recenter sform referential    
    nsform  = sform
    nsform[:3,3]=sform[:3,3]-new_zerocoords[:3]
    ima.set_sform(nsform)
    
    # Exact same thing with qform (with the same algebra...???)
    nqform  = qform
    nqform[:3,3]=qform[:3,3]-new_zerocoords[:3]
    ima.set_qform(qform)  
    
    # Write out    
    ima.to_filename(out_filename)
    out_file = os.path.abspath(out_filename)
    
    return out_file     
    
def topup_fix_orient(corrected,field,fieldcoef,orig,movpar):
    '''
    This function is meant to fix a problematic loss of orientation (sform/qform)
    in topup results with FSl 5.0.8 (origin is lost, as well as LR info).
    The best thing to do is probably to copy the orientation from the b0 before 
    correction.
    Now, the question is: what about the motion parameters ?
    For now we just copy them in the fixed folder for eddy to find them.
    
    NOT USEFUL ANYMORE WITH FSL 5.0.9
    '''
    
    import os
    import nibabel as nb
    import subprocess
    
    iorig = nb.load(orig)
    icorr= nb.load(corrected)
    ifield = nb.load(field)
    ifieldc= nb.load(fieldcoef)
    
    sform=iorig.get_sform()
    scode=iorig.header['sform_code']
    qform=iorig.get_qform()
    qcode=iorig.header['qform_code']
    
    icorr.set_sform(sform)
    icorr.set_qform(qform)
    icorr.header['sform_code']=scode
    icorr.header['qform_code']=qcode
    fix_corr= os.path.abspath('fixed_'+os.path.split(orig)[1])
    icorr.to_filename(fix_corr)

    ifield.set_sform(sform)
    ifield.set_qform(qform)
    ifield.header['sform_code']=scode
    ifield.header['qform_code']=qcode   
    fix_field = os.path.abspath('fixed_'+os.path.split(field)[1])
    ifield.to_filename(fix_field)
    
    fix_fieldc= os.path.abspath('fixed_'+os.path.split(fieldcoef)[1])  
    ifieldc.to_filename(fix_fieldc)
    
    # copy the motion parameters
    # Downstream, Eddy takes in a whole topup folder, 
    # so everything has to be in it.
    fix_movpar=os.path.abspath('fixed_'+os.path.split(movpar)[1])
    status=subprocess.call(['cp',movpar,fix_movpar]) 
    if status==1:
        raise IOError('Unable to write '+fix_movpar) 
    return(fix_corr, fix_field, fix_fieldc, fix_movpar)


def prepare_eddy(b0names, dwinifti,bvals,bvecs, thr):
    '''
    This functions in meant to produce concatenated bvecs and bvals,
    and the index file for the FSL Eddy script with a topup input.
    The inputs are essentially the ouput of Rorden's dcm2nii:
    Takes in a list of nifti 4D DWI series (with b0+actual DWI images), 
    the corresponding list of bval/bvecs (!! IN THE SAME ORDER !!), 
    and the b0 threshold (something like 6).
    The NIfTI are merely there for getting the number of frames.
    '''
    
    import numpy as np
    import nibabel as nb    
    import os
    import re
    
    index=[]
    bigbval=[]
    for (c,nifti_file) in enumerate(dwinifti):
        nifti_vol=nb.load(nifti_file)
        bval = np.genfromtxt(bvals[c],delimiter=" ")
        bvalues=bval.tolist()
        dims= nifti_vol.header.get_data_shape()
        series=int(re.findall('s(.*)a',nifti_file.split('.')[0][-9:])[0])
        idx = [i for (i, val) in enumerate(b0names) if int(val) == series] 
        dwidx = [i for (i, val) in enumerate(bvalues) if val < thr]
        bounds=[i-1 for i in dwidx]
        bounds.reverse()
        bounds.pop()
        bounds.reverse()
        steps=[dims[3]-1]
        bounds.append(steps[0])
        if len(bounds)>1:
            steps=[]
            for i in range(len(dwidx)):
                steps.append(bounds[i]-dwidx[i])
        for i in range(len(steps)):
            index = index + [idx[0]+i] * (steps[i]+1)
        bigbval = bigbval + bvalues
        bvec=  np.genfromtxt(bvecs[c],delimiter=" ")
        if c==0:
            bigbvec = bvec
        else: 
            bigbvec = np.concatenate((bigbvec,bvec),1)
    index=[i+1 for i in index] 
    fid=open('index','w')
    for i in index:
        fid.write("%s " % i)
    fid.write("\n")
    fid.close()
    b0_indices = [float(i) for i,bv in enumerate(bigbval) if bv < thr]
    index_out=os.path.abspath('index')  
    np.savetxt('bvec', bigbvec, fmt='%.8f')    
    bvec_out=os.path.abspath('bvec')
    np.savetxt('bval', bigbval, fmt='%.8f')    
    bval_out=os.path.abspath('bval')
    b0_topup_indices=list(set(index))
    return (index_out,bval_out,bvec_out,b0_topup_indices,b0_indices)


# for creating wide-format (as opposed to tall-format) numeric tables   
# No column names, no row names either, matrices have to be the same sizes.

def mergeTables(in_files):
    '''
    Merges two matrices by rows. 
    No column names, no row names either, matrices have to be the same sizes.
    '''
    import numpy as np
    import os as os
    
    vectors = []
    for idx, filename in enumerate(in_files):
        vectors.append(np.genfromtxt(filename))
        if (vectors[idx].shape).__len__() == 1:
            vectors[idx]=vectors[idx].reshape(vectors[idx].shape[0],1)
    cvectors = np.hstack(vectors)
    out_file = os.path.join(os.getcwd(), "merged_array.txt")
    np.savetxt(out_file, cvectors, fmt=str("%.10f"))
    return out_file

#%%  Actual scientific computations   
def extract_atlas_ROIs(timeseries_file,mask_file, atlas_file, atlas_xml,eig=True):
    """
    Extract average voxel time-courses for each ROI in an FSL/xml-style atlas

    Parameters
    ----------

    timeseries_file: a 4D Nifti file
    atlas_file: a 3D file containing rois in the same space/size of the 4D file
    atlas_xml: the list of regions in fslview xml format
    eig: compute SPM eigenvariates instead of mean timecourse over voxels

    Returns
    -------
    out_file: a text file containing time courses for each ROI
    
    from ginnipi.toolbox.computations import extract_atlas_ROIs
    extract_atlas_ROIs(atlas_file='/usr/local/fsl/data/atlases/AICHA/AICHA-2mm.nii.gz',atlas_xml='/usr/local/fsl/data/atlases/AICHA_Short.xml',timeseries_file='/netapp/vol6_gin/herve/iShareTestDrive/XNATABC2/ABC/_subject_id_MRi_ShareDB_E00681/func2/FuncPart1/fsl_merge/wvol0000_merged.nii.gz',mask_file='/netapp/vol6_gin/herve/iShareTestDrive/XNATABC2/ABC/_subject_id_MRi_ShareDB_E00681/func2/FuncPart1/spmWarpBrainMasktoStandard/wbrainmask_out_maths.nii',eig=True)
    """
    import nibabel as nb
    import os 
    import xml.etree.ElementTree as ET
    import numpy as np
    from scipy.stats import pearsonr
    
    img = nb.load(timeseries_file)
    data = img.get_data()
    roiimg = nb.load(atlas_file)
    rois   = roiimg.get_data()
    mskimg = nb.load(mask_file)
    msk = mskimg.get_data()
    atlas  = ET.parse(atlas_xml)
    labels = atlas.findall(".//label")
    name   = atlas.findall('.//shortname')[0].text.split()[0]
    if not eig:
        out_file = os.path.join(os.getcwd(), name+'_timeseries.csv')
    else:
        out_file = os.path.join(os.getcwd(), name+'_eigenvariate_timeseries.csv')
    print out_file
    with open(out_file, 'wt') as fp:
        for label in labels:
            print label.text
            ijk = np.nonzero((rois == int(label.get('index'))) & (msk>0))
            chunk =  np.transpose(data[ijk])
            chunk=chunk[np.isfinite(chunk).any(axis=1)]
            if not eig:
                ts = np.mean(chunk,axis=1)
            else:
                failed=False
                # SPM matlab to numpy code translation
                m,n   = np.shape(chunk)
                if (m > n):      # more observation than variables
                    print "first case"
                    print n
                    try:
                        V, s, Vp = np.linalg.svd(np.dot(np.transpose(chunk),chunk),full_matrices=True)
                        V = V[:,0]
                        U = np.dot(chunk,V/np.sqrt(s[0]))
                    except:
                        print('SVD failed to converge, taking mean instead')
                        print np.shape(V)
                        print np.shape(chunk)
                        ts = np.mean(chunk,axis=1)
                        failed=True
                else:      # more variables than observation
                    print "second case"
                    print n
                    try:
                        U, s, Up = np.linalg.svd(np.dot(chunk,np.transpose(chunk)),full_matrices=True)
                        U = U[:,0]
                        V = np.dot(np.transpose(chunk),U/np.sqrt(s[0]))
                    except:
                        print('SVD failed to converge, taking mean instead')
                        ts = np.mean(chunk,axis=1)
                        failed=True
                if not failed:
                    d = np.sign(np.sum(V))
                    U = np.dot(U,d);
                    ts = U*np.sqrt(s[0]/n) 
                    ts2 = np.mean(chunk,axis=1)
                    print pearsonr(ts,ts2)
                    print np.median(ts)
                    print np.median(ts2)
                print"-----------"
            fp.write('%s,' % (label.text) + ','.join(['%.10f' % val for val in ts]) + '\n')
    return out_file

def compute_correlation_matrix(in_file):
    """
    Compute a Pearson correlation matrix from a table of ROI time series 

    Parameters
    ----------
    in_file: a .csv with time series in rows and names of ROIs in first column.    
    
    Returns
    -------
    out_file : a correlation matrix in text format (semi-comma delimited, with row and column names)
    out_file2: a correlation matrix in text format (comma separated, without row and column names)
    """
    import numpy as np
    import os
    import pandas
    import csv
    ts = pandas.read_csv(in_file,sep=',',header=None,  index_col=0)
    cormat=np.corrcoef(ts.values[:,:])
    cmap = pandas.DataFrame(cormat,columns=list(ts.index.values),index=list(ts.index.values))  # Pearson...
    cmap2 = pandas.DataFrame(cormat)
    out_file = os.path.join(os.getcwd(), os.path.basename(in_file).split('.')[0]+'_cormat.txt')
    out_file2 = os.path.join(os.getcwd(), os.path.basename(in_file).split('.')[0]+'_cormat4jgex.txt')
    cmap.to_csv(out_file,sep=';',quotechar='"',decimal=".",quoting=csv.QUOTE_NONNUMERIC)
    cmap2.to_csv(out_file2,sep=',',quotechar='"',decimal=".",quoting=csv.QUOTE_NONNUMERIC,index=False,header=False)
    return (out_file,out_file2)
    
def compute_VMHC(in_file):
    """
    Computes homotopic pearson correlation coefficients (voxelwise method)
    This assumes the inter-hemispheric plane is in the sagittal plane, at 
    a particular voxel or exactly between two voxels.
    This function can deal with the plane being shifted from the middle of 
    the volume in voxel units, but not with subvoxel shifts.
     
    Parameters
    ----------
    in_file: a preprocessed time-series NIfTI volume
    
    Returns
    -------
    out_file: a 3D hemispheric map of homotopic pearson correlation coefficients  
    """  
    import nibabel as nb
    import numpy as np
    import os
    
    ima = nb.load(in_file)
    ima_data=ima.get_data()
    origin=np.linalg.inv(ima.affine) * np.transpose(np.matrix([0.,0.,0.,1.0]))
    x0=float(origin[0])
    xA=ima.affine * np.transpose(np.matrix([0.,0.,0.,1.0]))
    xB=ima.affine * np.transpose(np.matrix([float(ima.shape[0]-1),0.,0.,1.0]))
    xA=float(xA[0])
    xB=float(xB[0])
    if (xA+xB) != 0:
        print("Warning: the center of the volume does not coincide with the inter-hemispheric plane!" )
        print("Check results carefully.")
        if abs(xA)>abs(xB):
            A = np.linalg.inv(ima.affine) * np.transpose(np.matrix([-xB,0.,0.,0.,1.0]))
            A = int(A[0]) 
        else:
            B = np.linalg.inv(ima.affine) * np.transpose(np.matrix([-xA,0.,0.,1.0]))
            B = int(B[0])       
    else:
        A=0
        B=ima.shape[0]-1
    xlim=int(x0)
    if int(np.mod(x0,ima.affine[0,0])!=0.0): # odd case   ( 0 is middle of voxel, ex: -2;0;2, skip midline) 
        print('Skipping midline voxels (odd number of voxels)')
        offset=1
    else:               # even case (0 is between 2 voxels, ex: -3;-1;+1;+3, no midline)
        print('Center of volume between two voxels (even number of voxels)')        
        offset=0
    Amat = ima_data[range(A,xlim+1-offset),:,:,:]
    Ashape=Amat.shape
    Bmat = ima_data[range(xlim+offset,B+1),:,:,:]
    Bmat=Bmat[::-1,:,:,:]
    # Reshape and reduce
    Amat=Amat.reshape(np.prod(np.array(Amat.shape[0:3])),Amat.shape[3])
    Bmat=Bmat.reshape(np.prod(np.array(Bmat.shape[0:3])),Bmat.shape[3])
    idxes = np.where(np.any(Amat,axis=1) & np.any(Bmat,axis=1))[0].tolist()
    Amat=Amat[idxes,:]
    Bmat=Bmat[idxes,:]
    # compute VMHC
    Am=np.mean(Amat,1)
    Bm=np.mean(Bmat,1)
    Amat=(Amat.transpose()-Am).transpose()
    Bmat=(Bmat.transpose()-Bm).transpose()
    r=np.divide(np.sum(np.multiply(Amat,Bmat),1),np.multiply(np.sqrt(np.sum(np.multiply(Amat,Amat),1)),np.sqrt(np.sum(np.multiply(Bmat,Bmat),1))))
    #reshaping back
    Cmap=np.zeros([np.prod(np.array(Ashape[0:3])),1],dtype=float)
    Cmap[idxes,0]=r
    Cmap=Cmap.reshape(Ashape[0:3])
    Cmap2=np.ndarray(ima.shape[0:3])*0
    Cmap2[range(A,xlim+1-offset),:,:]=Cmap
    out = nb.Nifti1Image(Cmap2, ima.affine)
    outf = os.path.join(os.path.abspath('.'),os.path.split(in_file)[1].split('.')[0] + '_VMHC.nii.gz')
    out.to_filename(outf)
    return outf
           
def create_gexf_graph(in_file,threshold=None):
    """Creates a weighted or non-weighted GEXF graph from a correlation matrix 

    Parameters
    ----------
    in_file: a correlation matrix (csv file with row names and column names)  
    threshold: threshold applied to correlation matrix
    
    Returns
    -------
    out_file: a .gexf graph
    
    """
    import networkx as nx
    import os
    import pandas
    
    corr_mat = pandas.read_csv(in_file,sep=';',quotechar='"',decimal=".",index_col=0,header=0)
    if not threshold:
        # create the graph
        # we don't threshiold the values in the array.
        # Since we don't know if these are distances 
        # or proximities on the edges, it is better to
        #  keep metrics computations for latter...
        mygraph = nx.from_numpy_matrix(corr_mat.values)
    else: 
        # we threshold the graph and also compute
        # simple graph metrics, because there is no ambiguity
        mygraph = nx.from_numpy_matrix(corr_mat.values>threshold)
        bwc = nx.betweenness_centrality(mygraph)
        dgc = nx.degree_centrality(mygraph)
        clc = nx.closeness_centrality(mygraph)
        evc = nx.eigenvector_centrality(mygraph)
        labels = list(corr_mat.index.values)
        for idx,label in enumerate(labels):
            mygraph.node[idx]['label'] = label
            mygraph.node[idx]['betweenness_centrality'] = bwc[idx]
            mygraph.node[idx]['closeness_centrality'] = clc[idx]
            mygraph.node[idx]['degree_centrality'] = dgc[idx]
            mygraph.node[idx]['eigenvector_centrality'] = evc[idx]
    out_file=os.path.join(os.getcwd(), os.path.basename(in_file).split('.')[0]+'_graph.gexf')
    nx.write_gexf(mygraph,out_file)
    return out_file
        
def create_boldqc_assessor(motion,realigned,mask,rms_files):
    '''
    Creates an assessor suitable for the extended BOLD QC xnat module,
    along with a few plots.
    UNFINISHED
    '''
    import pandas
    import nibabel as nb
    import numpy as np
    import lxml.etree as ET
    import matplotlib.pyplot as plt 
    
    mpar = pandas.read_csv(motion,decimal=".",delim_whitespace=True,header=None)
    msk=nb.load(mask)
    imask=msk.get_data()
    idxes=imask==1
    ts= nb.load(realigned)
    its= ts.get_data()
    number_of_volumes = np.shape(ts)[3]
    number_of_slices = np.shape(ts)[2]
    number_of_motion_timepoints = np.shape(mpar)[0]
    grand_mean_within_mask = np.mean(its[idxes,:])    
    grand_sd_within_mask = np.std(its[idxes,:]) 
    mean_ts= [np.mean(its[idxes,i]) for i in range(number_of_volumes)]
    mean_ts_per_slice = np.zeros(shape=(number_of_slices,number_of_volumes), dtype=float)
    sd_per_slice = np.zeros(shape=(number_of_slices,0), dtype=float)
    for s in range(number_of_slices):     
        ts_slice = np.squeeze(its[:,:,s,:])
        idx_slice = np.squeeze(idxes[:,:,s]) 
        sd_per_slice[s] = np.std(ts_slice[idx_slice,:]) 
        if (np.sum(idx_slice) > 0):
            for n in range(number_of_volumes):  
                mean_ts_per_slice[s,n] = np.mean(ts_slice[idx_slice,n]) 
        else:
            mean_ts_per_slice[s,:]=0.0
    mean_sd_per_slice=np.mean(sd_per_slice[sd_per_slice>0])
    fig=plt.figure(1, figsize=(18,6))
    slicesplot=plt.plot(np.transpose(mean_ts_per_slice[np.sum(mean_ts_per_slice,1)>0,:]))
    plt.xlabel('Volume')
    plt.ylabel('Slice average intensity')
    fig.savefig('slice_average_tc.png')
    fig.clear()
    del fig
    ampar=mpar.as_matrix()
    rmpar = ampar[1:(number_of_motion_timepoints-1),:]-ampar[:(number_of_motion_timepoints-2),:]
    rmpar = np.concatenate((np.zeros(shape=(1,6)),rmpar))
    rel_means=np.mean(abs(rmpar),axis=0)
    rel_sd=np.std(abs(rmpar),axis=0)
    abs_means=np.mean(abs(ampar),axis=0)
    abs_sd=np.mean(abs(ampar),axis=0)
    fig2=plt.figure(2, figsize=(18,6))
    parplot=plt.plot(rmpar[:,:3])
    fig2.savefig('motion_parameters_pythonplot.png')
    fig2.clear()
    del fig2
    eqc=ET.Element('root') 
    return(eqc)    
    
def volume_sbc_maps(func, mask, seeds, rownames=[],ctype='pearson'):
    '''
    Correlation map from a seed time course file.
    seed timeseries should be a text file with values separated by spaces 
    and as many values as timepoints in the nifti file.
    
    Inputs
    ------
    func: filtered functional data
    mask: mask file in same geometry as func
    seeds: a text file containing the time-courses to regress, with time in colum and region in rows
           (rownames are necessary)
    rownames: a list of names that are actual row names in the seeds file (optional)
    ctype: pearson (default), spearman or kendall
    
    Outputs
    -------
    Writes a series of correlation maps and returns the file names
    '''
    import nibabel as nb
    import pandas  as pnd
    import numpy as np 
    import os
    from scipy.stats import pearsonr, spearmanr, kendalltau
    
    # load the complete seed data
    print 'loading seeds'
    try:
        # timeseries are in rows, with first columns for region names, comma separated, no header 
        seedc = pnd.read_csv(seeds,sep=',',header=None,  index_col=0)
    except:
        print('failed to load seed time course '+seeds)
        raise
    
    # select the seeds to use, otherwise process all seeds...
    # and this may generate a lot of data of course.
    print 'selecting seeds'
    if rownames:
        try:
            seedc = seedc.loc[rownames]
        except:
            print "one or more row names could not be found in supplied seed file"
            raise
    else:
        rownames=list(seedc.index)
        
    # load mask data
    print 'loading mask'
    try:
        ima=nb.load(mask)
        ijk= np.nonzero(ima.get_data()>0)
    except:
        print('failed to load analysis mask '+mask)
        raise 
    
    # load filtered func data
    print 'loading BOLD time series...this may take a while'
    try:
        fima=nb.load(func)
        its= fima.get_data()
    except:
        print('failed to load filtered functional data '+func+' within mask '+mask)
        raise
    
    # Compute pearson correlation maps
    out_files=[]
    for row in rownames:
        print 'computing seed based correlation map for seed ' + row + '...'
        rmap=np.zeros(shape=ima.header.get_data_shape())
        seed=np.array(seedc.loc[row])
        # compute correlation at each voxel
        if ctype=='pearson':
            rmap[ijk] = [pearsonr(seed,its[ijk[0][idx],ijk[1][idx],ijk[2][idx],:])[0] for idx in range(len(ijk[0]))]  # Banzaï !!
        elif ctype=='spearman':
            rmap[ijk] = [spearmanr(seed,its[ijk[0][idx],ijk[1][idx],ijk[2][idx],:])[0] for idx in range(len(ijk[0]))]  
        elif ctype=='kendall':
            rmap[ijk] = [kendalltau(seed,its[ijk[0][idx],ijk[1][idx],ijk[2][idx],:])[0] for idx in range(len(ijk[0]))]  
        else:
            raise ValueError('Wrong correlation type...either one of pearson, spearman or kendall')
        fname=os.path.abspath(row+'_sbcmap.nii.gz')
        print 'writing out: ' + fname
        nifti_out = nb.Nifti1Image(rmap, ima.affine)
        nifti_out.to_filename(fname)
        out_files.append(fname)
    print('done.')
    return out_files

    #  sbc_maps(func='/netapp/vol6_gin/herve/iShareTestDrive/XNATABC2/ABC/_subject_id_MRi_ShareDB_E00681/func2/FuncPart1/afniStereotaxicBlur/wvol0000_merged_blur.nii.gz',
    #           mask='/netapp/vol6_gin/herve/iShareTestDrive/XNATABC2/ABC/_subject_id_MRi_ShareDB_E00681/func2/FuncPart1/spmWarpBrainMasktoStandard/wbrainmask_out_maths.nii',
    #           seeds='/homes_unix/herve/AICHA_eigenvariate_timeseries.csv',
    #           rownames=['Frontal_Sup-3-L','Frontal_Sup-3-R'])


def freesurfer_sbc_maps(func, seeds, rownames=[],ctype='pearson'):
    '''
    Surface-based correlation map from a seed time course file.
    seed timeseries should be a text file with values separated by spaces 
    and as many values as timepoints in the .mgz file.
    
    Inputs
    ------
    func: filtered functional data in .mgz format
    seeds: a text file containing the time-courses to regress, with time in colum and region in rows
           (rownames are necessary)
    rownames: a list of names that are actual row names in the seeds file (optional)
    ctype: pearson (default), spearman or kendall
    
    Outputs
    -------
    Writes a series of correlation maps in .mgz format and returns the file names
    '''
    import nibabel as nb
    import pandas  as pnd
    import numpy as np 
    import os
    from scipy.stats import pearsonr, spearmanr, kendalltau
    
    # load the complete seed data
    print 'loading seeds'
    try:
        # timeseries are in rows, with first columns for region names, comma separated, no header 
        seedc = pnd.read_csv(seeds,sep=',',header=None,  index_col=0)
    except:
        print('failed to load seed time course '+seeds)
        raise
    
    # select the seeds to use, otherwise process all seeds...
    # and this may generate a lot of data of course.
    print 'selecting seeds'
    if rownames:
        try:
            seedc = seedc.loc[rownames]
        except:
            print "one or more row names could not be found in supplied seed file"
            raise
    else:
        rownames=list(seedc.index)
    
    # load surface-sampled filtered func data
    print 'loading BOLD time series...this may take a while'
    try:
        fima=nb.load(func)
        its= fima.get_data()
    except:
        print('failed to load filtered functional data '+func)
        raise
    
    # Compute pearson correlation maps
    out_files=[]
    nvertices=fima.header.get_data_shape()[0]
    for row in rownames:
        print 'computing seed based correlation map for seed ' + row + '...'
        rmap=np.zeros(shape=fima.header.get_data_shape()[:3],dtype=np.float32)
        seed=np.array(seedc.loc[row])
        # compute correlation at each voxel
        if ctype=='pearson':
            rmap[:,0,0] = [np.float32(pearsonr(seed,its[idx,0,0,:])[0]) for idx in range(nvertices)]  # Banzaï !!
        elif ctype=='spearman':
            rmap[:,0,0] = [np.float32(spearmanr(seed,its[idx,0,0,:])[0]) for idx in range(nvertices)]  
        elif ctype=='kendall':
            rmap[:,0,0] = [np.float32(kendalltau(seed,its[idx,0,0,:])[0]) for idx in range(nvertices)]  
        else:
            raise ValueError('Wrong correlation type...either one of pearson, spearman or kendall')
        fname=os.path.abspath(row+'_'+os.path.splitext(os.path.split(func)[1])[0]+'_sbcmap.mgz') # with the freesurfer lh and rh split...
        print 'writing out: ' + fname
        mgh_out = nb.MGHImage(rmap, fima.affine)
        mgh_out.to_filename(fname)
        out_files.append(fname)
    print('done.')
    return out_files

# def test_sbc():
#     '''
#     test function
#     '''
#     from ginnipi.toolbox.computations import freesurfer_sbc_maps
# 
#     fsfile='/data/vol6_gin/herve/iShareTestDrive/XNATABCTEST/ABC/_subject_id_ABACI_E00222/fsfunc/fsFunc/lhBOLDSurfaceSmooth/lh.vol0000_warp_merged_bp_smooth5.mgz'
#     aicha_ts='/data/vol6_gin/herve/iShareTestDrive/XNATABCTEST/ABC/_subject_id_ABACI_E00222/func2/FuncPart1/AICHA_ts_extract/AICHA_timeseries.csv'
#     regions=['Frontal_Sup-3-L','Frontal_Sup-3-R','Rolando-3-R','Rolando-3-L','Angular-2-L','Angular-2-R']
#     f=freesurfer_sbc_maps(func=fsfile,seeds=aicha_ts,rownames=regions)
#     return f

def generateAtlasAssessorXsd():
    '''
    From an FSK XML atlas definition, or JSON atlas definition gets an xsd for xnat
    '''
    pass

def generateFslAtlasAssessorXml():
    #import nibabel as nb
    #import numpy as np
    #import lxml.etree as et
    pass

def createVbmAssessorXml(gm,
                         wm,
                         csf,
                         aal,
                         template='TPM.nii',
                         algo='SPM12',
                         atlas='AAL for SPM12'):
    '''
    Computes total and tissue regional tissue volumes
    from jacobian modulated maps and packages this in an XNAT xml assessor
    '''
    import nibabel as nb
    import numpy as np
    import lxml.etree as et
    import os
    
    def get_integrals(imafile,mask=None,roi_idx=None):
        '''
        Computes a volume from a tissue classification map
        '''
        
        ima=nb.load(imafile)
        if len(np.shape(ima))!=3:
            raise Exception('image is not a 3D volume')
        vox_vol=np.prod(ima.header['pixdim'][1:4])
        if not mask:
            ima_volume = np.nansum(ima.get_data())*vox_vol          
            return ima_volume
        else:
            if not roi_idx:
                raise Exception('Please provide a ROI index')
            msk=nb.load(mask)
            mski = (msk.get_data()==np.int(roi_idx)).astype(np.float)
            try:
                ima_volume = np.nansum(np.multiply(ima.get_data(),mski))*vox_vol          
                return ima_volume
            except:
                return np.NaN
            
        
    # Compute the volumes
    # assuming mm units
    vol_grey = get_integrals(gm)
    print('gm: '+str(vol_grey/1000.0)+' CC')  
    
    vol_white = get_integrals(wm)
    print('wm: '+str(vol_white/1000.0)+ ' CC')
    
    vol_csf = get_integrals(csf)
    print('csf: '+str(vol_csf/1000.0)+ ' CC')
    
    vol_hipL = 0
    vol_hipR = 0
    try:
        vol_hipL = get_integrals(gm, mask=aal,roi_idx=4101)
        print('Hip L: '+str(vol_hipL/1000.0)+ ' CC')
        vol_hipR = get_integrals(gm, mask=aal,roi_idx=4102)
        print('Hip R: '+str(vol_hipR/1000.0)+ ' CC')
    except:
        print('Failed to compute hippocampal volumes. Check atlas size.')    
    # Now package that into a VBM assessor
    vbmass = et.ElementTree(et.XML('''<vbm:vbmAssessorData xmlns:vbm="http://nrg.wustl.edu/vbm" xmlns:xnat="http://nrg.wustl.edu/xnat" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"></vbm:vbmAssessorData>'''))
    
    namespaces={'vbm':'http://nrg.wustl.edu/vbm',
                'xnat':'http://nrg.wustl.edu/xnat',
                'xsi':'http://www.w3.org/2001/XMLSchema-instance'}
    
    GM_vol = et.Element(et.QName(namespaces['vbm'],"GM_Volume"),nsmap=namespaces)
    GM_vol.text = str(vol_grey) 
    vbmass.getroot().append(GM_vol)
    
    WM_vol = et.Element(et.QName(namespaces['vbm'],"WM_Volume"),nsmap=namespaces)
    WM_vol.text = str(vol_white)  
    vbmass.getroot().append(WM_vol)
    
    CSF_vol = et.Element(et.QName(namespaces['vbm'],"CSF_Volume"),nsmap=namespaces)
    CSF_vol.text = str(vol_csf)
    vbmass.getroot().append(CSF_vol)
    
    Stereo = et.Element(et.QName(namespaces['vbm'],"StereotaxicSpace"),nsmap=namespaces)      
    Stereo.text = template
    vbmass.getroot().append(Stereo)
    
    Algo = et.Element(et.QName(namespaces['vbm'],"AlgoNorm"),nsmap=namespaces)      
    Algo.text = algo
    vbmass.getroot().append(Algo)
    
    HIPL = et.Element(et.QName(namespaces['vbm'],"HYPL_Volume"), nsmap=namespaces)
    HIPL.text=str(vol_hipL)
    vbmass.getroot().append(HIPL)
    
    HIPR = et.Element(et.QName(namespaces['vbm'],"HYPR_Volume"),nsmap=namespaces)
    HIPR.text=str(vol_hipR)
    vbmass.getroot().append(HIPR)
   
    AtlasHip = et.Element(et.QName(namespaces['vbm'],"AtlasHyp"),nsmap=namespaces)      
    AtlasHip.text = atlas
    vbmass.getroot().append(AtlasHip)
    
    baseScan = et.Element(et.QName(namespaces['vbm'],"baseScanNumber"),nsmap=namespaces)
    baseScan.text='0'
    vbmass.getroot().append(baseScan)
    
    # done generating XML, write it...
    fname=os.path.abspath('vbmAssessor.xml')
    with open(fname,'w') as f:
        vbmass.write(f,pretty_print=True,encoding='utf-8', xml_declaration=True)
    f.close()
    return fname
    
def createFslAtlasAssessorXml(gm,wm,csf,aal):
    '''
    From an FSL XML atlas, creates an .xsd structure for extending the xnat data model.
    '''
    import nibabel as nb
    import numpy as np
    import lxml.etree as et
    pass 


def strip_multi_shell(dwi, 
                      bvals, 
                      bvecs,
                      out_fileroot='stripped_',
                      cutoff=1510.0, 
                      upper=True, 
                      b0_cutoff=-1.0):
    '''
    Sometimes it is necessary to remove a shell from an aggregated multi-shell DWI 
    4D NIfTI and associated files into their base components. 
    
    The default settings are for removing the frames with b-values above 1510,
    while keeping the b0 scans ('b0_cutoff' is negative). If you want to remove b0 scans too,
    set 'b0_cutoff' to above their b-value (about 6). 
    
    Conversely, in order to extract the b-values above the cutoff, set 'upper' to False .
    In this case, set 'b0_cutoff' appropriately (about 6) in order to keep the b0 scans, 
    or to a negative values if you want to remove b0s too.
    '''
    import nibabel as nb
    import numpy as np
    import os
    
    # read bvals
    try:
        bvalues=np.genfromtxt(bvals,delimiter=" ")
    except:
        print('Failed to load '+bvals)
        raise
    
    # read bvecs
    try:
        bvectors=np.genfromtxt(bvecs,delimiter=" ")
    except:
        print('Failed to load '+bvecs)
        raise
    
    # read image
    ima=nb.load(dwi)
    
    # perform a sanity check    
    dims=np.shape(bvectors)
    if dims[1] != np.shape(bvalues)[0]:
        raise Exception('bvalues and bvecs don''t have the same length!')
    if np.shape(ima)[3] != np.shape(bvalues)[0]:
        raise Exception('bvalues and bvecs don''t have the same length as the NIfTI file!')
    
    # select the frames to keep
    if upper:
        keep_idx=np.squeeze(np.nonzero(np.logical_and(bvalues<cutoff, bvalues > b0_cutoff)))
    else:        
        keep_idx=np.squeeze(np.nonzero(np.logical_or(bvalues>cutoff, bvalues < b0_cutoff)))
    
    # strip the 4D NIfTI
    out_ima=nb.Nifti1Image(ima.get_data()[:,:,:,keep_idx],affine=ima.affine)
    
    # save it
    out_dwi=os.path.abspath(out_fileroot+os.path.split(dwi)[1])
    out_ima.to_filename(out_dwi)
    
    # strip bval and bvec 
    out_bvals=os.path.abspath(out_fileroot+os.path.split(bvals)[1])
    out_bvecs=os.path.abspath(out_fileroot+os.path.split(bvecs)[1])
    
    # save  bval and bvec 
    np.savetxt(out_bvecs, bvectors[:,keep_idx], fmt='%.8f')  
    np.savetxt(out_bvals, bvalues[keep_idx], fmt='%.8f') 
        
    return (out_bvals, out_bvecs, out_dwi)

def combine_labels(atlas, combinations, 
                   out_filename_root='recombination',
                   obs_name=None, 
                   out_name=None, 
                   multi_label=True):
    '''
    Assembles a set of regions into a new label image
    
    Numbers should be a dictionary with keys being the anatomical label, 
    and with values a list of the regions that should be put together:
    
    [ ('left_frontal',[1,2,3,4,5,6]),('right_frontal',[7,8,9,10,11,12]) ]
     
    Will produce either a set of 3D images with single labels 
    or a single 3D image with an FSL-like XML look-up table
    if multi_label is set to True. 
    
    returns a list of path to NIfTI files, named after the 'keys' in the combinations
        or 
    returns a single NIfTI file path + path to xml FSL-style LUT (multi_label)
    '''
    import nibabel as nb
    import numpy as np
    import lxml.etree as et
    import os
    import re
    import pandas as pd
    # Load label image, convert it to int16
    ima = nb.load(atlas) 
    in_atlas = np.int16(ima.get_data())
    if not out_name:   
        out_name = 'custom recombination of '+ os.path.split(atlas)[1]
    
    # we need to handle multiple occurrences of combination labels
    # in the supplied list
    labels=dict()
    
    #%% multi-label case
    if multi_label:
        print('Producing a multi-label image')
        # create the stub for the atlas xml document

        fsl_xml = et.fromstringlist(['<?xml version="1.0" encoding="UTF-8"?>' , 
                                '<atlas version="1.0">', 
                                '<header>', 
                                '<name>'+out_name+'</name>', 
                                '<type>Label</type>', 
                                '<images>', 
                                '<imagefile>' + out_filename_root  + '</imagefile>', 
                                '<summaryimagefile>' + out_filename_root + '</summaryimagefile>', 
                                '</images>',  
                                '</header>', 
                                '<data>', 
                                '</data>', 
                                '</atlas>']
                                )
        data_elt = fsl_xml.find('./data')
        # preallocate 3D volume of the right shape
        rc_atlas = np.zeros(shape=ima.header.get_data_shape()[:3],dtype=np.int16)  
        # parse the combination tuples one after the other
        label_elt = et.Element('label', attrib={'index':'0'})
        label_elt.text = 'outside'
        data_elt.append(label_elt)

        # label per label
        if obs_name:
            volumes=[(re.sub('[:;,$\s\-\(\)\[\]]','_',out_name)+'/ROI/volume',obs_name)]
        else:
            volumes=[(re.sub('[:;,$\s\-\(\)\[\]]','_',out_name)+'/ROI/volume',ima)]
        for idx,combination in enumerate(combinations):
            # check if label wasalready seen:
            if not combination[0] in labels.keys():
                labels[combination[0]] = 1
            else:
                labels[combination[0]] += 1
                # if it this is the case append the numerical index
                combination[0]=combination[0] + '_'  + str(labels[combination[0]])
            
            # perform the combination
            volume = 0
            for intensity in combination[1]:
                ijk = np.nonzero(in_atlas == np.int16(intensity))
                if ijk:
                    rc_atlas[ijk] = np.int16(idx + 1)
                    volume= volume + np.shape(ijk)[1]*np.prod(ima.header['pixdim'][1:4])
            
            # add the label to the atlas xml LUT
            label_elt = et.Element('label', attrib={'index':str(idx + 1)})
            label_elt.text = combination[0]
            data_elt.append(label_elt)
            volumes.append( (label_elt.text,[volume]))
            
        # save multi-label 3D image + companion XML file
        out_fname =os.path.abspath(out_filename_root+'.nii.gz')
        out_xml = os.path.abspath(out_filename_root+'.xml')
        out_csv = os.path.abspath(out_filename_root+'_volumes.csv')
        print('writing out: ' +  out_fname)
        nifti_out = nb.Nifti1Image(rc_atlas, ima.get_affine())
        nifti_out.set_qform(ima.header.get_qform())
        nifti_out.to_filename(out_fname)
        print('writing out: ' +  out_xml)
        fsl_xml.getroottree().write(out_xml,pretty_print=True)
        
        # also add volumes a .csv line with volumes
        df_out=pd.DataFrame.from_items(volumes)
        df_out.to_csv(out_csv, index=False)
        
        return (out_fname, out_xml, out_csv)   
    else:
        #%% single-label case
        print('Producing a series of single-label images')
        out_files=list()
        #volumes=list()
        for idx,combination in enumerate(combinations):
            rc_atlas = np.zeros(shape=ima.header.get_data_shape()[:3],dtype=np.uint8)
            
            # check if label was already seen:
            if not combination[0] in labels.keys():
                labels[combination[0]] = 1
            else:
                labels[combination[0]] += 1
                # if it this is the case append the numerical index
                combination[0]=combination[0] + '_'  + str(labels[combination[0]])
            
            # perform the combination
            #volume = 0
            for intensity in combination[1]:
                ijk = np.nonzero(in_atlas == np.int16(intensity))
                if ijk:
                    rc_atlas[ijk] = np.uint8(1)
            #        volume= volume + np.shape(ijk)[2]*np.prod(ima.header['pixdim'][1:4])
            #volumes.append(volume)
            # save 3D image
            fname = os.path.abspath(out_filename_root + '_' + combination[0] + '.nii.gz').replace(' ','_')
            out_files.append(fname)
            print('writing out: ' + fname)
            nifti_out = nb.Nifti1Image(rc_atlas, ima.get_affine())
            nifti_out.set_qform(ima.header.get_qform())
            nifti_out.to_filename(fname)    
        return out_files

# the following two are helpful for ROI sampling
# can make the process independent on bounding-box.
def ijk2xyz(sform,ijk):
    '''
    Voxel numpy array indices (tuple with I, J and K index arrays) 
    to real-world coordinates
    '''
    import numpy as np
    ijk=np.array([ijk[i]  for i in range(3)]).astype(float).transpose()
    ijk=np.concatenate((ijk,np.ones((np.shape(ijk)[0],1))),1)
    xyz = np.dot(sform,ijk.transpose()).transpose()
    return xyz[:,range(3)]
    

def xyz2ijk(sform,xyz):
    '''
    Real-world coordinates array (X, Y, Z  as columns)
    to voxel numpy array indices
    '''
    import numpy as np
    xyz = np.concatenate((xyz,np.ones((np.shape(xyz)[0],1))),1).transpose()
    ijk = np.round(np.dot(np.linalg.inv(sform),xyz)).astype(int)
    return tuple([ ijk[i,:] for i in range(3)])


def stats_from_fslxml_atlas(ima, 
                            atlas_ima, 
                            atlas_xml,
                            meas,
                            stat='mean', 
                            obs_name=None, 
                            mask=None,
                            mask_thr=0,
                            out_fileroot='stats_table',
                            jacobian = None):
    '''
    Given a matching quantitative image and atlas image + fsl xml look-up-table file,
    gives a csv line with the requested descriptive statistic in 'ima' in each region 
    
    The ROI samplings are based on real-world coordinates, so the procedure is robust to
    bounding box and axis swap changes (unless parts of the brain are excluded from the box...).
    The affine matrices in all NIfTI headers thus NEED to be correct.
    
    ima: a NIfTI file to be analyzed.
    
    atlas_ima: the labelled image
    
    atlas_xml: the label information XML (see FSL atlases)
    
    stat: should be one of 'mean', 'median', 'sd', 'sum', 'volume','min'or 'max' (default: mean)
          Note that 'volume' is actually more than a simple voxel count times the
          voxel size,we use the sum of all voxel times the voxel size instead.
          So if you are not using jacobian maps or tissue probability maps, make
          sure that the supplied image contains zeros or ones and not actual 
          quantitative measurements. NaN are ignored.
    
    meas: is a string describing the quantitative measurement in the image
    
    mask, mask_thr: One may specify a mask to apply to the image (mask), 
                   with a threshold (mask_thr) if the mask image is not binary.
    
    out_fileroot: prefix for the output csv file
     
    jacobian: an optional jacobian modulation map to weight the average measurements with voxel volumes.
              Use only with 'mean' measurement.
    
    Warning: Does not take in account scl_slope scl_intercept, uses raw values.
    '''
    import nibabel as nb
    import pandas as pd
    import lxml.etree as et
    import numpy as np
    import os
    import re
    from ginnipi.toolbox.computations import ijk2xyz,xyz2ijk
           
    if stat in ['mean', 'median', 'sd', 'sum','volume','min','max']:
        # ...then we can get working
        # load atlas
        roi_atlas = nb.load(atlas_ima)
        rois = roi_atlas.get_data()
        
        # load the FSL XML Look-Up-Table
        fsl_xml = et.parse(atlas_xml)
        labels  = fsl_xml.findall(".//label")
        atlas_name = re.sub('[:;,$\s\-\(\)\[\]]','_',fsl_xml.findall('./header/name')[0].text)
        # load quantitative image
        img = nb.load(ima)
        data = img.get_data()  
        
        if mask:
            mask_ima = nb.load(mask)
            msk = mask_ima.get_data()
            
        # Load jacobian image if present
        if jacobian:
            jcb_ima = nb.load(jacobian)  
            jcb = jcb_ima.get_data()
            
        # Select statistic computation function
        if stat=='mean':
            comp = np.nanmean
        elif stat=='sd':
            comp = np.nanstd
        elif stat == 'median':
            comp = np.nanmedian
        elif stat == 'sum':
            comp = np.nansum
        elif stat == 'min':
            comp = np.nanmin
        elif stat == 'max':
            comp = np.nanmax
        elif stat == 'volume':
            # Jacobian determinants volumetry 
            # (jacobians are the in_file, they do not serve as weights in this case)
            comp = lambda x : np.nansum(x)*np.prod(img.header['pixdim'][1:4])
        
        # label per label
        if obs_name:
            values=[(atlas_name+'/'+stat+'/'+meas,obs_name)]
        else:
            values=[(atlas_name+'/'+stat+'/'+meas,ima)]
        for label in labels:
            print(label.text)
            if int(label.get('index')) > 0:
                ijk = np.nonzero(rois == int(label.get('index')))
                xyz = ijk2xyz(roi_atlas.get_affine(),ijk)
                if mask:
                    # mask needs not be in same orientation and shape as label
                    in_idx = msk[xyz2ijk(mask_ima.get_affine(),xyz)] > mask_thr
                    ijk=map(lambda x: x[in_idx],ijk)
                    xyz=ijk2xyz(roi_atlas.get_affine(),ijk)
                if ijk:
                    # we have data, now compute stat
                    # ignoring neuro or radio and shape problems by converting to and from real world coordinates
                    chunk = data[xyz2ijk(img.get_affine(),xyz)]
                    if stat=='mean' and jacobian:
                        # special case where we weight the mean by jacobian determinants
                        jcbians = jcb[xyz2ijk(jcb_ima.get_affine(),xyz)]
                        values.append( (label.text, [np.sum(np.multiply(chunk, jcbians))/np.sum(jcbians)]) ) 
                    else:
                        chunk=chunk[np.isfinite(chunk)]
                        values.append( (label.text, [np.float(comp(chunk))]) )
                else:
                    # No data, NaN
                    values.append(label.text,[np.NaN])      
       
        # Write csv table with pandas
        if obs_name:
            out_filename = out_fileroot + '_' + obs_name
        else:
            out_filename = out_fileroot
        if meas:
            out_filename = out_filename + '_' + meas
        out_filename = out_filename + '.csv'
        out_file=os.path.abspath(out_filename)
        df_out=pd.DataFrame.from_items(values)
        df_out.to_csv(out_file, index=False)
        return out_file
    else:
        print("Please supply a valid statistic name:  'mean', 'median', 'sd', 'sum','min','max', 'volume'")
    

     
# def parse_eprime2(in_file,
#                   ckeys=['TrialCondition'],
#                   bkey=None,
#                   per_event=False,
#                   ignore=None,
#                   infer_durations=False):
#     '''
#     Takes in an E-prime 2.x txt log-file, parses it and outputs nipype bunches
#     
#     Inputs:
#     ------
#         ckeys: a list of strings with condition keys (fields
#                of the log file containing condition information)
#         
#         bkey: a string with the key for the blocks
#         
#         in_file: the path to the eprime text log-file to parse
#         
#         per_event: event-level output
#         
#         ignore: do not use this condition key when reporting
#         
#         infer_durations: when durations are not present, infer duration 
#                          from onset time differences.
#     
#     Outputs:
#     --------
#         bunch: the nipype bunch object for fMRI modeling
#     '''
#                 
#     from ginnipi.toolbox.psych import PsychRun
#     
#     ps=PsychRun()
#     if bkey:
#         ps.fromEprime(in_file, 
#                       ckeys=ckeys,
#                       bkey=bkey)
#     else:
#         ps.fromEprime(in_file, 
#                       ckeys=ckeys)
# 
#     # Now we need to sort this huge mess of a log file
#     # In case durations were not available
#     if infer_durations:
#         ps.durationsFromOnsets()
#     
#     # events
#     ps.getEventOnsetDurationTable()
#     if per_event:
#         if ignore:
#             event_bunch = ps.getNipypeBunches('event',ignore=ignore)
#         else:
#             event_bunch = ps.getNipypeBunches('event')
#         return event_bunch
#    
#     if bkey: 
#         # block level
#         ps.getBlockOnsetDurationTable()
#         if ignore:
#             block_bunch= ps.getNipypeBunches('block',ignore=ignore)
#         else:
#             block_bunch= ps.getNipypeBunches('block')
#         return block_bunch
#     else:   
#         # trial level    
#         ps.getTrialOnsetDurationTable()
#         if ignore:
#             trial_bunch = ps.getNipypeBunches('trial',ignore=ignore)
#         else:
#             trial_bunch = ps.getNipypeBunches('trial')
#         return trial_bunch
    
    
def parse_eprime2(in_file,
                  out_filename='out_psychrun.pkl',
                  procedure=None,
                  onset_key=None,
                  duration_key=None,
                  RT_key=None,
                  response_key=None,
                  rttime_keys=None,
                  aliasing=None,
                  dirac=250,
                  ckeys=['TrialCondition'],
                  bkey=None,
                  per_event=False,
                  ignore=None,
                  infer_durations=False):
    '''
    Takes in an E-prime 2.x txt log-file, parses it and outputs nipype bunches
    
    Inputs:
    ------
        ckeys: a list of strings with condition keys (fields
               of the log file containing condition information)
        
        bkey: a string with the key for the blocks
        
        in_file: the path to the eprime text log-file to parse
        
        per_event: event-level output
        
        ignore: do not use this condition key when reporting
        
        infer_durations: when durations are not present, infer duration 
                         from onset time differences.
    
    Outputs:
    --------
        bunch: the nipype bunch object for fMRI modeling
    '''
                
    from ginnipi.toolbox.psych import PsychRun
    import os 
    
    ps=PsychRun()
    if procedure is not None:
        kwargs=dict(procedure=procedure,
                    onset_key=onset_key, 
                    duration_key=duration_key)
        
            
        if response_key is not None: 
            kwargs['response_key']=response_key
        if rttime_keys is not None:
            kwargs['rttime_keys']=rttime_keys 
            if dirac is not None:
                kwargs['dirac']=dirac
            else:
                kwargs['dirac']=250 # default value
        if RT_key is not None:
            kwargs['RT_key']=RT_key
    else:
        kwargs={}
    
    if bkey:
        ps.fromEprime(in_file, 
                      ckeys=ckeys,
                      bkey=bkey,
                      **kwargs)
    else:
        ps.fromEprime(in_file, 
                      ckeys=ckeys,
                      **kwargs) 

    # Now we need to sort this huge mess of a log file
    # In case durations were not available
    if infer_durations:
        ps.durationsFromOnsets()
    # save PsychRun object
    out_file=os.path.join(os.path.abspath('.'),out_filename)
    ps.save(out_file)
    # get bunches
    # events
    ps.getEventOnsetDurationTable()
    if per_event:
        if ignore:
            event_bunch = ps.getNipypeBunches('event',ignore=ignore)
        else:
            event_bunch = ps.getNipypeBunches('event')
        return event_bunch, out_file
   
    if bkey: 
        # block level
        ps.getBlockOnsetDurationTable()
        if ignore:
            block_bunch= ps.getNipypeBunches('block',ignore=ignore)
        else:
            block_bunch= ps.getNipypeBunches('block')
        return block_bunch, out_file
    else:   
        # trial level    
        ps.getTrialOnsetDurationTable()
        if ignore:
            trial_bunch = ps.getNipypeBunches('trial',ignore=ignore)
        else:
            trial_bunch = ps.getNipypeBunches('trial')
        return trial_bunch, out_file
    
def get_parse_rules(rules_xml, sequence):
    '''
    Convenience tool for psychological log file parsing
    Parses an omnibus xml rules files and looks for the 
    entry corresponding to the supplied eprime sequence 
    ID.
    
    Works with the parse_eprime2 function in this module
    
    'conditions keys' are the names of the log fields
    where to find the experimental condition information.
    
    Inputs:
    -------
        sequence: a string that should match a run_id attribute in the XML
        
        rules_xml: see below how the rules_xml file should be structured
    
    <?xml version="1.0" encoding="UTF-8"?>
    <experiment name="archeoneuro">
        <sequence run_id="VISULOCA" infer_durations="True">
            <condition_key>TrialCondition</condition_key>
            <condition_key ignore="True">Running</condition_key>
            <block_key>TrialCondition</block_key>
        </sequence>
            <sequence run_id="GRAVURE1">
            <condition_key>TrialCondition</condition_key>
        </sequence>
        <sequence run_id="TRACE1">
            <condition_key>EventCondition</condition_key>
        </sequence>
    </experiment>
    
        
    
    Outputs:
    --------
        ckeys: condition keys
        
        bkey: block key
        
        ignore: a condition key used for parsing, not reporting
                A single 'ignore' is allowed per sequance at the moment
                
        infer_durations: estimate durations from onset time differences.
                         Used when no duration was logged.
    '''
    import lxml.etree as et
    import os
    
    ignore=None
    if os.path.exists(rules_xml):
        try:
            xml = et.parse(rules_xml)
        except:
            print('Failed to parse the %s xml eprime-log parsing rules file:' % rules_xml)
            raise
    else:
        raise IOError('file %s does not exist, or we don''t have read permissions' % rules_xml)
    
    seq=xml.find("./sequence[@run_id='"+sequence+"']")
    if seq is None:
        raise KeyError("The xml rules file does not contain an entry for the eprime sequence %s" % sequence)
    else:
        # let's parse the information...
        # should we infer the durations?
        if 'infer_durations' in seq.attrib.keys():
            if seq.attrib['infer_durations']=='True':
                infer_durations=True
            else:
                infer_durations=False
        else:
            infer_durations=False
            
        if 'save' in seq.attrib.keys():
            if seq.attrib['save']=='True':
                to_save=True
            else:
                to_save=False
        else:
            to_save=False
            
        if 'managed' in seq.attrib.keys():
            if seq.attrib['managed']=='True':
                managed=True
            else:
                managed=False
        else:
            managed=False
            
        if 'per_event' in seq.attrib.keys():
            if seq.attrib['per_event']=='True':
                per_event=True
            else:
                per_event=False
        else:
            per_event=False
            
        # extract condition keys
        ckey_nodes=seq.findall('./condition_key')
        ckeys = [ckey_node.text for ckey_node in ckey_nodes]
        for ckey_node in ckey_nodes:
            if 'ignore' in ckey_node.attrib.keys():
                if ckey_node.attrib['ignore']=='True':
                    ignore=ckey_node.text
                    # one ignore only
                    break
        # block key
        bkey_node= seq.find('./block_key')
        if bkey_node is not None:
            bkey=bkey_node.text
        else:
            bkey=None
        
        out_dict={'ckeys':ckeys, 
                  'bkey':bkey, 
                  'ignore':ignore, 
                  'infer_durations':infer_durations,
                  'per_event':per_event}
        if to_save:
            try:
                out_dict['out_filename']=os.path.abspath(seq.find('./filename').text)
            except:
                out_dict['out_filename']=os.path.abspath('psych_run.pkl')
        
        if managed:
            procedure_node= seq.find('./procedure')
            out_dict['procedure']=procedure_node.text
            onset_node= seq.find('./onset_key')
            out_dict['onset_key']=onset_node.text
            duration_node= seq.find('./duration_key')
            out_dict['duration_key']=duration_node.text
            try:
                response_key_node= seq.find('./response_key')
                out_dict['response_key']=response_key_node.text
            except:
                pass
            try:
                RT_key_node= seq.find('./RT_key')
                out_dict['RT_key']=RT_key_node.text
            except:
                pass
            try:
                dirac_node= seq.find('./dirac')
                out_dict['dirac']=int(dirac_node.text())
            except:
                pass
            try:
                rttime_key_nodes= seq.findall('./rttime_key')
                if rttime_key_nodes:
                    out_dict['rttime_keys']=[rttime_key_node.text for rttime_key_node in rttime_key_nodes]
            except:
                pass            
    return out_dict
        
                
def combined_parse_eprime2(in_file, rules_xml, sequence): 
    '''
    combines get_parse_rules and parse_eprime2 functions from this module
    
    Inputs:
    -------
        in_file: the path to the eprime2 text log file to parse
        
        per_event: event-level output
        
        sequence: a string that should match a run_id attribute in the XML
        
        rules_xml: see below how the rules_xml file should be structured
    
    <?xml version="1.0" encoding="UTF-8"?>
    <experiment name="archeoneuro">
        <sequence run_id="VISULOCA" infer_durations="True">
            <condition_key>TrialCondition</condition_key>
            <condition_key ignore="True">Running</condition_key>
            <block_key>TrialCondition</block_key>
        </sequence>
            <sequence run_id="GRAVURE1">
            <condition_key>TrialCondition</condition_key>
        </sequence>
        <sequence run_id="TRACE1">
            <condition_key>EventCondition</condition_key>
        </sequence>
    </experiment>
    
    
    Outputs:
    --------
        bunch: the nipype bunch object for fMRI linear modeling   
    '''   
    from ginnipi.toolbox.computations import parse_eprime2, get_parse_rules
    
    kwargs = get_parse_rules(rules_xml, sequence)
    return parse_eprime2(in_file,**kwargs)

def psych_concatenate(in_list,out_filename):
    '''
    Loads and concatenate pickled PsychRun objects
    Saves a pickled concatenated Psych object for future use
    '''
    from ginnipi.toolbox.psych import PsychRun
    import os
    
    ps=PsychRun()
    ps.load(in_list[0])
    if len(in_list)>0:
        for eprime in in_list[1:]:
            other=PsychRun()
            other.load(eprime)
            ps.appendPsychRun(other)
    ps.save(os.path.abspath(out_filename))
    return os.path.abspath(out_filename)

def get_bunches(psych,mode='Trial',ignore=[],msec2sec=True):
    '''
    Gets a bunch from a psych object
    '''
    from ginnipi.toolbox.psych import PsychRun
    
    ps=PsychRun()
    ps.load(psych)
    return ps.getNipypeBunches(mode, 
                        ignore,
                        msec2sec) 

def get_trace_bunches(trace):
    '''
    Gets a bunch from a psych object
    Totally ad-hoc for archeoneuro
    '''
    from ginnipi.toolbox.psych import PsychRun
    import re
    from nipype.interfaces.base import Bunch
     
    
    ps=PsychRun()
    ps.load(trace)
    table=ps.getEventOnsetDurationTable(msec2sec=True)
    ponsets=[]
    pdurations=[]
    pconditions=[]
    ref_onsets=[]
    ref_durations=[]
    ref_levels=['vertical', 'nonvertical']
    ignore=['Img','Condition','Stimulation']
    for condition, levels in ps.conditions.iteritems():
        if condition not in ignore: 
            for level in levels:
                if level not in ref_levels:
                    pconditions.append(level)
                    try:
                        subset=table[table[condition].str.match('^'+level+'$') & table['name'].str.match('customEventReport')]
                    except:
                        table[condition]=table[condition].astype(str)
                        subset=table[table[condition].str.match('^'+str(level)+'$') & table['name'].str.match('customEventReport')]
                    ponsets.append(subset['onset'].tolist())
                    pdurations.append(subset['duration'].tolist())
                    imgs= map(lambda x: re.sub('[E]','',x), subset['Img'].tolist())
                    pref_onsets=[]
                    pref_durations=[]
                    for img in imgs:
                        ref_rows = table[table['Img'].str.match('C'+img)] 
                        pref_onsets.append(float(ref_rows['onset'].get_values()[0]))
                        pref_durations.append(float(ref_rows['duration'].get_values()[0]))
                    ref_onsets.append(pref_onsets)
                    ref_durations.append(pref_durations)
    man_made=['humain']
    onsets = [[],[],[],[]]
    durations = [[],[],[],[]]
    conditions = ['Attribution_H', 
                  'Attribution_NH', 
                  'Orientation_refH',
                  'Orientation_refNH']
    for idx,condition in enumerate(pconditions):
        if condition in man_made: 
            onsets[0]+=ponsets[idx]
            durations[0]+=pdurations[idx]
            onsets[2]+=ref_onsets[idx]
            durations[2]+=ref_durations[idx]
        else:
            onsets[1]+=ponsets[idx]
            durations[1]+=pdurations[idx]
            onsets[3]+=ref_onsets[idx]
            durations[3]+=ref_durations[idx]
                
    for name in list(set(table['name'].tolist())):
        if not name=='customEventReport':
            conditions.append(name)
            subset=table[table['name'].str.match('^'+name+'$')]
            onsets.append(subset['onset'].tolist())
            durations.append(subset['duration'].tolist())   
    out_bunch=Bunch(conditions=conditions,
                    onsets=onsets,
                    durations=durations)
    return out_bunch

def csv_append(csv_list, out_filename):
    '''
    Takes in a list of csv (w. similar formatting),
    outputs a merged csv (i.e. with rows concatenated).
    '''
    import pandas as pd
    import os
    

    merged_csv = pd.concat([pd.read_csv(to_append, delim_whitespace=True, header=None) for to_append in csv_list], axis=0)
    merged_csv.to_csv(os.path.abspath(out_filename), sep=' ', header=None, index=False)
    return os.path.abspath(out_filename)
    
    
    
    
    
    
