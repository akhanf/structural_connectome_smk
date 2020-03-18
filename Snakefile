from os.path import join
from snakemake.utils import format
from bids import BIDSLayout

configfile: "config.yaml"


# bids_dir and out_dir set in json file.
# can also override at command-line with e.g.:  --config bids_dir='path/to/dir'  or --configfile ...
bids_dir = config['bids_dir']
prepdwi_dir = join(config['prepdwi_dir'],'prepdwi')
fmriprep_dir = join(config['fmriprep_dir'],'fmriprep')
out_dir = config['out_dir']

work_dir = config['work_dir']

# this uses pybids to get info about the bids dataset (validate=False allows for non-validated tags)
layout = BIDSLayout(prepdwi_dir, validate=False)

#get entities from cfg file to select input files from prepdwi
entities = config['entities']
subjects = layout.get_subjects(**entities)
sessions = layout.get_sessions(**entities)

print(subjects)
print(sessions)


#create strings including wildcards for subj_sess_dir and subj_sess_prefix
if len(sessions) > 0:
    subj_sess_dir = join('sub-{subject}','ses-{session}')
    subj_sess_prefix = 'sub-{subject}_ses-{session}'
else:
    subj_sess_dir = 'sub-{subject}'
    subj_sess_prefix = 'sub-{subject}'

datatype = config['entities']['datatype']
space = config['entities']['space']
suffix = config['entities']['suffix']


rule all:
    input:
        tracts = expand(join(work_dir,subj_sess_dir,'tracts.tck'),subject=subjects,session=sessions)

 
rule convert_to_mif:
    input:
        dwi_nii = join(prepdwi_dir,subj_sess_dir,datatype,subj_sess_prefix + '_dwi_space-T1w_preproc.nii.gz'),
        dwi_bval = join(prepdwi_dir,subj_sess_dir,datatype,subj_sess_prefix + '_dwi_space-T1w_preproc.bval'),
        dwi_bvec = join(prepdwi_dir,subj_sess_dir,datatype,subj_sess_prefix + '_dwi_space-T1w_preproc.bvec')
    output:
        #use simplified naming for mrtrix work files
        dwi_mif = join(work_dir,subj_sess_dir,'dwi.mif'),
        b0_mif = join(work_dir,subj_sess_dir,'b0.mif')
    envmodules:
        'mrtrix'
    shell:
        "mrconvert {input.dwi_nii} {output.dwi_mif} -fslgrad {input.dwi_bvec} {input.dwi_bval} -datatype float32 -strides 0,0,0,1 &&"
        "dwiextract {output.dwi_mif} - -bzero | mrmath - mean {output.b0_mif} -axis 3"







rule gen_5tt:
    input: 
        lut_input = 'lut/FreeSurferColorLUT.txt',
        lut_output = 'lut/FreeSurfer2ACT.txt',
        aparcaseg_nii = join(fmriprep_dir,subj_sess_dir,'anat',subj_sess_prefix + '_desc-aparcaseg_dseg.nii.gz')
    output:
        indices_mif = join(work_dir,subj_sess_dir,'indices.mif'),
        cgm_mif = temp(join(work_dir,subj_sess_dir,'cgm.mif')),
        sgm_mif = temp(join(work_dir,subj_sess_dir,'sgm.mif')),
        wm_mif = temp(join(work_dir,subj_sess_dir,'wm.mif')),
        csf_mif = temp(join(work_dir,subj_sess_dir,'csf.mif')),
        path_mif = temp(join(work_dir,subj_sess_dir,'path.mif')),
        seg_5tt = join(work_dir,subj_sess_dir,'seg_5tt.mif')
    envmodules:
        'mrtrix'
    shell:
        "labelconvert {input.aparcaseg_nii} {input.lut_input} {input.lut_output} {output.indices_mif} && "
        "mrcalc {output.indices_mif} 1 -eq {output.cgm_mif} &&"
        "mrcalc {output.indices_mif} 2 -eq {output.sgm_mif} &&"
        "mrcalc {output.indices_mif} 3 -eq {output.wm_mif} &&"
        "mrcalc {output.indices_mif} 4 -eq {output.csf_mif} &&"
        "mrcalc {output.indices_mif} 5 -eq {output.path_mif} &&"
        "mrcat {output.cgm_mif} {output.sgm_mif} {output.wm_mif} {output.csf_mif} {output.path_mif} - -axis 3 | mrconvert - {output.seg_5tt} -datatype float32 && "
        "5ttcheck {output.seg_5tt}"

rule estimate_response:
    input:         
        dwi_mif = join(work_dir,subj_sess_dir,'dwi.mif'),
        seg_5tt = join(work_dir,subj_sess_dir,'seg_5tt.mif')
    params:
        tempdir = join(work_dir,subj_sess_dir,'tmp_dwi2response')
    #can add additional params here (lmax)
    output:
        rf_wm = join(work_dir,subj_sess_dir,'rf_wm.txt'),
        rf_gm = join(work_dir,subj_sess_dir,'rf_gm.txt'),
        rf_csf = join(work_dir,subj_sess_dir,'rf_csf.txt'),
        rf_voxels = join(work_dir,subj_sess_dir,'rf_voxels.mif')
    envmodules:
        'mrtrix'
    shell:
        "dwi2response -tempdir {params.tempdir} msmt_5tt {input.dwi_mif} {input.seg_5tt} {output.rf_wm} {output.rf_gm} {output.rf_csf} -voxels {output.rf_voxels}"
    
rule compute_fod:
    input:
        dwi_mif = join(work_dir,subj_sess_dir,'dwi.mif'),
        rf_wm = join(work_dir,subj_sess_dir,'rf_wm.txt'),
        rf_gm = join(work_dir,subj_sess_dir,'rf_gm.txt'),
        rf_csf = join(work_dir,subj_sess_dir,'rf_csf.txt'),
        mask_nii = join(prepdwi_dir,subj_sess_dir,datatype,subj_sess_prefix + '_dwi_space-T1w_brainmask.nii.gz')
    output:
        fod_wm = join(work_dir,subj_sess_dir,'fod_wm.mif'),
        fod_gm = join(work_dir,subj_sess_dir,'fod_gm.mif'),
        fod_csf = join(work_dir,subj_sess_dir,'fod_csf.mif')
    envmodules:
        'mrtrix'
    shell:
        "dwi2fod msmt_csd {input.dwi_mif} {input.rf_wm}  {output.fod_wm} {input.rf_gm}  {output.fod_gm} {input.rf_csf}  {output.fod_csf}"
       
rule gen_tracts:
    input:
        fod_wm = join(work_dir,subj_sess_dir,'fod_wm.mif'),
    output:
        tracts = join(work_dir,subj_sess_dir,'tracts.tck')
    envmodules:
        'mrtrix'
    shell:
        "tckgen {input.fod_wm} {output.tracts} -select config['num_tracts'] -seed_dynamic {input.fod_wm} -backtrack -crop_at_gmwmi -maxlength 250 -cutoff 0.06"
   

