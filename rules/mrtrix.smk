rule convert_to_mif:
    input:
        dwi_nii = join(prepdwi_dir,subj_sess_dir,datatype,f'{subj_sess_prefix}_dwi_space-T1w_preproc.nii.gz'),
        dwi_bval = join(prepdwi_dir,subj_sess_dir,datatype,f'{subj_sess_prefix}_dwi_space-T1w_preproc.bval'),
        dwi_bvec = join(prepdwi_dir,subj_sess_dir,datatype,f'{subj_sess_prefix}_dwi_space-T1w_preproc.bvec'),
        mask_nii = join(prepdwi_dir,subj_sess_dir,datatype,f'{subj_sess_prefix}_dwi_space-T1w_brainmask.nii.gz')
    output:
        #use simplified naming for mrtrix work files
        dwi_mif = join(work_dir,subj_sess_dir,'dwi.mif'),
        b0_mif = join(work_dir,subj_sess_dir,'b0.mif'),
        mask_mif = join(work_dir,subj_sess_dir,'mask.mif')
    log:
        f'logs/convert_to_mif/{subj_sess_prefix}.log'
    shell:
        "(mrconvert {input.dwi_nii} {output.dwi_mif} -fslgrad {input.dwi_bvec} {input.dwi_bval} -datatype float32 -strides 0,0,0,1 &&"
        "mrconvert {input.mask_nii} {output.mask_mif} &&"
        "dwiextract {output.dwi_mif} - -bzero | mrmath - mean {output.b0_mif} -axis 3) &> {log}"




rule gen_5tt:
    input: 
        lut_input = 'lut/FreeSurferColorLUT.txt',
        lut_output = 'lut/FreeSurfer2ACT.txt',
        aparcaseg_nii = join(fmriprep_dir,subj_sess_dir,'anat',f'{subj_sess_prefix}_desc-aparcaseg_dseg.nii.gz')
    output:
        indices_mif = join(work_dir,subj_sess_dir,'indices.mif'),
        cgm_mif = temp(join(work_dir,subj_sess_dir,'cgm.mif')),
        sgm_mif = temp(join(work_dir,subj_sess_dir,'sgm.mif')),
        wm_mif = temp(join(work_dir,subj_sess_dir,'wm.mif')),
        csf_mif = temp(join(work_dir,subj_sess_dir,'csf.mif')),
        path_mif = temp(join(work_dir,subj_sess_dir,'path.mif')),
        seg_5tt = join(work_dir,subj_sess_dir,'seg_5tt.mif')
    log:
        f'logs/gen_5tt/{subj_sess_prefix}.log'
    shell:
        "(labelconvert {input.aparcaseg_nii} {input.lut_input} {input.lut_output} {output.indices_mif} && "
        "mrcalc {output.indices_mif} 1 -eq {output.cgm_mif} &&"
        "mrcalc {output.indices_mif} 2 -eq {output.sgm_mif} &&"
        "mrcalc {output.indices_mif} 3 -eq {output.wm_mif} &&"
        "mrcalc {output.indices_mif} 4 -eq {output.csf_mif} &&"
        "mrcalc {output.indices_mif} 5 -eq {output.path_mif} &&"
        "mrcat {output.cgm_mif} {output.sgm_mif} {output.wm_mif} {output.csf_mif} {output.path_mif} - -axis 3 | mrconvert - {output.seg_5tt} -datatype float32 && "
        "5ttcheck {output.seg_5tt}) &> {log}"

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
    threads: 8
    log:
        f'logs/estimate_response/{subj_sess_prefix}.log'

    shell:
        "dwi2response -tempdir {params.tempdir} msmt_5tt {input.dwi_mif} {input.seg_5tt} {output.rf_wm} {output.rf_gm} {output.rf_csf} -voxels {output.rf_voxels} -nthreads {threads} &> {log}"
    
rule compute_fod:
    input:
        dwi_mif = join(work_dir,subj_sess_dir,'dwi.mif'),
        rf_wm = join(work_dir,subj_sess_dir,'rf_wm.txt'),
        rf_gm = join(work_dir,subj_sess_dir,'rf_gm.txt'),
        rf_csf = join(work_dir,subj_sess_dir,'rf_csf.txt'),
        mask_mif = join(work_dir,subj_sess_dir,'mask.mif')
    output:
        fod_wm = join(work_dir,subj_sess_dir,'fod_wm.mif'),
        fod_gm = join(work_dir,subj_sess_dir,'fod_gm.mif'),
        fod_csf = join(work_dir,subj_sess_dir,'fod_csf.mif')
    resources:
        time = 6*60, #in minutes
    threads: 8
    log:
        f'logs/compute_fod/{subj_sess_prefix}.log'
    shell:
        "dwi2fod -mask {input.mask_mif} msmt_csd {input.dwi_mif} {input.rf_wm}  {output.fod_wm} {input.rf_gm}  {output.fod_gm} {input.rf_csf}  {output.fod_csf} &> {log}"

rule mtnormalise:
    input:
        fod_wm = join(work_dir,subj_sess_dir,'fod_wm.mif'),
        fod_gm = join(work_dir,subj_sess_dir,'fod_gm.mif'),
        fod_csf = join(work_dir,subj_sess_dir,'fod_csf.mif'),
        mask_mif = join(work_dir,subj_sess_dir,'mask.mif')
    output:
        fodn_wm = join(work_dir,subj_sess_dir,'fodn_wm.mif'),
        fodn_gm = join(work_dir,subj_sess_dir,'fodn_gm.mif'),
        fodn_csf = join(work_dir,subj_sess_dir,'fodn_csf.mif')
    log:
        f'logs/mtnormalise/{subj_sess_prefix}.log'
    shell:
        "mtnormalise -mask {input.mask_mif} {input.fod_wm} {output.fodn_wm}  {input.fod_gm} {output.fodn_gm}  {input.fod_csf} {output.fodn_csf}  &> {log}"

rule gen_tracts:
    input:
        fodn_wm = join(work_dir,subj_sess_dir,'fodn_wm.mif'),
    output:
        tracts = join(work_dir,subj_sess_dir,'tracts.tck')
    params:
        ntracts = f'{ntracts_million}M'
    resources:
        time = get_walltime_tracts, #  24*60, #in minutes
        mem_mb = 4000 #in mb
    threads: 8 # num cores
    log:
        f'logs/gen_tracts/{subj_sess_prefix}.log'
    shell:
        "tckgen {input.fodn_wm} {output.tracts} -select {params.ntracts} -seed_dynamic {input.fodn_wm} -backtrack -crop_at_gmwmi -maxlength 250 -cutoff 0.06 -nthreads {threads} &> {log}"
   

rule run_sift2:
    input:
        tracts = join(work_dir,subj_sess_dir,'tracts.tck'),
        fodn_wm = join(work_dir,subj_sess_dir,'fodn_wm.mif'),
        seg_5tt = join(work_dir,subj_sess_dir,'seg_5tt.mif')
       
    output:
        weights = join(work_dir,subj_sess_dir,'sift2_weights.txt')

    threads: 8 # num cores
    log:
        f'logs/run_sift2/{subj_sess_prefix}.log'
    shell:
        "tcksift2 -act {input.seg_5tt} {input.tracts} {input.fodn_wm} {output.weights} -nthreads {threads} &> {log}"


rule convert_atlas_labels:
    input:
        aparcaseg_nii = join(fmriprep_dir,subj_sess_dir,'anat',f'{subj_sess_prefix}_desc-aparcaseg_dseg.nii.gz'),
        lut_in = 'lut/FreeSurferColorLUT.txt',
        lut_out = 'lut/fs_default.txt'
    output:
        atlas_labels = join(work_dir,subj_sess_dir,'atlas_aparcaseg.nii.gz')
    log:
        f'logs/convert_atlas_labels/{subj_sess_prefix}.log'
    shell:
        "labelconvert {input.aparcaseg_nii} {input.lut_in} {input.lut_out} {output.atlas_labels} &> {log}"


rule gen_connectome:
    input:  
        tracts = join(work_dir,subj_sess_dir,'tracts.tck'),
        atlas_labels = join(work_dir,subj_sess_dir,'atlas_aparcaseg.nii.gz'),
        weights = join(work_dir,subj_sess_dir,'sift2_weights.txt')
    output:
        connectome = join(work_dir,subj_sess_dir,'connectome_aparcaseg.csv')
    log:
        f'logs/gen_connectome/{subj_sess_prefix}.log'
    shell:
        "tck2connectome -tck_weights_in {input.weights} {input.tracts} {input.atlas_labels} {output.connectome} &> {log}"

rule visualize_connectome:
    input:
        connectome = join(work_dir,subj_sess_dir,'connectome_aparcaseg.csv')
    output:
        report(join('plots',f'{subj_sess_prefix}_connectome_vis.png'),caption='report/connectome_vis.rst',category='Connectome Visualization')
    log:
        f'logs/visualize_connectome/{subj_sess_prefix}.log'
    notebook:
        'notebooks/vis_connectome.ipynb'


