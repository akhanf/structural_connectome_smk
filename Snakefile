from os.path import join
from snakemake.utils import format
from bids import BIDSLayout
import pandas as pd

configfile: "config.yaml"

#load participants.tsv file, and strip of sub- from participant_id column
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.str.slice(4).to_list() 


# can also override at command-line with e.g.:  --config bids_dir='path/to/dir'  or --configfile ...
prepdwi_dir = join(config['prepdwi_dir'],'prepdwi')
fmriprep_dir = join(config['fmriprep_dir'],'fmriprep')


#set container for --use-singularity
container: config['singularity']

# this uses pybids to get info about the bids dataset (validate=False allows for non-validated tags)
layout = BIDSLayout(prepdwi_dir, validate=False)

#get entities from cfg file to select input files from prepdwi
entities = config['entities']
#subjects = layout.get_subjects(**entities)  #subjects from participants.tsv
sessions = layout.get_sessions(**entities)


ntracts_million = config['ntracts_million']

#create strings including wildcards for subj_sess_dir and subj_sess_prefix
if len(sessions) > 0:
    subj_sess_dir = join('sub-{subject}','ses-{session}')
    subj_sess_prefix = 'sub-{subject}_ses-{session}'
else:
    subj_sess_dir = 'sub-{subject}'
    subj_sess_prefix = 'sub-{subject}'


rule all:
    input:
        vis = expand(join('plots',f'{subj_sess_prefix}_connectome_vis.png'),subject=subjects,session=sessions),



include: 'rules/common.smk'
include: 'rules/mrtrix.smk'

