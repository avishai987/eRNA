#old format:
# run with sbatch manually, without profile.
sbatch_params = "--mail-user=avishai.wizel@mail.huji.ac.il --export=ALL,conda_path=/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/bin/,lmod=/etc/profile.d/huji-lmod.sh"
#TODO: covert to use slurm profile


# ---------------------------------------------------------------------------- #
#                              GSE126074_SNARE_seq-  BABEL                     #
# ---------------------------------------------------------------------------- #
# use BABEL to predict chromatin accessability 

# create input
rule create_h5_files:
    input:
        script = home_dir +"/scripts/run_tsv_to_h5.sh",
        py_script = home_dir + "/scripts/tsv_to_h5.py",
        path_to_rna_tsv = home_dir + "/01_raw_data/GSE126074/GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv",
        path_to_atac_tsv =home_dir + "/01_raw_data/GSE126074/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv"
    shell:
        '''
        sbatch  {sbatch_params}   --output={home_dir}/logs/create_h5_files_%j.log  \
	    {input.script} {input.py_script} {input.path_to_rna_tsv} {input.path_to_atac_tsv}  {home_dir}/BABEL/01_h5_files/GSE126074_SNARE_seq.h5
        '''
        
        

       
rule run_babel_SNARE_seq:
    input:
        script = home_dir +"/scripts/run_babel.sh",
        py_script = home_dir + "/BABEL/babel_repo/bin/train_model.py",
        path_to_rna_atac_h5 = home_dir + "/BABEL/01_h5_files/GSE126074_SNARE_seq.h5"
    params:
        out_dir = home_dir + "/BABEL/02_run/"
    shell:
        '''
        sbatch  {sbatch_params}  --output={home_dir}/logs/{rule}_%j.log  \
	    {input.script} {input.py_script} {input.path_to_rna_atac_h5}   {input.params.out_dir}
        '''

