#!/usr/bin/python3

import argparse
import sys
import os
import subprocess

import pfpd_const as pfpd

######################################################################################################
"""WE are using SLURM workload manager. If you are using something else - change the commands below"""
######################################################################################################

SBATCH_PREP_FPD_INPUT = "#!/bin/sh\n" \
                        "#SBATCH --ntasks=1\n" \
                        "#SBATCH --time=10:00:00\n" \
                        "#SBATCH --get-user-env\n" \
                        + pfpd.PREPARE_FPD_IN
SBATCH_FPD = '#!/bin/bash\n' \
             '#SBATCH --ntasks=300\n' \
             '#SBATCH --time=50:00:00\n' \
             '#SBATCH --get-user-env\n' \
             '{fpd_version}'
SBATCH_EXTRACT_TOP_MODEL = '#!/bin/bash\n' \
                           '#SBATCH --ntasks=1\n' \
                           '#SBATCH --time=1:00:00\n' \
                           '#SBATCH --get-user-env\n' \
                           + pfpd.EXTRACT_MODEL
SBATCH_RESCORING = '#!/bin/bash\n' \
                   '#SBATCH --ntasks=300\n' \
                   '#SBATCH --time=30:00:00\n' \
                   '#SBATCH --get-user-env\n' \
                   + pfpd.RESCORING
SBATCH_CLUSTERING = '#!/bin/bash\n' \
                    '#SBATCH --ntasks=1\n' \
                    '#SBATCH --time=1:00:00\n' \
                    '#SBATCH --get-user-env\n' \
                    + pfpd.CLUSTERING

# 'dependency' in these commands will be changed to specific 'dependency=aftercorr$job_num' before run

RUN_PREP_FPD_INPUTS = ['sbatch', '--mem-per-cpu=1500m', 'run_prepare_fpd_inputs']
RUN_REFINEMENT = ['sbatch', '--mem-per-cpu=1600m', 'run_refinement']
RUN_EXTRACT_TOP_MODEL = ['sbatch', '--mem-per-cpu=1600m', 'extract_model']
RUN_RESCORING = ['sbatch', '--mem-per-cpu=1600m', 'rescoring']
RUN_CLUSTERING = ['sbatch', '--mem-per-cpu=1600m', 'run_clustering']

#############################################################################################
"""Attention! Next 2 functions send all the jobs to cluster. They are also SLURM dependent"""
#############################################################################################


def create_batch(receptor, run, **kwargs):
    """Create batch scripts for jobs that will be sent to cluster, such as
    PIPER docking, models extraction, fpd input preparation and fpd refinement"""
    if run == 'prepare_inputs':
        ppk_receptor = os.path.join(prepack_dir, receptor)
        with open('run_prepare_fpd_inputs', 'w') as inputs:
            inputs.write(SBATCH_PREP_FPD_INPUT % (ppk_receptor, refinement_dir, refinement_dir))
    elif run == 'refinement':
        with open(os.path.join(refinement_dir, 'run_refinement'), 'w') as refinement:
            if talaris:
                refinement.write(SBATCH_FPD.format(fpd_version=pfpd.FPD_TALARIS.format(flags='refine_flags')))
            else:
                refinement.write(SBATCH_FPD.format(fpd_version=pfpd.FPD.format(flags='refine_flags')))
    elif run == 'clustering':
        with open(os.path.join(clustering_dir, 'run_clustering'), 'w') as cluster:
            if not talaris:
                cluster.write(SBATCH_CLUSTERING.format(native=os.path.join(prepack_dir, 'start.pdb'),
                                                       decoys=os.path.join(refinement_dir, kwargs['silent']),
                                                       sc_func='ref2015'))
            else:
                cluster.write(SBATCH_CLUSTERING.format(native=os.path.join(prepack_dir, 'start.pdb'),
                                                       decoys=os.path.join(refinement_dir, kwargs['silent']),
                                                       sc_func='talaris14'))


def run_fpd_jobs(receptor, native_path, silent_f):
    """Run FPD and clustering"""
# Prepare inputs
    receptor_name = os.path.basename(receptor).lower()
    ppk_receptor = os.path.splitext(receptor_name)[0] + '.ppk.pdb'
    jobs_list = []
    frags_num = pfpd.count_dirs(piper_dir) - 1  # There is also 'ligands' directory
    os.chdir(piper_dir)
    for i in range(1, frags_num + 1):
        run_dir = "{:02}".format(i)
        os.chdir(run_dir)
        create_batch(ppk_receptor, 'prepare_inputs')
        RUN_PREP_FPD_INPUTS.insert(1, '--dependency=aftercorr%s' % after_job)
        run_prepare_fpd_inputs = str(subprocess.check_output(RUN_PREP_FPD_INPUTS))
        jobs_list.append(''.join(d for d in run_prepare_fpd_inputs if d.isdigit()))
        os.chdir(piper_dir)

# FPD refinement
    if not os.path.exists(refinement_dir):
        os.makedirs(refinement_dir)
    os.chdir(refinement_dir)
    create_batch(receptor, 'refinement')
    refine_flags_file(native_path)
    all_prep_jobs = ''
    for job_id in jobs_list:
        all_prep_jobs += ':' + job_id  # there are multiple jobs that need to be separated by semicolon
    RUN_REFINEMENT.insert(1, '--dependency=aftercorr%s' % all_prep_jobs)
    run_refinement = str(subprocess.check_output(RUN_REFINEMENT))

    refinement_id = ''.join(d for d in run_refinement if d.isdigit())

    # Rescoring and Clustering
    if not os.path.exists(clustering_dir):
        os.makedirs(clustering_dir)
    create_batch(receptor, 'clustering', silent=silent_f)
    if not native:  # If native structure was not provided, top scoring structure will be taken as native for rescoring
        with open(os.path.join(refinement_dir, 'extract_model'), 'w') as extract_pdb:
            extract_pdb.write(SBATCH_EXTRACT_TOP_MODEL)
        with open(os.path.join(refinement_dir, 'rescoring'), 'w') as rescore:
            if talaris:
                rescore.write(SBATCH_RESCORING.format(sc_func='talaris14'))
            else:
                rescore.write(SBATCH_RESCORING.format(sc_func='ref2015'))
        RUN_EXTRACT_TOP_MODEL.insert(1, '--dependency=aftercorr:%s' % refinement_id)
        extract_model = str(subprocess.check_output(RUN_EXTRACT_TOP_MODEL))

        extract_model_id = ''.join(d for d in extract_model if d.isdigit())
        RUN_RESCORING.insert(1, '--dependency=aftercorr:%s' % extract_model_id)
        rescoring = str(subprocess.check_output(RUN_RESCORING))

        rescoring_id = ''.join(d for d in rescoring if d.isdigit())
        os.chdir(clustering_dir)
        RUN_CLUSTERING.insert(1, '--dependency=aftercorr:%s' % rescoring_id)
        subprocess.call(RUN_CLUSTERING)
    else:
        os.chdir(clustering_dir)
        RUN_CLUSTERING.insert(1, '--dependency=aftercorr:%s' % refinement_id)
        subprocess.call(RUN_CLUSTERING)
    os.chdir(root)


######################################################################
"""It is not recommended to change the flags, unless you know what you 
 are doing. If you do, flags can be changed in the next 2 functions"""
######################################################################


def prepack_flags_file(receptor):
    with open('prepack_flags', 'w') as flags:
        flags.write('-s start.pdb\n'
                    '-out:pdb\n'
                    '-scorefile ppk.score.sc\n'
                    '-nstruct 1\n'
                    '-flexpep_prepack\n'
                    '-ex1\n'
                    '-ex2aro\n'
                    '-use_input_sc\n'
                    '-unboundrot ' + receptor + '\n'
                    '-mute protocols.moves.RigidBodyMover\n'
                    '-mute core.chemical\n'
                    '-mute core.scoring.etable\n'
                    '-mute protocols.evalution\n'
                    '-mute core.pack.rotamer_trials\n'
                    '-mute protocols.abinitio.FragmentMover\n'
                    '-mute core.fragment\n'
                    '-mute protocols.jd2.PDBJobInputter')


def refine_flags_file(native_path):
    with open(os.path.join(refinement_dir, 'refine_flags'), 'w') as flags:
        flags.write('-in:file:l input_list\n'
                    '-scorefile score.sc\n'
                    '-out:pdb_gz\n'
                    '-out:file:silent_struct_type binary\n'
                    '-out:file:silent decoys.silent\n'
                    '-lowres_preoptimize\n'
                    '-flexPepDocking:pep_refine\n'
                    '-flexPepDocking:flexpep_score_only\n'
                    '-ex1\n'
                    '-ex2aro\n'
                    '-use_input_sc\n'
                    '-unboundrot {receptor}\n'
                    '-mute protocols.moves.RigidBodyMover\n'
                    '-mute core.chemical\n'
                    '-mute core.scoring.etable\n'
                    '-mute protocols.evalution\n'
                    '-mute core.pack.rotamer_trials\n'
                    '-mute protocols.abinitio.FragmentMover\n'
                    '-mute core.fragment\n'
                    '-mute protocols.jd2.PDBJobInputter'.format(receptor=receptor_path))
        if minimization:
            flags.write('\n-min_receptor_bb')
        if native:
            flags.write('\n-native {native}'.format(native=native_path))

####################################################
"""ATTENTION!!!!! Do not change the code below!!!"""
####################################################


def check_native_structure(native_path):
    """ Check the length of native struct with the length of the receptor + pep"""
    with open(native_path, 'r') as n:
        native_calphas = 0
        # native_sequence = ''
        for line in n:
            if line[13:15] == 'CA':
                native_calphas += 1
            #   native_sequence += line[17:20]  # TODO: you are not comparing the sequence !!!
    with open(receptor_path, 'r') as rec:
        rec_calphas = 0
        for line in rec:
            if line[13:15] == 'CA':
                rec_calphas += 1
        complex_calphas = rec_calphas + pep_length
    if native_calphas != complex_calphas:
        print(pfpd.BAD_NATIVE)
        sys.exit()
    else:
        return


def build_peptide(pep):
    """Create directory for prepacking and, extended peptide and change its chain id to 'B'"""
    print("Building extended peptide for prepacking")
    if not os.path.exists(prepack_dir):
        os.makedirs(prepack_dir)
    os.chdir(prepack_dir)
    # Build extended peptide
    os.system(pfpd.BUILD_PEPTIDE.format(pep))

    # Change chain ID to 'B'
    pfpd.rename_chain('peptide.pdb', 'B')
    os.chdir(root)


def combine_receptor_peptide(receptor, ligand, out):
    """Put together receptor and ligand"""
    with open(out, 'w') as combined:
        with open(receptor, 'r') as rcptr:
            for line in rcptr:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                    combined.write(line)
        with open(ligand, 'r') as lig:
            for line in lig:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                    combined.write(line)


def prepack_receptor(processed_receptor):
    """Prepack receptor for FlexPepDock run"""
    print("**************Prepacking receptor**************")
    os.chdir(prepack_dir)
    ppk_receptor = os.path.splitext(processed_receptor)[0] + '.ppk.pdb'
    receptor = os.path.join(piper_dir, processed_receptor)
    combine_receptor_peptide(receptor, 'peptide.pdb', 'start.pdb')
    prepack_flags_file(receptor)
    if os.path.isfile(ppk_receptor):            # No need to prepack the same receptor again
        print("Prepacked receptor already exists")
        return
    if talaris:
        os.system(pfpd.PREPACK_TALARIS)
    else:
        os.system(pfpd.PREPACK)
    with open(ppk_receptor, 'w') as renamed_rec:
        with open('start_0001.pdb', 'r') as start:
            for line in start:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                    if line[21] == 'A':
                        renamed_rec.write(line)
                else:
                    continue
    print("Prepack done!")


def arg_parser():
    parser = argparse.ArgumentParser(description='You need to provide a pdb file for receptor '
                                                 'and a text file with peptide sequence (or FASTA file)\n'
                                                 '\nIf you want to run it with talaris2014, add '
                                                 '"--restore_talaris_behavior" option and make sure you have'
                                                 ' both - 2016 and current versions of Rosetta.\n'
                                                 'For running FlexPepDock with receptor minimization add '
                                                 '"--receptor_min" flag.\nIf available, provide a complex structure '
                                                 'with "--native" flag for RMSD calculation in FlexPepDock step.')

    parser.add_argument('receptor')
    parser.add_argument('peptide_sequence')
    parser.add_argument('--restore_talaris_behavior', dest='talaris', action='store_true', default=False)
    parser.add_argument('--receptor_min', dest='minimize_receptor', action='store_true', default=False)
    parser.add_argument('--native', dest='native_structure', default=None)
    parser.add_argument('--after', dest='last_job', default=None)

    return parser


def run_protocol(receptor):

    if native:
        native_path = os.path.abspath(native)
        check_native_structure(native_path)
        silent_file = 'decoys.silent'
    else:
        native_path = None
        silent_file = 'decoys_rescored.silent'

    build_peptide(os.path.abspath(sys.argv[2]))  # build extended peptide and rename it's chain id to 'B'

    prepack_receptor(os.path.basename(receptor).lower())

    # run FlexPepDock, clustering, rescoring (if not native) and get top 10 models

    run_fpd_jobs(receptor, native_path, silent_file)


if __name__ == "__main__":

    arguments = arg_parser().parse_args()

    with open(arguments.peptide_sequence, 'r') as peptide:
        peptide_seq = peptide.readlines()
        if peptide_seq[0][0] == '>':
            peptide_seq = peptide_seq[1].strip()
        else:
            peptide_seq = peptide_seq[0].strip()

    pep_length = len(peptide_seq)
    receptor_path = os.path.abspath(arguments.receptor)

    talaris = arguments.talaris
    minimization = arguments.minimize_receptor
    native = arguments.native_structure
    after_job = arguments.last_job

    # Define all the directories that will be created:
    root = os.getcwd()

    piper_dir = os.path.join(root, 'piper')
    prepack_dir = os.path.join(root, 'prepacking')
    refinement_dir = os.path.join(root, 'refinement')
    clustering_dir = os.path.join(refinement_dir, 'clustering')

    final_dir = 'FINAL_RESULTS'  # top 10 models and score file

    run_protocol(receptor_path)
