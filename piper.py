#!/usr/bin/python3

import argparse
import os
import subprocess

import pfpd_const as pfpd


######################################################################################################
"""WE are using SLURM workload manager. If you are using something else - change the commands below"""
######################################################################################################

SBATCH_PIPER = '#!/bin/sh\n' \
               '#SBATCH --ntasks=1\n' \
               '#SBATCH --time=20:00:00\n' \
               '#SBATCH --get-user-env\n' \
               + pfpd.PIPER_DOCKING
SBATCH_EXTRACT_TOP_DECOYS = "#!/bin/sh\n" \
                            "#SBATCH --nodes 1\n" \
                            "#SBATCH --ntasks=1\n" \
                            "#SBATCH --time=10:00:00\n" \
                            "#SBATCH --get-user-env\n" \
                            + pfpd.EXTRACT_PIPER_MODELS
SBATCH_PREP_FPD_INPUT = "#!/bin/sh\n" \
                        "#SBATCH --ntasks=1\n" \
                        "#SBATCH --time=10:00:00\n" \
                        "#SBATCH --get-user-env\n" \
                        + pfpd.PREPARE_FPD_IN

# 'dependency' in these commands will be changed to specific 'dependency=aftercorr$job_num' before run

RUN_PIPER = ['sbatch', '--mem=1500m', 'run_piper']
RUN_EXTRACT_DECOYS = ['sbatch', '--mem-per-cpu=1500m', 'run_extract_decoys']

#############################################################################################
"""Attention! Next 3 functions send all the jobs to cluster. They are also SLURM dependent"""
#############################################################################################


def create_batch(receptor, run, i):
    """Create batch scripts for jobs that will be sent to cluster, such as
    PIPER docking, models extraction, fpd input preparation and fpd refinement"""
    rec_name = os.path.join(piper_dir, receptor.lower() + '_nmin.pdb')
    lig_name = 'lig.' + "{:04}".format(i) + '_nmin.pdb'
    if run == 'piper':
        with open('run_piper', 'w') as piper:
            piper.write(SBATCH_PIPER.format(decoys=pfpd.N_ROTS, r=rec_name, l=lig_name))
    elif run == 'decoys':
        with open('run_extract_decoys', 'w') as extract_decoys:
            extract_decoys.write(SBATCH_EXTRACT_TOP_DECOYS % (pfpd.PIPER_MODELS_EXTRACTED, lig_name))


def send_piper_jobs(processed_receptor, n_frags):
    print("**************Sending PIPER jobs**************")
    receptor_name = os.path.splitext(os.path.basename(receptor_path))[0]
    array_for_flexpepdock = ''

    # PIPER step
    os.chdir(piper_dir)
    ppk_receptor = os.path.splitext(processed_receptor)[0] + '.ppk.pdb'
    for i in range(1, n_frags + 1):
        run_dir = "{:02}".format(i)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        os.chdir(run_dir)
        os.system(pfpd.COPY.format(os.path.join(piper_dir, 'ligands', 'lig.' + "{:04}".format(i) + '_nmin.pdb'),
                                   'lig.' + "{:04}".format(i) + '_nmin.pdb'))
        create_batch(receptor_name, 'piper', i)
        create_batch(receptor_name, 'decoys', i)
        create_batch(ppk_receptor, 'prepare_inputs', i)

        # run PIPER docking, extract top 250 decoys and prepare input for refinement
        run_piper = str(subprocess.check_output(RUN_PIPER))
        piper_id = ''.join(d for d in run_piper if d.isdigit())
        RUN_EXTRACT_DECOYS.insert(1, '--dependency=aftercorr:%s' % piper_id)
        run_extract_decoys = str(subprocess.check_output(RUN_EXTRACT_DECOYS))
        array_for_flexpepdock += ':' + ''.join(d for d in run_extract_decoys if d.isdigit())
        os.chdir(piper_dir)
    print(array_for_flexpepdock)


####################################################
"""ATTENTION!!!!! Do not change the code below!!!"""
####################################################


def process_for_piper(receptor):
    """Prepare inputs for piper run"""
    frags_list = []
    if not os.path.exists(piper_dir):  # Create directory
        os.makedirs(piper_dir)
    print("**************Preparing PIPER inputs**************")
    for frag in os.listdir(fixbb_dir):
        if os.path.splitext(frag)[1] == '.pdb':   # Rename chain ID to 'B'
            pfpd.rename_chain(os.path.join(fixbb_dir, frag), 'B')
            frags_list.append(os.path.basename(frag))
    n_frags = pfpd.count_pdbs(fixbb_dir)
    with open(os.path.join(piper_dir, 'runs_list'), 'w') as runs_list:  # Create list 1 - 50
        for i in range(1, n_frags + 1):
            runs_list.write("{:02}\n".format(i))

    # Create symlinks for piper run
    if not os.path.exists(ligands_dir):
        os.makedirs(ligands_dir)
    lig_inx = 1
    for frag_name in frags_list:
        os.system('ln -s {} {}/lig.{:04}.pdb'.format(os.path.join(fixbb_dir, frag_name),
                                                     ligands_dir, lig_inx))
        lig_inx += 1

    # process ligands for PIPER
    os.chdir(ligands_dir)
    for lig in os.listdir(ligands_dir):
        os.system(pfpd.PDBPREP.format(lig))
        os.system(pfpd.PDBNMD.format(lig))

    # prepare receptor for piper
    os.system(pfpd.COPY.format(receptor, piper_dir))
    os.chdir(piper_dir)

    os.system(pfpd.CLEAN_PDB.format(receptor, 'nochain'))  # Clean the receptor
    pfpd.rename_chain(os.path.basename(receptor)[:-4] + '_nochain.pdb', 'A')
    name_for_piper = os.path.basename(receptor).lower()
    os.rename(os.path.splitext(os.path.basename(receptor))[0] + '_nochain.pdb',
              name_for_piper)
    os.system(pfpd.PDBPREP.format(name_for_piper))
    os.system(pfpd.PDBNMD.format(name_for_piper))
    os.chdir(root)
    return n_frags


def arg_parser():
    parser = argparse.ArgumentParser(description='Please provide a pdb file for receptor '
                                                 'and a text file with peptide sequence (or FASTA file)\n')

    parser.add_argument('receptor')
    parser.add_argument('peptide_sequence')

    return parser


def run_protocol(receptor):

    # process ligands and receptor for piper run
    n_frags = process_for_piper(receptor)  # only basename

    # run piper docking, extract top 250 models from each run

    send_piper_jobs(os.path.basename(receptor).lower(), n_frags)

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

    # Define all the directories that will be created:
    root = os.getcwd()
    fragments_dir = os.path.join(root, 'top_50_frags')
    fixbb_dir = os.path.join(fragments_dir, 'fixbb')
    piper_dir = os.path.join(root, 'piper')
    ligands_dir = os.path.join(piper_dir, 'ligands')

    run_protocol(receptor_path)
