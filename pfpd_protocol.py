#!/usr/bin/python3

import argparse
import sys
import os
import subprocess

#########################
""" Change paths here """
#########################

ROSETTA_DIR = '/vol/ek/Home/alisa/rosetta/Rosetta/'
ROSETTA_2016_DIR = '/vol/ek/share/rosetta/rosetta_src_2016.20.58704_bundle/'

ROSETTA_DB = os.path.join(ROSETTA_DIR, 'main/database')
ROSETTA_2016_DB = os.path.join(ROSETTA_2016_DIR, 'main/database')

ROSETTA_BIN = os.path.join(ROSETTA_DIR, 'main/source/bin/')
ROSETTA_2016_BIN = os.path.join(ROSETTA_2016_DIR, 'main/source/bin/')

ROSETTA_TOOLS = os.path.join(ROSETTA_DIR, 'tools/')

# Path to this script (and also make_fragments.pl, clustering.py)
# Change paths also in make_fragments.pl script itself
PFPD_SCRIPTS = '/vol/ek/Home/alisa/scripts/piper-fpd/'

PIPER_DIR = '/vol/ek/Home/alisa/PIPER/'
PIPER_BIN = os.path.join(PIPER_DIR, 'bin/')

######################################
"""DO NOT change following commands"""
######################################

# Commands (ROSETTA)

CLEAN_PDB = os.path.join(ROSETTA_TOOLS, 'protein_tools/scripts/clean_pdb.py >> pdb_log ') + ' {} {}'

FIXBB_JD3 = 'mpirun -n 6 ' + os.path.join(ROSETTA_BIN, 'fixbb_jd3.mpiserialization.linuxgccrelease') + \
            ' -database ' + ROSETTA_DB + ' -in:file:job_definition_file {} > fixbb.log'
FIXBB_JD3_TALARIS = 'mpirun -n 6 ' + os.path.join(ROSETTA_BIN, 'fixbb_jd3.mpiserialization.linuxgccrelease') + \
                    ' -database ' + ROSETTA_DB + ' -restore_talaris_behavior ' \
                     '-in:file:job_definition_file {} > fixbb.log'
FIXBB = os.path.join(ROSETTA_BIN, 'fixbb.linuxgccrelease') + ' -database ' + ROSETTA_DB + \
        ' -in:file:s {frag} -resfile {resfile} -ex1 -ex2 -use_input_sc -scorefile design_score.sc >design.log'
FIXBB_TALARIS = os.path.join(ROSETTA_2016_BIN, 'fixbb.linuxgccrelease') + ' -database ' + ROSETTA_2016_DB + \
                ' -in:file:s {frag} -resfile {resfile} -ex1 -ex2 -use_input_sc -scorefile design_score.sc >design.log'

BUILD_PEPTIDE = os.path.join(ROSETTA_BIN, 'BuildPeptide.linuxgccrelease') + ' -in:file:fasta {}' \
                ' -database ' + ROSETTA_DB + ' -out:file:o peptide.pdb > build_peptide.log'

MAKE_FRAGMENTS = 'perl ' + os.path.join(PFPD_SCRIPTS, 'make_fragments.pl') + \
                 ' -verbose -id xxxxx {} 2>log'

FRAG_PICKER = os.path.join(ROSETTA_BIN, 'fragment_picker.linuxgccrelease') + \
              ' -database ' + ROSETTA_DB + ' @flags >makeFrags.log'

PREPACK = os.path.join(ROSETTA_BIN, 'FlexPepDocking.default.linuxgccrelease') + \
          ' -database ' + ROSETTA_DB + ' @prepack_flags >ppk.log'
PREPACK_TALARIS = 'mpirun -n 5 ' + os.path.join(ROSETTA_2016_BIN, 'FlexPepDocking.mpi.linuxgccrelease') + \
                  ' -database ' + ROSETTA_2016_DB + ' @prepack_flags >ppk.log'

FPD = 'ls *gz >input_list\n' \
      'mpirun ' + os.path.join(ROSETTA_BIN, 'FlexPepDocking.mpiserialization.linuxgccrelease') + \
      ' -database ' + ROSETTA_DB + ' @{flags} >>refinement_log'

FPD_TALARIS = 'ls *gz >input_list\n' \
              'mpirun ' + os.path.join(ROSETTA_2016_BIN, 'FlexPepDocking.mpi.linuxgccrelease') + \
              ' -database ' + ROSETTA_2016_DB + ' @{flags} >>refinement_log'
EXTRACT_MODEL = 'python ' + PFPD_SCRIPTS + 'extract_top_model.py'
RESCORING = 'python ' + PFPD_SCRIPTS + 'rescoring.py {sc_func}'
CLUSTERING = 'python ' + PFPD_SCRIPTS + 'clustering.py 2.0 {native} {decoys}'

# Commands (PIPER)

PDBPREP = 'perl ' + PIPER_BIN + 'phplibbin/pdbprep.pl {}'
PDBNMD = 'perl ' + PIPER_BIN + 'phplibbin/pdbnmd.pl "{}"' + ' ?'

PIPER_DOCKING = PIPER_BIN + 'piper -vv -c1.0 -k4 --msur_k=1.0' \
                ' --maskr=1.0 -T FFTW_EXHAUSTIVE -R {decoys} -t 1 -p ' + PIPER_DIR + \
                'prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec -f ' + PIPER_DIR + \
                'prms/coeffs.0.0.6.motif -r ' + PIPER_DIR + 'prms/rot70k.0.0.4.prm ' \
                                                            '{r} {l} >piper.log'

EXTRACT_PIPER_MODELS = "for f in `awk '{print $1}' ft.000.00 | head -%s`;" \
                       "do if [ ! -f {f}.pdb ]; then " + PIPER_DIR + "apply_ftresult.py " \
                       "-i $f ft.000.00 " + PIPER_DIR + "prms/rot70k.0.0.4.prm %s " \
                       "--out-prefix $f;fi;done"
APPLY_FTRESULTS = 'python ' + PIPER_DIR + 'apply_ftresult.py -i {model} ft.000.00 ' \
                  + PIPER_DIR + 'prms/rot70k.0.0.4.prm {lig} --out-prefix {out}'

# Other commands

COPY = 'cp {} {}'

PREPARE_FPD_IN = "piper_run=`pwd | awk -F'/' '{print $NF}'`\nfor f in `ls [0-9]*.pdb`;" \
                 "do cat %s ${f} | grep ATOM >%s/${piper_run}.${f};gzip %s/${piper_run}.${f};done"

# Other constants

NUM_OF_FRAGS = 50  # should be 50
PIPER_MODELS_EXTRACTED = str(250)  # should be 250
N_ROTS = str(70000)  # should be 70000

FRAGS_FILE = 'frags.100.{}mers'
THREE_TO_ONE_AA = {'G': 'GLY',
                   'A': 'ALA',
                   'V': 'VAL',
                   'L': 'LEU',
                   'I': 'ILE',
                   'P': 'PRO',
                   'C': 'CYS',
                   'M': 'MET',
                   'H': 'HIS',
                   'F': 'PHE',
                   'Y': 'TYR',
                   'W': 'TRP',
                   'N': 'ASN',
                   'Q': 'GLN',
                   'S': 'SER',
                   'T': 'THR',
                   'K': 'LYS',
                   'R': 'ARG',
                   'D': 'ASP',
                   'E': 'GLU'}

PSIPRED_OUTPUT = ['C', 'E', 'H']

PDB = '{}_{}.pdb'
FASTA = '{}_{}.fasta'
BAD_NATIVE = "The native structure and the receptor should have the same sequence and the same length.\n" \
             "Please provide a valid native structure"

######################################################################################################
"""WE are using SLURM workload manager. If you are using something else - change the commands below"""
######################################################################################################

SBATCH_PIPER = '#!/bin/sh\n' \
               '#SBATCH --ntasks=1\n' \
               '#SBATCH --time=20:00:00\n' \
               '#SBATCH --get-user-env\n'\
               + PIPER_DOCKING
SBATCH_EXTRACT_TOP_DECOYS = "#!/bin/sh\n" \
                            "#SBATCH --nodes 1\n" \
                            "#SBATCH --ntasks=1\n" \
                            "#SBATCH --time=10:00:00\n" \
                            "#SBATCH --get-user-env\n" \
                            + EXTRACT_PIPER_MODELS
SBATCH_PREP_FPD_INPUT = "#!/bin/sh\n" \
                        "#SBATCH --ntasks=1\n" \
                        "#SBATCH --time=10:00:00\n" \
                        "#SBATCH --get-user-env\n" \
                        + PREPARE_FPD_IN
SBATCH_FPD = '#!/bin/bash\n' \
             '#SBATCH --ntasks=300\n' \
             '#SBATCH --time=50:00:00\n' \
             '#SBATCH --get-user-env\n' \
             '{fpd_version}'
SBATCH_EXTRACT_TOP_MODEL = '#!/bin/bash\n' \
                           '#SBATCH --ntasks=1\n' \
                           '#SBATCH --time=1:00:00\n' \
                           '#SBATCH --get-user-env\n' \
                           + EXTRACT_MODEL
SBATCH_RESCORING = '#!/bin/bash\n' \
                   '#SBATCH --ntasks=300\n' \
                   '#SBATCH --time=30:00:00\n' \
                   '#SBATCH --get-user-env\n' \
                   + RESCORING
SBATCH_CLUSTERING = '#!/bin/bash\n' \
                    '#SBATCH --ntasks=1\n' \
                    '#SBATCH --time=1:00:00\n' \
                    '#SBATCH --get-user-env\n' \
                    + CLUSTERING

# 'dependency' in these commands will be changed to specific 'dependency=aftercorr$job_num' before run

RUN_PIPER = ['sbatch', '--mem=1500m', 'run_piper']
RUN_EXTRACT_DECOYS = ['sbatch', 'dependency', '--mem-per-cpu=1500m', 'run_extract_decoys']
RUN_PREP_FPD_INPUTS = ['sbatch', 'dependency', '--mem-per-cpu=1500m', 'run_prepare_fpd_inputs']
RUN_REFINEMENT = ['sbatch', 'dependency', '--mem-per-cpu=1600m', 'run_refinement']
RUN_EXTRACT_TOP_MODEL = ['sbatch', 'dependency', '--mem-per-cpu=1600m', 'extract_model']
RUN_RESCORING = ['sbatch', 'dependency', '--mem-per-cpu=1600m', 'rescoring']
RUN_CLUSTERING = ['sbatch', 'dependency', '--mem-per-cpu=1600m', 'run_clustering']

#############################################################################################
"""Attention! Next 3 functions send all the jobs to cluster. They are also SLURM dependent"""
#############################################################################################


def create_batch(receptor, run, i=0):
    """Create batch scripts for jobs that will be sent to cluster, such as
    PIPER docking, models extraction, fpd input preparation and fpd refinement"""
    rec_name = os.path.join(piper_dir, receptor.lower() + '_nmin.pdb')
    lig_name = 'lig.' + "{:04}".format(i) + '_nmin.pdb'
    if run == 'piper':
        with open('run_piper', 'w') as piper:
            piper.write(SBATCH_PIPER.format(decoys=N_ROTS, r=rec_name, l=lig_name))
    elif run == 'decoys':
        with open('run_extract_decoys', 'w') as extract_decoys:
            extract_decoys.write(SBATCH_EXTRACT_TOP_DECOYS % (PIPER_MODELS_EXTRACTED, lig_name))
    elif run == 'prepare_inputs':
        ppk_receptor = os.path.join(prepack_dir, receptor)
        with open('run_prepare_fpd_inputs', 'w') as inputs:
            inputs.write(SBATCH_PREP_FPD_INPUT % (ppk_receptor, refinement_dir, refinement_dir))
    elif run == 'refinement':
        with open(os.path.join(refinement_dir, 'run_refinement'), 'w') as refinement:
            if talaris:
                refinement.write(SBATCH_FPD.format(fpd_version=FPD_TALARIS.format(flags='refine_flags')))
            else:
                refinement.write(SBATCH_FPD.format(fpd_version=FPD.format(flags='refine_flags')))
    elif run == 'clustering':
        with open(os.path.join(clustering_dir, 'run_clustering'), 'w') as cluster:
            if not talaris:
                cluster.write(SBATCH_CLUSTERING.format(native=os.path.join(prepack_dir, 'start.pdb'),
                                                       decoys=os.path.join(refinement_dir, silent_file),
                                                       sc_func='ref2015'))
            else:
                cluster.write(SBATCH_CLUSTERING.format(native=os.path.join(prepack_dir, 'start.pdb'),
                                                       decoys=os.path.join(refinement_dir, silent_file),
                                                       sc_func='talaris14'))


def send_piper_jobs(processed_receptor):
    print("**************Sending PIPER, FlexPepDock and clustering jobs**************")
    receptor_name = os.path.splitext(os.path.basename(receptor_path))[0]
    jobs_list = []
    # PIPER step
    os.chdir(piper_dir)
    ppk_receptor = os.path.splitext(processed_receptor)[0] + '.ppk.pdb'
    for i in range(1, NUM_OF_FRAGS + 1):
        run_dir = "{:02}".format(i)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        os.chdir(run_dir)
        os.system(COPY.format(os.path.join(piper_dir, 'ligands', 'lig.' + "{:04}".format(i) + '_nmin.pdb'),
                              'lig.' + "{:04}".format(i) + '_nmin.pdb'))
        create_batch(receptor_name, 'piper', i)
        create_batch(receptor_name, 'decoys', i)
        create_batch(ppk_receptor, 'prepare_inputs', i)

        # run PIPER docking, extract top 250 decoys and prepare input for refinement
        run_piper = str(subprocess.check_output(RUN_PIPER))
        piper_id = ''.join(d for d in run_piper if d.isdigit())
        RUN_EXTRACT_DECOYS[1] = '--dependency=aftercorr:%s' % piper_id
        run_extract_decoys = str(subprocess.check_output(RUN_EXTRACT_DECOYS))
        extract_decoys_id = ''.join(d for d in run_extract_decoys if d.isdigit())
        RUN_PREP_FPD_INPUTS[1] = '--dependency=aftercorr:%s' % extract_decoys_id
        run_prepare_fpd_inputs = str(subprocess.check_output(RUN_PREP_FPD_INPUTS))
        jobs_list.append(''.join(d for d in run_prepare_fpd_inputs if d.isdigit()))
        os.chdir(piper_dir)

    return jobs_list, receptor_name


def run_piper_fpd_jobs(processed_receptor):
    """Run PIPER, FPD and clustering"""
    jobs_list, receptor_name = send_piper_jobs(processed_receptor)
    # FPD refinement
    if not os.path.exists(refinement_dir):
        os.makedirs(refinement_dir)
    os.chdir(refinement_dir)
    create_batch(receptor_name, 'refinement')
    refine_flags_file()
    all_prep_jobs = ''
    for job_id in jobs_list:
        all_prep_jobs += ':' + job_id  # there are multiple jobs and a semicolon needs to be written before each of them
    RUN_REFINEMENT[1] = '--dependency=aftercorr%s' % all_prep_jobs
    run_refinement = str(subprocess.check_output(RUN_REFINEMENT))

    refinement_id = ''.join(d for d in run_refinement if d.isdigit())

    # Rescoring and Clustering
    if not os.path.exists(clustering_dir):
        os.makedirs(clustering_dir)
    create_batch(receptor_path, 'clustering')
    if not native:  # If native structure was not provided, top scoring structure will be taken as native for rescoring
        with open(os.path.join(refinement_dir, 'extract_model'), 'w') as extract_pdb:
            extract_pdb.write(SBATCH_EXTRACT_TOP_MODEL)
        with open(os.path.join(refinement_dir, 'rescoring'), 'w') as rescore:
            if talaris:
                rescore.write(SBATCH_RESCORING.format(sc_func='talaris14'))
            else:
                rescore.write(SBATCH_RESCORING.format(sc_func='ref2015'))
        RUN_EXTRACT_TOP_MODEL[1] = '--dependency=aftercorr:%s' % refinement_id
        extract_model = str(subprocess.check_output(RUN_EXTRACT_TOP_MODEL))

        extract_model_id = ''.join(d for d in extract_model if d.isdigit())
        RUN_RESCORING[1] = '--dependency=aftercorr:%s' % extract_model_id
        rescoring = str(subprocess.check_output(RUN_RESCORING))

        rescoring_id = ''.join(d for d in rescoring if d.isdigit())
        os.chdir(clustering_dir)
        RUN_CLUSTERING[1] = '--dependency=aftercorr:%s' % rescoring_id
        subprocess.call(RUN_CLUSTERING)
    else:
        os.chdir(clustering_dir)
        RUN_CLUSTERING[1] = '--dependency=aftercorr:%s' % refinement_id
        subprocess.call(RUN_CLUSTERING)
    os.chdir(root)


######################################################################
"""It is not recommended to change the flags, unless you know what you 
 are doing. If you do, flags can be changed in the next 3 functions"""
######################################################################


def fragments_flags_and_cfg(psipred='xxxxx.psipred_ss2', checkpoint='xxxxx.checkpoint'):
    # Create psi_L1.cfg file:
    with open('psi_L1.cfg', 'w') as scores:
        scores.write('#score\tname\tpriority\twght\tmin_allowed\textras\n'
                     'SecondarySimilarity\t350\t2.0\t-\tpsipred\n'
                     'ProfileScoreL1\t200\t1.0\t-\n')
    # Write flags files
    with open('flags', 'w') as flags_file:
        flags_file.write('-in:file:vall\t' + ROSETTA_TOOLS +
                         'fragment_tools/vall.jul19.2011.gz\n'
                         '-in:file:checkpoint\t{check}\n'
                         '-frags:describe_fragments\tfrags.fsc\n'
                         '-frags:frag_sizes\t{len}\n'
                         '-frags:n_candidates\t2000\n'
                         '-frags:n_frags\t100\n'
                         '-out:file:frag_prefix\tfrags\n'
                         '-frags:ss_pred\t{psi} psipred\n'
                         '-frags:scoring:config\tpsi_L1.cfg\n'
                         '-frags:bounded_protocol\ttrue\n'
                         '-mute\tcore.util.prof\n'
                         '-mute\tcore.conformation\n'
                         '-mute\tcore.chemical\n'
                         '-mute\tprotocols.jumping'.format(check=checkpoint, len=pep_length, psi=psipred))


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


def refine_flags_file():
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


def check_native_structure():
    """ Check the length of native struct with the length of the receptor + pep"""
    with open(native_path, 'r') as n:
        native_calphas = 0
        native_sequence = ''
        for line in n:
            if line[13:15] == 'CA':
                native_calphas += 1
#               native_sequence += line[17:20]  # TODO: you are not comparing the sequence !!!
    with open(receptor_path, 'r') as rec:
        rec_calphas = 0
        for line in rec:
            if line[13:15] == 'CA':
                rec_calphas += 1
        complex_calphas = rec_calphas + pep_length
    if native_calphas != complex_calphas:
        print(BAD_NATIVE)
        sys.exit()
    else:
        return

    # with open(native_path) as nat:
    #     residues = ''
    #     cur_line = nat.readline()
    #     for i in range(receptor_len + pep_length):
    #         if cur_line[:4] == 'ATOM':
    #             if cur_line[17:20] == THREE_TO_ONE_AA[complex_seq[i]]:
    #                 residues += complex_seq[i]
    #                 cur_resi = cur_line[22:27].strip()
    #                 while cur_resi == cur_line[22:27].strip():
    #                     cur_line = nat.readline()
    #             else:
    #                 print(BAD_NATIVE)
    #                 sys.exit()
    #     if residues == complex_seq:
    #         return
    # print(BAD_NATIVE)
    # sys.exit()


def cut_file(col_num):
    if col_num == 1:
        file_name = 'psipred_ss2_pep'
        original_file = 'xxxxx.psipred_ss2'
    else:
        file_name = 'checkpoint'
        original_file = 'xxxxx.checkpoint'
    with open(original_file, 'r') as f:
        f_lines = f.readlines()
        i = 0
        for line_num in range(1, len(f_lines) + 1):
            if f_lines[line_num].split()[col_num] == peptide_seq[i]:
                i += 1
                if i == len(peptide_seq):
                    with open(file_name, 'w') as cut_f:
                        if file_name == 'psipred_ss2_pep':
                            cut_f.write(f_lines[0])
                            for j in range(pep_length - 1, -1, -1):
                                new_line = f_lines[line_num - j].split()
                                new_line[0] = str(i - j)
                                cut_f.write('\t'.join(new_line) + '\n')
                        else:
                            cut_f.write(str(pep_length) + '\n')
                            for j in range(pep_length - 1, -1, -1):
                                cut_f.write(f_lines[line_num - j])
                    break
            else:
                i = 0
    return file_name


def create_psipred_from_full_protein(full_pep_fasta):
    """Reads fasta and calls make_pick fragments"""
    if not os.path.exists(frag_picker_dir):
        os.makedirs(frag_picker_dir)
    os.chdir(frag_picker_dir)
    with open(full_pep_fasta) as fasta:
        full_seq = fasta.readlines()
        if full_seq[0][0] == '>':
            full_seq = full_seq[1:]
            full_seq = "".join(line.strip() for line in full_seq)
        else:
            full_seq = "".join(line.strip() for line in full_seq)

    make_pick_fragments(full_seq)


def make_pick_fragments(pep_seq):
    """Run fragment picker"""
    if not os.path.exists(frag_picker_dir):
        os.makedirs(frag_picker_dir)
    # Create fasta file:
    with open(os.path.join(frag_picker_dir, 'xxxxx.fasta'), 'w') as fasta_file:
        fasta_file.write('>|' + pep_seq + '\n' + pep_seq + '\n')
    os.chdir(frag_picker_dir)

    os.system(MAKE_FRAGMENTS.format('xxxxx.fasta'))  # Run make_fragments.pl script
    if sec_struct:
        os.rename('xxxxx.psipred_ss2', 'xxxxx.psipred_ss2_orig')
        with open('xxxxx.psipred_ss2', 'w') as psi_new:
            with open('xxxxx.psipred_ss2_orig', 'r') as psipred:
                psipred_lines = psipred.readlines()
                for i, line in enumerate(psipred_lines):
                    new_line = line.split()
                    new_line[2] = ss_pred[i]
                    if ss_pred[i] == 'C':
                        new_line[3] = '0.700'
                        new_line[4] = '0.290'
                        new_line[5] = '0.010'
                    elif ss_pred[i] == 'H':
                        new_line[3] = '0.290'
                        new_line[4] = '0.700'
                        new_line[5] = '0.010'
                    elif ss_pred[i] == 'E':
                        new_line[3] = '0.290'
                        new_line[4] = '0.010'
                        new_line[5] = '0.700'
                    psi_new.write('\t'.join(new_line) + '\n')
    elif full_p:
        psi_p_file = cut_file(1)  # cut psipred_ss2
        checkpoint_file = cut_file(0)  # cut checkpoint
        fragments_flags_and_cfg(psi_p_file, checkpoint_file)  # Write flags files
    else:
        fragments_flags_and_cfg()
    print("**************Picking fragments**************")
    os.system(FRAG_PICKER)  # Run fragment picker
    os.system(COPY.format(FRAGS_FILE.format(pep_length), root))  # Copy fragments file (frags.100.nmers)
    os.chdir(root)


def create_params_file(frags):
    """Read only needed values from frags_file and store them in frags_parameters file"""
    if not os.path.isfile('frags_parameters'):
        parameters_sets = []
        frags_parameters = open('frags_parameters', 'w+')
        with open(frags) as frags_file:
            all_file_lines = frags_file.readlines()
        # Get parameters and save them in parameters_sets list
        j = 0
        for i in range(2, len(all_file_lines)):
            if j % (pep_length + 1) == 0:
                first_line_in_set = all_file_lines[i]
                last_line_in_set = all_file_lines[i+(pep_length-1)]
                seq = ""
                for k in range(i, i+pep_length):
                    line = all_file_lines[k].split()
                    seq += line[3]
                first_line_words = first_line_in_set.split()
                last_line_words = last_line_in_set.split()
                pdb = first_line_words[0]
                chain = first_line_words[1]
                start_res = first_line_words[2]
                end_res = last_line_words[2]
                parameters_sets.append([pdb, chain, start_res, end_res, seq])
            j += 1
        for item in parameters_sets:
            for par in item:
                frags_parameters.write("%s\t" % par)
            frags_parameters.write("\n")
        frags_parameters.close()
    with open('frags_parameters', 'r') as params:
        all_frags = params.readlines()
    return all_frags


def count_pdbs(folder):
    """Count files in a directory"""
    return len([frag for frag in os.listdir(folder) if
                os.path.isfile(os.path.join(folder, frag)) and
                os.path.splitext(os.path.basename(frag))[1] == '.pdb'])


def extract_frag(pdb, start, end, outfile):
    """Extract fragment from a full chain"""
    with open(outfile, 'w') as frag:
        with open(pdb, 'r') as full_pdb:
            cur_line = full_pdb.readline()
            while cur_line[22:27].strip() != start:
                cur_line = full_pdb.readline()
                if not cur_line:
                    return
            while cur_line[22:27].strip() != str(int(end) + 1):
                frag.write(cur_line)
                cur_line = full_pdb.readline()
                if not cur_line:
                    return


def bad_frag(fragment):
    """Move bad fragment to separate directory"""
    print("Bad fragment. it will be saved in separate directory 'bad_fragments'")
    if not os.path.exists(bad_frags_dir):
        os.makedirs(bad_frags_dir)
    os.rename(fragment, os.path.join(bad_frags_dir, fragment))
    return False


def review_frag(outfile, sequence):
    """Check whether fragment is of a right length, right sequence and
    does not contain zero occupancy atoms. Otherwise move to 'bad_fragments' directory."""
    if os.path.getsize(outfile) <= 0:
        bad_frag(outfile)
    with open(outfile) as frag:
        residues = ''
        cur_line = frag.readline()
        for i in range(pep_length):
            if cur_line[:4] == 'ATOM' and cur_line[17:20] == THREE_TO_ONE_AA[sequence[i]]:
                if cur_line[54:60].strip() == 0.00:
                    print("Zero occupancy atoms!")
                    return bad_frag(outfile)
                residues += sequence[i]
                cur_resi = cur_line[22:27].strip()
                while cur_resi == cur_line[22:27].strip():
                    cur_line = frag.readline()
            else:
                return bad_frag(outfile)
        if residues == sequence:
            return True
    return bad_frag(outfile)


# TODO: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def renumber_frag(fragment):
    renumbered = []
    with open(fragment, 'r') as frag:
        pdb_lines = frag.readlines()
        i = 1  # residues counter

    for j in range(len(pdb_lines)):
        if pdb_lines[j][0:4] != 'ATOM' and pdb_lines[j][0:6] != 'HETATM':
            continue
        elif pdb_lines[j][17:20] != pdb_lines[j + 1][17:20]:
            i += 1
        elif i < 10:
            new_line = list(pdb_lines[j])
            new_line[22:26] = str(i) + '   '
            renumbered.append(''.join(new_line))
        else:
            new_line = list(pdb_lines[j])
            new_line[22:26] = str(i) + '  '
            renumbered.append(''.join(new_line))
    os.remove(fragment)
    with open(fragment, 'w') as new_structure:
        for renumbered_line in renumbered:
            new_structure.write(renumbered_line)


def process_frags(pep_sequence, fragments, add_frags_num=0):
    """Process and extract frags, filter bad fragments"""
    pdb_resfiles_dict = dict()  # pdbs and their resfiles
    # # in case there are fragments
    # if os.path.isdir(fragments_dir) and os.path.isdir(resfiles_dir):
    #     if count_pdbs(fragments_dir) == 50 and len([resf for resf in os.listdir(resfiles_dir) if
    #                                                os.path.isfile(os.path.join(resfiles_dir, resf))]):
    #         for frag in os.listdir(fragments_dir):
    #             if os.path.splitext(os.path.basename(frag))[1] == '.pdb':
    #                 pdb_resfiles_dict[os.path.basename(frag)] = 'resfile_'+os.path.splitext(os.path.basename(frag))[0]
    #         return pdb_resfiles_dict
    #     else:
    #         os.removedirs(fragments_dir)
    #         os.removedirs(resfiles_dir)

    # # create directory for top 50 frags
    # if not os.path.exists(fragments_dir):
    #     os.makedirs(fragments_dir)
    #
    # # create directory for storing resfiles
    # if not os.path.exists(resfiles_dir):
    #     os.makedirs(resfiles_dir)
    #
    # pdb_resfiles_dict = dict()  # pdbs and their resfiles
    # os.chdir(fragments_dir)

    # extract and append parameters to different lists
    pdbs = []
    chains = []
    start_res = []
    end_res = []
    sequences = []
    # Start and end residues numbers will not be used for fragment extraction - the fragments will be extracted
    # from fasta of sequentially renumbered pdb. Both, fasta file and renumbered pdb, are output of clean_pdb.py
    # However, the original numbers will be used in fragment and it's resfile names
    for frag in fragments:
        pdbs.append(frag.split()[0])
        chains.append(frag.split()[1])
        start_res.append(frag.split()[2])
        end_res.append(frag.split()[3])
        sequences.append(frag.split()[4])

    # create directory for top 50 frags
    if not os.path.exists(fragments_dir):
        os.makedirs(fragments_dir)
    # create directory for storing resfiles
    if not os.path.exists(resfiles_dir):
        os.makedirs(resfiles_dir)

    os.chdir(fragments_dir)
    print("**************Extracting fragments**************")
    # Fetch PDBs, extract fragments and create resfiles
    for pdb, chain, start, end, sequence in zip(pdbs, chains, start_res, end_res, sequences):
        fragment_name = pdb + '.' + chain + '.' + start + '.' + end
        outfile = fragment_name + '.pdb'
        if chain == '_':
            chain = 'A'
        os.system(CLEAN_PDB.format(pdb, chain))  # get clean pdb and it's fasta

        pdb_full = PDB.format(pdb.upper(), chain)  # names of clean_pdb output files
        fasta_name = FASTA.format(pdb.upper(), chain)

        if os.path.exists(pdb_full):
            print("Extracting fragment")

            with open(fasta_name, 'r') as f:
                fasta = f.read()
            clean_fasta = fasta[fasta.find('\n') + 1:]
            # These fasta_start and fasta_end numbers are only temporary numbers for fragments extraction
            fasta_start = clean_fasta.find(sequence) + 1  # +1 because of zero-based numbering
            if fasta_start == 0:  # -1 would mean 'sequence doesn't exist', but we added 1
                print("no matching sequence")
                continue
            fasta_end = fasta_start + pep_length - 1

            extract_frag(pdb_full, str(fasta_start), str(fasta_end), outfile)

            is_frag_ok = review_frag(outfile, sequence)
            if is_frag_ok:
                os.remove(pdb_full)
                os.remove(fasta_name)
                # renumber_frag(outfile)      # todo: uncomment after it is finished
                frags_count = count_pdbs(fragments_dir)
                print("creating resfile")
                create_resfile(pep_sequence, chain, fasta_start, sequence, fragment_name)
                pdb_resfiles_dict[outfile] = 'resfile_' + fragment_name
                if frags_count >= NUM_OF_FRAGS + add_frags_num:
                    print("**************Finished with fragments**************")
                    break
            else:
                print("Failed to extract fragment")  # or the fragment went to bad_fragments
                os.remove(pdb_full)
                os.remove(fasta_name)
                continue
        else:
            print("Failed to fetch pdb")  # The PDB can be obsolete
            continue  # if failed to fetch PDB
    os.chdir(root)
    return pdb_resfiles_dict


def create_resfile(ori_seq, chain, start, sequence, fragment_name):
    """Create resfiles for each fragment individually"""
    cur_resfile = os.path.join(resfiles_dir, 'resfile_%s')
    # Create resfile for each fragment
    resfile = open(cur_resfile % fragment_name, 'w')
    resfile.write('NATRO\nstart')
    if chain == '_':
        chain = 'A'
    for i, res in enumerate(sequence):
        if res == ori_seq[i]:
            resfile.write('\n' + str(i + int(start)) + ' ' + chain + ' NATRO')
        else:
            resfile.write('\n' + str(i + int(start)) + ' ' + chain + ' PIKAA ' + ori_seq[i] +
                          ' EX 1 EX 2')
    resfile.close()


def create_xml(pdb_resfile_dict):
    """Create xml file for jd3_fixbb with 50 jobs with different pdbs and refiles"""
    job_string = '<Job>\n' \
                 '\t<Input>\n' \
                 '\t\t<PDB filename="../{}"/>\n' \
                 '\t</Input>\n' \
                 '\t<TASKOPERATIONS>\n' \
                 '\t\t<ReadResfile name="read_resfile" filename="../../resfiles/{}"/>\n' \
                 '\t</TASKOPERATIONS>\n' \
                 '\t<PackRotamersMover name="mover" scorefxn="{}" task_operations="read_resfile"/>\n' \
                 '</Job>\n'
    if not os.path.exists(fixbb_dir):
        os.makedirs(fixbb_dir)
    with open(fixbb_dir + '/design.xml', 'w') as xml_file:
        xml_file.write('<JobDefinitionFile>\n')
        if talaris:
            xml_file.write('<Common>\n'                          
                           '\t<SCOREFXNS>\n'
                           '\t\t<ScoreFunction name="Talaris14" weights="talaris2014.wts"/>\n'
                           '\t</SCOREFXNS>\n'
                           '</Common>\n')
            for pdb, resfile in pdb_resfile_dict.items():
                xml_file.write(job_string.format(pdb, resfile, 'Talaris14'))
        else:
            xml_file.write('<Common>\n'
                           '\t<SCOREFXNS>\n'
                           '\t\t<ScoreFunction name="ref2015" weights="ref2015.wts"/>\n'
                           '\t</SCOREFXNS>\n'
                           '</Common>\n')
            for pdb, resfile in pdb_resfile_dict.items():
                xml_file.write(job_string.format(pdb, resfile, 'ref2015'))
        xml_file.write('</JobDefinitionFile>')


def check_designed_frags():
    """Review fragments after design"""
    for frag in os.listdir('.'):
        if os.path.splitext(os.path.basename(frag))[1] == '.pdb':
            review_frag(frag, peptide_seq)
    frags_count = count_pdbs(fixbb_dir)
    if frags_count < NUM_OF_FRAGS:
        return NUM_OF_FRAGS - frags_count
    else:
        return False


def extract_more_frags(n_frags, defective_from_fixbb):
    """If there were wrong fragments after fixbb design"""
    os.chdir(fragments_dir)
    with open(os.path.join(root, 'frags_parameters'), 'r') as param_file:
        add_frags = param_file.readlines()
    if os.path.isdir(bad_frags_dir):
        bad_frags_shift = count_pdbs(bad_frags_dir) + defective_from_fixbb
    else:
        bad_frags_shift = defective_from_fixbb
    add_frags = add_frags[NUM_OF_FRAGS + bad_frags_shift:]
    return process_frags(peptide_seq, add_frags, n_frags)


def run_fixbb(pdb_and_resfiles):
    """Run fixbb design (option to restore talaris behaviour)"""
    if not os.path.exists(fixbb_dir):
        os.makedirs(fixbb_dir)
    os.chdir(fixbb_dir)
    print("**************Fixbb design...**************")
    if talaris:
        if jd3:
            create_xml(pdb_and_resfiles)
            os.system(FIXBB_JD3_TALARIS.format('design.xml'))
        else:
            for pdb, resfile in pdb_and_resfiles.items():
                os.system(FIXBB_TALARIS.format(frag=os.path.join(fragments_dir, pdb),
                                               resfile=os.path.join(resfiles_dir, resfile)))
    else:
        if jd3:
            create_xml(pdb_and_resfiles)
            os.system(FIXBB_JD3.format('design.xml'))
        else:
            for pdb, resfile in pdb_and_resfiles.items():
                os.system(FIXBB.format(frag=os.path.join(fragments_dir, pdb),
                                       resfile=os.path.join(resfiles_dir, resfile)))

    print("Done!")
    # If we need to extract additional fragments for more then once, we need to add bad fragments
    # to frags_file shift also (but only starting from second time)
    if os.path.isdir(bad_frags_dir):
        already_defective = count_pdbs(bad_frags_dir)
    else:
        already_defective = 0
    # check frags and return False if all of them are OK, or number of defective fragments
    fragments_needed = check_designed_frags()
    if not fragments_needed:
        os.chdir(root)
        return
    else:
        print("**************Some fragments were defective. Extracting more fragments**************")
        # extract more frags and create a new dictionary for creating an xml
        new_frag_resfile_dict = extract_more_frags(fragments_needed, already_defective)
        create_xml(new_frag_resfile_dict)
        run_fixbb(new_frag_resfile_dict)


def rename_chain(structure, chain_id):
    """Rename chain id of a structure"""
    renamed_struct = []
    with open(structure, 'r') as pdb:
        pdb_lines = pdb.readlines()
    for line in pdb_lines:
        if line[0:4] != 'ATOM' and line[0:6] != 'HETATM':
            continue
        else:
            new_line = list(line)
            new_line[21] = chain_id
            renamed_struct.append(''.join(new_line))
    os.remove(structure)
    with open(structure, 'w') as new_structure:
        for new_chain_line in renamed_struct:
            new_structure.write(new_chain_line)


def process_for_piper(receptor):
    """Prepare inputs for piper run"""
    frags_list = []
    if not os.path.exists(piper_dir):  # Create directory
        os.makedirs(piper_dir)
    print("**************Preparing PIPER inputs**************")
    for frag in os.listdir(fixbb_dir):
        if os.path.splitext(frag)[1] == '.pdb':   # Rename chain ID to 'B'
            rename_chain(os.path.join(fixbb_dir, frag), 'B')
            frags_list.append(os.path.basename(frag))
    with open(os.path.join(piper_dir, 'runs_list'), 'w') as runs_list:  # Create list 1 - 50
        for i in range(1, NUM_OF_FRAGS + 1):
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
        os.system(PDBPREP.format(lig))
        os.system(PDBNMD.format(lig))

    # prepare receptor for piper
    os.system(COPY.format(receptor, piper_dir))
    os.chdir(piper_dir)

    os.system(CLEAN_PDB.format(receptor, 'nochain'))  # Clean the receptor
    rename_chain(os.path.basename(receptor)[:-4] + '_nochain.pdb', 'A')
    name_for_piper = os.path.basename(receptor).lower()
    os.rename(os.path.splitext(os.path.basename(receptor))[0] + '_nochain.pdb',
              name_for_piper)
    os.system(PDBPREP.format(name_for_piper))
    os.system(PDBNMD.format(name_for_piper))
    os.chdir(root)


def build_peptide(pep):
    """Create directory for prepacking and, extended peptide and change its chain id to 'B'"""
    print("Building extended peptide for prepacking")
    if not os.path.exists(prepack_dir):
        os.makedirs(prepack_dir)
    os.chdir(prepack_dir)
    # Build extended peptide
    os.system(BUILD_PEPTIDE.format(pep))

    # Change chain ID to 'B'
    rename_chain('peptide.pdb', 'B')
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
        os.system(PREPACK_TALARIS)
    else:
        os.system(PREPACK)
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
    parser.add_argument('--jd3', dest='job_distributor', action='store_true', default=False)
    parser.add_argument('--sec_struct', dest='ss_pred', default=None)
    parser.add_argument('--pep_from_fasta', dest='full_prot_seq', default=None)  # pdb_name,chain (e.g. 1abc,A)

    return parser


def run_protocol(peptide_sequence, receptor):

    if native:
        check_native_structure()

    if full_p:
        create_psipred_from_full_protein(os.path.abspath(full_p))
    else:
        make_pick_fragments(peptide_sequence)

    all_frags = create_params_file(FRAGS_FILE.format(str(pep_length)))

    # extract fragments, create resfiles and return a dictionary of fragments names and matching resfiles names
    pdb_and_resfiles = process_frags(peptide_seq, all_frags)

    # run fixbb
    run_fixbb(pdb_and_resfiles)

    # process ligands and receptor for piper run
    process_for_piper(receptor)  # only basename

    build_peptide(os.path.abspath(sys.argv[2]))  # build extended peptide and rename it's chain id to 'B'

    prepack_receptor(os.path.basename(receptor).lower())

    # run piper docking, extract top 250 models from each run, run FlexPepDock, clustering,
    # rescoring and get top 10 models

    run_piper_fpd_jobs(os.path.basename(receptor).lower())


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
    jd3 = arguments.job_distributor
    sec_struct = arguments.ss_pred
    full_p = arguments.full_prot_seq

    # Define all the directories that will be created:
    root = os.getcwd()
    frag_picker_dir = os.path.join(root, 'frag_picker')
    fragments_dir = os.path.join(root, 'top_50_frags')
    resfiles_dir = os.path.join(root, 'resfiles')
    fixbb_dir = os.path.join(fragments_dir, 'fixbb')
    piper_dir = os.path.join(root, 'piper')
    ligands_dir = os.path.join(piper_dir, 'ligands')
    prepack_dir = os.path.join(root, 'prepacking')
    refinement_dir = os.path.join(root, 'refinement')
    clustering_dir = os.path.join(refinement_dir, 'clustering')
    bad_frags_dir = 'bad_fragments'

    final_dir = 'FINAL_RESULTS'  # top 10 models and score file

    # TODO: !!!!!!!!!!!!!!!!! I need to pull this up from the main
    if native:
        native_path = os.path.abspath(native)
        silent_file = 'decoys.silent'
    else:
        silent_file = 'decoys_rescored.silent'

    if sec_struct:
        with open(sec_struct, 'r') as ss_h:
            ss_pred = ss_h.readline().strip()
        for char in ss_pred:
            if char not in PSIPRED_OUTPUT:
                print('Wrong secondary structure file format. A valid file should contain only 1 line with C, H or E '
                      'letters')
                sys.exit()

    run_protocol(peptide_seq, receptor_path)
