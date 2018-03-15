#!/usr/bin/python

import sys
import os
import time

#########################
""" Change paths here """
#########################

ROSETTA_DIR = '/vol/ek/Home/alisa/rosetta/Rosetta/'
ROSETTA_2016_DIR = '/vol/ek/share/rosetta/rosetta_src_2016.20.58704_bundle/'

ROSETTA_DB = ROSETTA_DIR + 'main/database'
ROSETTA_2016_DB = ROSETTA_2016_DIR + 'main/database'

ROSETTA_BIN = ROSETTA_DIR + 'main/source/bin/'
ROSETTA_2016_BIN = ROSETTA_2016_DIR + 'main/source/bin/'

ROSETTA_TOOLS = ROSETTA_DIR + 'tools/protein_tools/scripts/'

MAKE_FRAGMENTS_DIR = '/vol/ek/share/scripts/global_pep_dock/fragpicker_setup/'

PIPER_DIR = '/vol/ek/Home/alisa/PIPER/'
PIPER_BIN = PIPER_DIR + 'bin/'


##################################
"""DO NOT change these commands"""
##################################

# Commands (ROSETTA)

GET_PDB = ROSETTA_TOOLS + 'clean_pdb.py {} {}'

FIXBB_JD3 = 'mpirun -n 6 ' + ROSETTA_BIN + 'fixbb_jd3.mpiserialization.linuxgccrelease' \
            ' -database ' + ROSETTA_DB + ' -in:file:job_definition_file {} > fixbb.log'
FIXBB_JD3_TALARIS = 'mpirun -n 6 ' + ROSETTA_BIN + 'fixbb_jd3.mpiserialization.linuxgccrelease' \
                     ' -database ' + ROSETTA_DB + ' -restore_talaris_behavior' \
                     ' -in:file:job_definition_file {} > fixbb.log'

BUILD_PEPTIDE = ROSETTA_BIN + 'BuildPeptide.linuxgccrelease -in:file:fasta {}' \
                ' -database ' + ROSETTA_DB + ' -out:file:o peptide.pdb ' \
                '> build_peptide.log'

MAKE_FRAGMENTS = 'perl ' + MAKE_FRAGMENTS_DIR + 'make_fragments.pl' \
                 ' -verbose -id xxxxx {} 2>log'

FRAG_PICKER = ROSETTA_BIN + 'fragment_picker.linuxgccrelease' \
              ' -database ' + ROSETTA_DB + ' @flags >makeFrags.log'

PREPACK = ROSETTA_BIN + 'FlexPepDocking.default.linuxgccrelease -database ' + \
          ROSETTA_DB + ' @prepack_flags >ppk.log'
PREPACK_TALARIS = ROSETTA_2016_BIN + 'FlexPepDocking.mpi.linuxgccrelease -database ' + \
                  ROSETTA_DB + ' @prepack_flags >ppk.log'

FPD_REFINEMENT = 'mpirun ' + ROSETTA_BIN + 'FlexPepDocking.mpi.linuxgccrelease -database' + ROSETTA_DB + \
                 '@refine_flags >refinement_log'
FPD_REFINEMENT_TALARIS = 'mpirun ' + ROSETTA_2016_BIN + 'FlexPepDocking.mpi.linuxgccrelease -database' +\
                         ROSETTA_2016_DB + ' @refine_flags >refinement_log'

# Commands (PIPER)

PDBPREP = 'perl ' + PIPER_BIN + 'phplibbin/pdbprep.pl {}'
PDBNMD = 'perl ' + PIPER_BIN + 'phplibbin/pdbnmd.pl "{}"' + ' ?'

PIPER_DOCKING = PIPER_BIN + 'piper.acpharis.omp.20120803 -vv -c1.0 -k4 --msur_k=1.0 --maskr=1.0 ' \
                '-T FFTW_EXHAUSTIVE -R 70000 -t 1 -p ' + PIPER_DIR + 'prms/atoms.0.0.4.prm.ms.3cap+0.5ace.Hr0rec -f ' \
                + PIPER_DIR + 'prms/coeffs.0.0.4.motif -r ' + PIPER_DIR + 'prms/rot70k.0.0.4.prm {r} {l} >piper.log'

#####################################################################################
"""WE are using SLURM workload manager. If you are not - change the commands below"""
#####################################################################################

RUN_PIPER_SLURM = 'sbatch --mem=1500m --nice=10000 run_piper'

# This will be created under the name run_piper. PIPER_DOCKING is the only unchangeable command here
SBATCH_PIPER = '#!/bin/sh\n' \
               '#SBATCH --ntasks=1\n' \
               '#SBATCH --time=20:00:00\n' \
               '#SBATCH --get-user-env\n'\
               + PIPER_DOCKING
APPLY_FTRESULTS = 'python ' + PIPER_DIR + 'apply_ftresult.py -i {model} ft.000.00 ' \
                  + PIPER_DIR + 'prms/rot70k.0.0.4.prm {lig} --out-prefix {out}'


# SBATCH_EXTRACT_TOP_DECOYS = "#!/bin/sh\n#SBATCH --nodes 1\n#SBATCH --ntasks=1\n#SBATCH --time=10:00:00\n" \
#                             "#SBATCH --get-user-env\nfor f in `awk '{print $1}' ft.000.00 | head -250`;" \
#                             "do if [ ! -f {f}.pdb ]; then " + PIPER_DIR + "apply_ftresult.py -i $f ft.000.00 "\
#                             + PIPER_DIR + "prms/rot70k.0.0.4.prm %s --out-prefix $f;fi;done"

# RUN_PIPER = "job_id_PIPER_dock=`sbatch --mem=1500m --nice=10000 run_piper| awk '{print $NF}' | grep [0-9]` && " \
#             "job_id_extraction=$(sbatch --dependency=afterany:$job_id_PIPER_dock --mem-per-cpu=1500m " \
#             "run_extract_decoys)"  # This one will run SBATCH_PIPER and SBATCH_EXTRACT TOP_DECOYS

RUN_FPD_SLURM = 'sbatch --mem-per-cpu=1600m run_refinement'
SBATCH_FPD = '#!/bin/bash\n' \
             '#SBATCH --ntasks=300\n' \
             '#SBATCH --time=50:00:00\n' \
             '#SBATCH --get-user-env\n' \


#######################################################################################

# Other commands

COPY = 'cp {} {}'

# Other constants

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


######################################################################
"""It is not recommended to change the flags, unless you know what you 
 are doing. If you do, flags can be changed in the next 2 functions"""
######################################################################


def make_pick_fragments(pep_seq):
    """Run fragment picker"""

    if not os.path.exists(frag_picker_dir):
        os.makedirs(frag_picker_dir)
    # Create fasta file:
    with open(os.path.join(frag_picker_dir, 'xxxxx.fasta'), 'w') as fasta_file:
        fasta_file.write('>|' + pep_seq + '\n' + pep_seq + '\n')
    os.chdir(frag_picker_dir)
    os.system(MAKE_FRAGMENTS.format('xxxxx.fasta'))  # Run make_fragments.pl script
    # Create psi_L1.cfg file:
    with open('psi_L1.cfg', 'w') as scores:
        scores.write('#score\tname\tpriority\twght\tmin_allowed\textras\n'
                     'SecondarySimilarity\t350\t2.0\t-\tpsipred\n'
                     'ProfileScoreL1\t200\t1.0\t-\n')
    # Write flags
    with open('flags', 'w') as flags_file:
        flags_file.write('-in:file:vall\t' + ROSETTA_DB + '/sampling/'
                         'filtered.vall.dat.2006-05-05.gz\n'
                         '-in:file:checkpoint\txxxxx.checkpoint\n'
                         '-frags:describe_fragments\tfrags.fsc\n'
                         '-frags:frag_sizes\t' + str(pep_length) + '\n'
                         '-frags:n_candidates\t2000\n'
                         '-frags:n_frags\t100\n'
                         '-out:file:frag_prefix\tfrags\n'
                         '-frags:ss_pred\txxxxx.psipred_ss2 psipred\n'
                         '-frags:scoring:config\tpsi_L1.cfg\n'
                         '-frags:bounded_protocol\ttrue\n'
                         '-mute\tcore.util.prof\n'
                         '-mute\tcore.conformation\n'
                         '-mute\tcore.chemical\n'
                         '-mute\tprotocols.jumping')
    os.system(FRAG_PICKER)
    os.system(COPY.format(FRAGS_FILE.format(pep_length), root))
    os.chdir(root)


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
                    '-mute protocols.jd2.PDBJobInputter\n')


def refine_flags_file():
    with open('refine_flags', 'w') as flags:
        flags.write('-in:file:l input_list\n'
                    '-scorefile score.sc\n'
                    '-out:pdb_gz\n'
                    '-out:file:silent_struct_type binary\n'
                    '-out:file:silent decoys.silent\n'
                    '-min_receptor_bb\n'
                    '-lowres_preoptimize\n'
                    '-flexPepDocking:pep_refine\n'
                    '-flexPepDocking:flexpep_score_only\n'
                    '-ex1\n'
                    '-ex2aro\n'
                    '-use_input_sc\n'
                    '-unboundrot receptor.pdb\n'
                    '-mute protocols.moves.RigidBodyMover\n'
                    '-mute core.chemical\n'
                    '-mute core.scoring.etable\n'
                    '-mute protocols.evalution\n'
                    '-mute core.pack.rotamer_trials\n'
                    '-mute protocols.abinitio.FragmentMover\n'
                    '-mute core.fragment\n'
                    '-mute protocols.jd2.PDBJobInputter')


######################################
"""End of flags creating functions"""
######################################

################################################
"""ATTENTION! Do not change the code below!!!"""
################################################


def create_params_file(frags):
    """Read only needed values from frags_file and store them in frags_parameters file"""

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


def bad_frag(fragment):  # remove bad fragments
    print("Bad fragment. it will be saved in separate directory 'bad_fragments'")
    if not os.path.exists('bad_fragments'):
        os.makedirs('bad_fragments')
    os.rename(fragment, 'bad_fragments/' + fragment)
    return False


def review_frags(outfile, start, end):
    """Check if fragment length is correct and there are no zero occupancy atoms"""

    with open(outfile) as frag:
        residues = set()
        cur_line = frag.readline()
        while 'ATOM' not in cur_line[0:4]:  # find ATOM lines
            cur_line = frag.readline()

            # if there there no atoms...
            if not cur_line:
                return bad_frag(outfile)

        cur_line = frag.readline()  # line with a first atom

        while cur_line[22:27].strip() != start:
            cur_line = frag.readline()

            # if there is no start residue
            if not cur_line:
                return bad_frag(outfile)

        cur_line = frag.readline()  # first atom of the start residue
        while cur_line[22:27].strip() != end:

            # check occupancy
            if cur_line[54:60].strip() == 0.00:
                print("Zero occupancy atoms")
                return bad_frag(outfile)

            residues.add(cur_line[22:27])
            cur_line = frag.readline()

            # if there is no end residue
            if not cur_line:
                return bad_frag(outfile)

        residues.add(cur_line[22:27])
    if len(residues) != (int(end) - int(start) + 1):
        bad_frag(outfile)
        return False

    else:
        return True


def review_fasta_frag(outfile, sequence):
    """Review fragments that have different numbering and were extracted with fasta"""

    with open(outfile) as frag:
        residues = ''
        cur_line = frag.readline()
        for i in range(pep_length):
            if cur_line[:4] == 'ATOM' and cur_line[17:20] == THREE_TO_ONE_AA[sequence[i]]:
                residues += sequence[i]
                cur_resi = cur_line[22:27].strip()
                while cur_resi == cur_line[22:27].strip():
                    cur_line = frag.readline()
            else:
                return bad_frag(outfile)
        if residues == sequence:
            return True
        return False


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


def process_frags(pep_sequence):

    # Open the frags_parameters, extract and append parameters to different lists
    with open('frags_parameters', 'r') as f:
        fragments = f.readlines()
    pdbs = []
    chains = []
    start_res = []
    end_res = []
    sequences = []
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

    pdb_resfiles_dict = dict()  # names of pdbs and their resfiles

    os.chdir('top_50_frags')

    # Fetch PDBs, extract fragments and create resfiles
    for pdb, chain, start, end, sequence in zip(pdbs, chains, start_res, end_res, sequences):
        fragment_name = pdb + '.' + chain + '.' + start + '.' + end
        outfile = fragment_name + '.pdb'
        if chain == '_':
            chain = 'A'
        os.system(GET_PDB.format(pdb, chain))  # get and clean pdb and fasta

        pdb_full = '{}_{}.pdb'.format(pdb.upper(), chain)
        fasta_name = '{}_{}.fasta'.format(pdb.upper(), chain)

        if os.path.exists(pdb_full):

            print("Extracting fragment")

            extract_frag(pdb_full, start, end, outfile)
            is_frag_ok = False

            if not os.path.exists(outfile) or os.path.getsize(outfile) == 0:
                print("Trying to extract from FASTA")

                with open(fasta_name, 'r') as f:
                    fasta = f.read()
                clean_fasta = fasta[fasta.find('\n') + 1:]
                fasta_start = clean_fasta.find(sequence) + 1  # +1 because of zero-based numbering
                if fasta_start == 0:  # -1 would mean 'sequence doesn't exist', but we added 1
                    print("no matching sequence")
                    continue
                fasta_end = fasta_start + pep_length - 1

                extract_frag(pdb_full, str(fasta_start), str(fasta_end), outfile)

                is_frag_ok = review_fasta_frag(outfile, sequence)
                if is_frag_ok:
                    print("success!")
                else:
                    print("Failed to extract fragment")
                    os.remove(pdb_full)
                    os.remove(fasta_name)
                    continue
            else:
                is_frag_ok = review_frags(outfile, start, end)

            os.remove(pdb_full)
            os.remove(fasta_name)

            if is_frag_ok:
                cur_dir = os.getcwd()
                frags_count = len([frag for frag in os.listdir('.') if
                                   os.path.isfile(os.path.join(cur_dir, frag))])
                print("creating resfile")
                create_resfile(pep_sequence, chain, start, sequence, fragment_name)
                pdb_resfiles_dict[outfile] = 'resfile_' + fragment_name
                if frags_count >= 50:
                    print("Got top 50 fragments")
                    break
            else:
                continue
        else:
            print("Failed to fetch pdb")  # The PDB could be obsolete
            continue  # if failed to fetch PDB

    os.chdir(root)
    return pdb_resfiles_dict


def create_resfile(ori_seq, chain, start, sequence, fragment_name):

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
    job_string = '<Job>\n' \
                 '\t<Input>\n' \
                 '\t\t<PDB filename="../{}"/>\n' \
                 '\t</Input>\n' \
                 '\t<TASKOPERATIONS>\n' \
                 '\t\t<ReadResfile name="read_resfile" filename="../../resfiles/{}"/>\n' \
                 '\t</TASKOPERATIONS>\n' \
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
            xml_file.write(job_string.format(pdb, resfile))
        xml_file.write('</JobDefinitionFile>')


def run_fixbb():
    os.chdir(fixbb_dir)
    print("running fixbb design")
    if talaris:
        os.system(FIXBB_JD3_TALARIS.format('design.xml'))
    else:
        os.system(FIXBB_JD3.format('design.xml'))
    os.chdir(root)


def rename_chain(structure, chain_id):
    renamed_struct = []
    with open(structure, 'r') as pdb:
        pdb_lines = pdb.readlines()
    for line in pdb_lines:
        if line[0:4] != 'ATOM' and line[0:6] != 'HETATM':
            continue
        if line[21].isalpha():
            new_line = list(line)
            new_line[21] = chain_id
            renamed_struct.append("".join(new_line))
    os.remove(structure)
    with open(structure, 'w') as new_structure:
        for new_chain_line in renamed_struct:
            new_structure.write(new_chain_line)


def process_for_piper(pdb_dict, receptor):
    # Create directory
    if not os.path.exists(piper_dir):
        os.makedirs(piper_dir)

    for frag in os.listdir(fixbb_dir):
        if os.path.splitext(frag)[1] == '.pdb':   # Rename chain ID to 'B'
            rename_chain(os.path.join(fixbb_dir, frag), 'B')
            os.system(COPY.format(os.path.join(fixbb_dir, frag), piper_dir))

    with open(os.path.join(piper_dir, 'runs_list'), 'w') as runs_list:  # Create list 1 - 50
        for i in range(1, 51):
            runs_list.write("{:02}\n".format(i))

    # Create symlinks for piper run
    ligands_dir = os.path.join(piper_dir, 'ligands/')
    if not os.path.exists(ligands_dir):
        os.makedirs(ligands_dir)
    lig_inx = 1
    for frag_name in pdb_dict.keys():
        designed_name = frag_name[:-4] + '_0001.pdb'
        os.system('ln -s {} {}/lig.{:04}.pdb'.format(os.path.join(piper_dir, designed_name), ligands_dir, lig_inx))
        lig_inx += 1

    # process ligands for PIPER
    os.chdir(ligands_dir)
    for lig in os.listdir(ligands_dir):
        os.system(PDBPREP.format(lig))
        os.system(PDBNMD.format(lig))

    # prepare receptor for piper
    os.system(COPY.format(receptor, piper_dir))
    os.chdir(piper_dir)

    rename_chain(os.path.basename(receptor), 'A')
    os.system(PDBPREP.format(os.path.basename(receptor)))
    os.system(PDBNMD.format(os.path.basename(receptor)))
    os.chdir(root)


def build_peptide(pep):
    """Create directory for prepacking and, extended peptide and change its chain id to 'B'"""
    if not os.path.exists(prepack_dir):
        os.makedirs(prepack_dir)
    os.chdir(prepack_dir)
    # Build extended peptide
    os.system(BUILD_PEPTIDE.format(pep))

    # Change chain ID to 'B'
    rename_chain('peptide.pdb', 'B')
    os.chdir(root)


def create_batch(receptor, i, script):
    rec_name = os.path.join(piper_dir, receptor + '_nmin.pdb')
    lig_name = 'lig.' + "{:04}".format(i) + '_nmin.pdb'

    if script == 'piper':
        with open('run_piper', 'w') as piper:
            piper.write(SBATCH_PIPER.format(r=rec_name, l=lig_name))
    # elif script == 'decoys':
    #     with open('run_extract_decoys', 'w') as extract_decoys:
    #         extract_decoys.write(SBATCH_EXTRACT_TOP_DECOYS % lig_name)


def run_piper_fpd(receptor):
    os.chdir(piper_dir)
    receptor_name = os.path.splitext(os.path.basename(receptor))[0]
    list_of_runs = []
    for i in range(1, 51):
        run_dir = "{:02}".format(i)
        list_of_runs.append(os.path.basename(run_dir))
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
            os.chdir(run_dir)
            os.system(COPY.format(os.path.join(piper_dir, 'ligands', 'lig.' + "{:04}".format(i) + '_nmin.pdb'),
                                  'lig.' + "{:04}".format(i) + '_nmin.pdb'))

            create_batch(receptor_name, i, 'piper')
            # create_batch(receptor_name, i, 'decoys')
            os.system(RUN_PIPER_SLURM)  # run PIPER docking and extract 250 top decoys
        os.chdir(piper_dir)

    while len(list_of_runs) > 0:
        for results_dir in os.listdir(piper_dir):
            if os.path.basename(results_dir) in list_of_runs:
                os.chdir(results_dir)
                if os.path.exists('ft.000.00'):
                    with open('ft.000.00') as trans_rot_file:
                        i = 1
                        for line in trans_rot_file:
                            piper_model = line.split()[0]
                            os.system(APPLY_FTRESULTS.format(model=piper_model,
                                                             lig='lig.*.pdb', out=line.split()[0]))
                            if i >= 250:
                                break
                            i += 1
                        input_list = prepare_fpd_input(receptor,
                                                       os.path.join(refinement_dir, os.path.basename(results_dir) +
                                                                    '.' + piper_dir))
                        run_refinement(input_list)
                        list_of_runs.remove(os.path.basename(results_dir))
                        os.chdir(piper_dir)

                else:
                    os.chdir(piper_dir)
                    time.sleep(3)
                    continue
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


def prepack_receptor(receptor):
    """Prepack receptor for FlexPepDock run"""
    print("Prepacking receptor")
    os.chdir(prepack_dir)
    combine_receptor_peptide(receptor, 'peptide.pdb', 'start.pdb')
    prepack_flags_file(receptor)
    if talaris:
        os.system(PREPACK_TALARIS)
    else:
        os.system(PREPACK)
    ppk_receptor = os.path.splitext(os.path.basename(receptor))[0] + '.ppk.pdb'
    with open(ppk_receptor, 'w') as renamed_rec:
        with open('start_0001.pdb', 'r') as start:
            for line in start:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                    if line[21] == 'A':
                        renamed_rec.write(line)
                else:
                    continue
    print("Prepack done!")


def prepare_fpd_input(receptor, output):
    """For each model in current dir replace receptor with prepacked one"""
    list_of_models = []  # Here store all the names of the models for fpd run
    if not os.path.exists(refinement_dir):  # Create directory for FPD
        os.makedirs(refinement_dir)
    os.chdir(refinement_dir)  # Enter the directory
    # Define a name of prepacked repceptor
    ppk_receptor = os.path.join(prepack_dir, os.path.splitext(os.path.basename(receptor))[0] + '.ppk.pdb')
    for model in os.listdir('.'):
        if os.path.splitext(model)[1] == '.pdb':
            combine_receptor_peptide(ppk_receptor, model, output)
            os.system('gzip ' + output)
            list_of_models.append(output + '.gz')
    return list_of_models


def run_refinement(input_list):
    """Prepare fpd input: replace receptor in piper models with prepacked one, create .gz files,
    create refinement flags, run refinement"""
    with open('input_list', 'w') as i_list:
        for model in input_list:
            i_list.write(model + '\n')
    refine_flags_file()
    if talaris:
        os.system(FPD_REFINEMENT_TALARIS)
    else:
        os.system(FPD_REFINEMENT)


def run_protocol(peptide_sequence, receptor):

    # make_pick_fragments(peptide_sequence)
    #
    # create_params_file(FRAGS_FILE.format(str(pep_length)))
    #
    # # extract fragments, create resfiles and return a dictionary of fragments names and matching resfiles names
    # pdb_and_resfiles = process_frags(peptide_seq)
    #
    # create_xml(pdb_and_resfiles)  # create xml for running fixbb with JD3
    #
    # # run fixbb
    # run_fixbb()
    #
    # # process ligands and receptor for piper run
    # process_for_piper(pdb_and_resfiles, receptor)

    # build_peptide(os.path.abspath(sys.argv[2]))  # build extended peptide and rename it's chain id to 'B'
    #
    # prepack_receptor(receptor)

    # run piper docking and extract top 250 models
    run_piper_fpd(receptor)

    # TODO: clustering


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print('Usage: \n [receptor.pdb] [peptide_sequence] optional: [restore_talaris_behaviour]'
              '\nYou need to provide a pdb file for receptor and a text file with peptide sequence '
              '(up to 15 amino acids)\n'
              'If you want to run it with talaris2014, add "restore_talaris_behaviour" option and make'
              ' sure you have both 2016 and 2018 versions of Rosetta')
        sys.exit()

    if len(sys.argv) == 4 and sys.argv[3] == 'restore_talaris_behaviour':
        talaris = True
    else:
        talaris = False

    with open(sys.argv[2], 'r') as peptide:
        peptide_seq = peptide.readline().strip()

    root = os.getcwd()

    # Define all the directories that will be created:
    frag_picker_dir = os.path.join(root, 'frag_picker')
    fragments_dir = os.path.join(root, 'top_50_frags')
    resfiles_dir = os.path.join(root, 'resfiles')
    fixbb_dir = os.path.join(root, 'top_50_frags/fixbb')
    piper_dir = os.path.join(root, 'piper')
    prepack_dir = os.path.join(root, 'prepacking')
    refinement_dir = os.path.join(root, 'refinement')

    pep_length = len(peptide_seq)

    receptor_path = os.path.abspath(sys.argv[1])

    run_protocol(peptide_seq, receptor_path)
