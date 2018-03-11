#!/usr/bin/python

import sys
import os

# Rosetta directories

##################### Change the ROSETTA_DIR here #####################
ROSETTA_DIR = '/vol/ek/Home/alisa/rosetta/Rosetta/'

ROSETTA_DATABASE = ROSETTA_DIR + 'main/database'

# Commands

GET_PDB = os.path.join(ROSETTA_DIR, 'tools/protein_tools/scripts/clean_pdb.py {} {}')

################## Change number of nodes if needed ####################
FIXBB_JD3 = 'mpirun -n 6 ' + ROSETTA_DIR + 'main/source/bin/fixbb_jd3.mpiserialization.linuxgccrelease' \
            ' -database ' + ROSETTA_DATABASE + ' -restore_talaris_behavior' \
            ' -in:file:job_definition_file {} > fixbb.log'

BUILD_PEPTIDE = ROSETTA_DIR + 'main/source/bin/BuildPeptide.linuxgccrelease -in:file:fasta {}' \
                                  ' -database ' + ROSETTA_DATABASE + ' -out:file:o peptide.pdb ' \
                                  '> build_peptide.log'

MAKE_FRAGMENTS = 'perl /vol/ek/share/scripts/global_pep_dock/fragpicker_setup/make_fragments.pl -verbose' \
                 ' -id xxxxx {} 2>log'

FRAG_PICKER = ROSETTA_DIR + 'main/source/bin/fragment_picker.linuxgccrelease' \
              ' -database ' + ROSETTA_DATABASE + ' @flags >makeFrags.log'

PDBPREP = 'perl /vol/ek/share/scripts/global_pep_dock/PIPER/bin/phplibbin/pdbprep.pl {}'
PDBNMD = 'perl /vol/ek/share/scripts/global_pep_dock/PIPER/bin/phplibbin/pdbnmd.pl "{}"' + ' ?'

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


def make_pick_fragments(pep_seq):
    frags_dir = os.path.join(root, 'frag_picker')
    if not os.path.exists(frags_dir):
        os.makedirs(frags_dir)
    # Create fasta file:
    with open(os.path.join(frags_dir, 'xxxxx.fasta'), 'w') as fasta_file:
        fasta_file.write('>|' + pep_seq + '\n' + pep_seq + '\n')
    os.chdir(frags_dir)
    os.system(MAKE_FRAGMENTS.format('xxxxx.fasta'))  # Run make_fragments.pl script
    # Create psi_L1.cfg file:
    with open('psi_L1.cfg', 'w') as scores:
        scores.write('#score\tname\tpriority\twght\tmin_allowed\textras\n'
                     'SecondarySimilarity\t350\t2.0\t-\tpsipred\n'
                     'ProfileScoreL1\t200\t1.0\t-\n'
                     '#SequenceIdentity\t150\t0.0\t-\n'
                     '#RamaScore\t100\t6.0\t-\tpsipred\n'
                     '#FragmentCrmsd\t30\t0.0\t-\n'
                     '#FragmentAllAtomCrmsd\t20\t0.0\t-')
    # Write flags
    with open('flags', 'w') as flags_file:
        flags_file.write('-in:file:vall\t' + ROSETTA_DATABASE + '/sampling/'
                         'filtered.vall.dat.2006-05-05.gz\n'
                         '-in:file:checkpoint\txxxxx.checkpoint\n'
                         '-frags:describe_fragments\tfrags.fsc\n'
                         '-frags:frag_sizes\t' + str(len(pep_seq)) + '\n'
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
    """Check if fragment length is correct and there are no zero occupancy atoms."""
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
    frags_dir = os.path.join(root, 'top_50_frags')
    if not os.path.exists(frags_dir):
        os.makedirs(frags_dir)

    # create directory for storing resfiles
    if not os.path.exists('resfiles'):
        os.makedirs('resfiles')

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

    path_to_resfile = os.path.join(root, 'resfiles/resfile_%s')

    # Create resfile for each fragment
    resfile = open(path_to_resfile % fragment_name, 'w')
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


def create_xml(pdb_resfile_dict, path_to_fixbb):
    job_string = '<Job>\n\t<Input>\n\t\t<PDB filename="../{}"/>\n\t</Input>\n' \
                 '\t<TASKOPERATIONS>\n\t\t<ReadResfile name="read_resfile" filename="../../resfiles/{}"/>\n' \
                 '\t</TASKOPERATIONS>\n</Job>\n'
    if not os.path.exists(path_to_fixbb):
        os.makedirs(path_to_fixbb)
    with open(path_to_fixbb + '/design.xml', 'w') as xml_file:
        xml_file.write('<JobDefinitionFile>\n')
        xml_file.write('<Common>\n\t<SCOREFXNS>\n\t\t<ScoreFunction name="Talaris14" weights="talaris2014.wts"/>'
                       '\n\t</SCOREFXNS>\n</Common>\n')
        for pdb, resfile in pdb_resfile_dict.items():
            xml_file.write(job_string.format(pdb, resfile))
        xml_file.write('</JobDefinitionFile>')


def run_fixbb(path_to_fixbb):
    os.chdir(path_to_fixbb)
    print("running fixbb design")
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


def process_for_piper(fix_bb_dir, pdb_dict, receptor):
    piper_dir = os.path.join(root, 'PIPER')  # Create directory
    if not os.path.exists(piper_dir):
        os.makedirs(piper_dir)

    for frag in os.listdir(fix_bb_dir):
        if os.path.splitext(frag)[1] == '.pdb':   # Rename chain ID to 'B'
            rename_chain(os.path.join(fix_bb_dir, frag), 'B')
            os.system(COPY.format(os.path.join(fix_bb_dir, frag), piper_dir))

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
    os.system(COPY.format(os.path.join(root, receptor), piper_dir))
    os.chdir(piper_dir)

    rename_chain(receptor, 'A')
    os.system(PDBPREP.format(receptor))
    os.system(PDBNMD.format(receptor))
    os.chdir(root)

    return(piper_dir)


def build_peptide(pep_seq):

    # Build extended peptide
    os.system(BUILD_PEPTIDE.format(pep_seq))

    # Change chain ID to 'B'
    rename_chain('peptide.pdb', 'B')


def run_piper(piper_dir):
    pass


def run_protocol(peptide_sequence, receptor):

    # make_pick_fragments(peptide_sequence)

    # create_params_file(FRAGS_FILE.format(str(pep_length)))

    # extract fragments, create resfiles and return a dictionary of fragments names and matching resfiles names
    pdb_and_resfiles = process_frags(peptide_seq)

    path_to_fixbb = os.path.join(root, 'top_50_frags/fixbb')

    create_xml(pdb_and_resfiles, path_to_fixbb)  # create xml for running fixbb with JD3

    # run fixbb
    run_fixbb(path_to_fixbb)

    # process ligands and receptor for piper run
    piper_dir = process_for_piper(path_to_fixbb, pdb_and_resfiles, receptor)

    run_piper(piper_dir)


if __name__ == "__main__":

    with open(sys.argv[1], 'r') as peptide:
        peptide_seq = peptide.readline().strip()

    root = os.getcwd()
    pep_length = len(peptide_seq)

    run_protocol(peptide_seq, sys.argv[2])

    # TODO: run PIPER
    # TODO: extract top 250 models

    build_peptide(sys.argv[1])  # build extended peptide and rename it's chain id to 'B'
    # TODO: prepack receptor
    # TODO: run refinement
    # TODO: clustering
