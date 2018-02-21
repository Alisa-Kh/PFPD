import sys
import os

# Rosetta directories

ROSETTA_DIR = '/vol/ek/share/rosetta/'
ROSETTA_VERSION = 'rosetta_src_2017.45.59812_bundle'
PATH_TO_ROSETTA = ROSETTA_DIR + ROSETTA_VERSION
ROSETTA_DATABASE = PATH_TO_ROSETTA + '/main/database'


# Commands

PRODY = '/vol/ek/share/csbw/tools/env2017/bin/prody fetch %s'  # used to fetch PDB, but any other method can be used

RENUMBERING = PATH_TO_ROSETTA + '/tools/protein_tools/scripts/pdb_renumber.py {} {}'

EXCISE_PDB = ROSETTA_DIR + '/pdbUtil/excisePdb_v2.pl {} {} {} {} {}'

FIXBB = PATH_TO_ROSETTA + '/main/source/bin/fixbb.linuxgccrelease' \
                          ' -database ' + ROSETTA_DATABASE + ' -in:file:s %s' \
                          ' -resfile resfiles/resfile_%s -ex1 -ex2 -use_input_sc' \
                          ' -scorefile design_score.sc -ignore_zero_occupancy false' \
                          ' >>design.log'

FIXBB_SCRIPTS = PATH_TO_ROSETTA + '/main/source/bin/rosetta_scripts.default.linuxgccrelease -database '\
                + ROSETTA_DATABASE + ' -parser:protocol {}'

BUILD_PEPTIDE = PATH_TO_ROSETTA + '/main/source/bin/BuildPeptide.linuxgccrelease -in:file:fasta {}' \
                                  ' -database ' + ROSETTA_DATABASE + ' -out:file:o peptide.pdb ' \
                                  '> build_peptide.log'

MAKE_FRAGMENTS = 'perl make_fragments.pl -verbose -id xxxxx {} 2>log'
FRAG_PICKER = PATH_TO_ROSETTA + '/main/source/bin/fragment_picker.linuxgccrelease' \
                                ' -database ' + ROSETTA_DATABASE + ' @flags >makeFrags.log'

COPY = 'cp {} {}'

# Other constants

FRAGS_FILE = 'frags.100.{}mers'
CUR_DIR = os.getcwd()


def make_pick_fragments(pep_seq):
    frags_dir = 'frag_picker'
    if not os.path.exists(frags_dir):
        os.makedirs(frags_dir)
    os.system(COPY.format('make_fragments.pl', frags_dir))
    with open(os.path.join(frags_dir, 'xxxxx.fasta'), 'w') as fasta_file:
        fasta_file.write('>|' + pep_seq + '\n' + pep_seq + '\n')
    os.chdir(frags_dir)
    os.system(MAKE_FRAGMENTS.format('xxxxx.fasta'))
    with open('psi_L1.cfg', 'w') as scores:
        scores.write('#score\tname\tpriority\twght\tmin_allowed\textras\n'
                     'SecondarySimilarity\t350\t2.0\t-\tpsipred\n'
                     'ProfileScoreL1\t200\t1.0\t-\n'
                     '#SequenceIdentity\t150\t0.0\t-\n'
                     '#RamaScore\t100\t6.0\t-\tpsipred\n'
                     '#FragmentCrmsd\t30\t0.0\t-\n'
                     '#FragmentAllAtomCrmsd\t20\t0.0\t-')
    with open('flags', 'w') as flags_file:
        flags_file.write('-in:file:vall\t' + ROSETTA_DIR +
                         'rosetta_fragments_latest/nnmake_database/vall.dat.2006-05-05\n'
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
    os.system(COPY.format(FRAGS_FILE.format(pep_length), '../'))
    os.chdir('../')


def create_params_file(frags):

    # Read only needed values from frags_file and store them in frags_file_values
    parameters_sets = []
    frags_parameters = open('50_frags_parameters', 'w+')
    with open(frags) as frags_file:
        all_file_lines = frags_file.readlines()

    top_50_frags_block = 40 + pep_length * 50

    # Get parameters and save them in parameters_sets list
    j = 0
    for i in range(2, top_50_frags_block):
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


def delete_frag(fragment):  # delete bad fragments
    os.remove(fragment)
    print("Wrong length fragment has been deleted")


def review_frags(outfile, start, end):
    # check correctness of fragments
    with open(outfile) as frag:
        residues = set()
        cur_line = frag.readline().split()
        while cur_line[0] != 'ATOM':  # find ATOM lines
            cur_line = frag.readline().split()

            # if there is no atoms...
            if not cur_line:
                delete_frag(outfile)
                return False

        cur_line = frag.readline()
        while cur_line[22:27].strip() != start:
            cur_line = frag.readline()

            # if there is no start residue
            if not cur_line:
                delete_frag(outfile)
                return False

        cur_line = frag.readline()
        while cur_line[22:27].strip() != end:
            residues.add(cur_line[22:27])
            cur_line = frag.readline()

            # if there is no end residue
            if not cur_line:
                delete_frag(outfile)
                return False

        residues.add(cur_line[22:27])
    if len(residues) != (int(end) - int(start) + 1):
        delete_frag(outfile)
        return False

    else:
        return True


def extract_fragments(pep_sequence):

    # Open the frags_parameters, extract and append parameters to different lists
    with open('50_frags_parameters', 'r') as f:
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
    if not os.path.exists('top_50_frags'):
        os.makedirs('top_50_frags')

    # create directory for storing resfiles
    if not os.path.exists('resfiles'):
        os.makedirs('resfiles')

    # for storing resfiles_names and pdbs
    pdb_resfiles_dict = dict()

    os.chdir('top_50_frags')
    # Fetch PDBs (here done with prody), call excisePdb and then delete the pdb file
    for pdb, chain, start, end, sequence in zip(pdbs, chains, start_res, end_res, sequences):
        fragment_name = pdb + '.' + chain + '.' + start + '.' + end
        outfile = fragment_name + '.pdb'
        pdb_full = pdb + '.pdb'

        os.system(PRODY % pdb)

        if os.path.exists(pdb_full):
            os.system(EXCISE_PDB.format(pdb_full, chain, start, end, outfile))
            print("Extracting fragment")
        else:
            print("Failed to fetch pdb")
            continue  # if failed to fetch PDB

        if not os.path.exists(outfile):
            new_pdb = pdb_full
            print("The numeration is incorrect, renumbering...:")
            os.system(RENUMBERING.format(pdb_full, new_pdb))  # Renumber PDB if it doesn't start from 1

            os.system(EXCISE_PDB.format(pdb_full, chain, start, end, outfile))
        os.remove(pdb_full)

        is_frag_ok = review_frags(outfile, start, end)

        if is_frag_ok:
            create_resfile(pep_sequence, chain, start, sequence, fragment_name)
            pdb_resfiles_dict[outfile] = 'resfile_' + fragment_name

    os.chdir('../')
    return pdb_resfiles_dict


def create_resfile(ori_seq, chain, start, sequence, fragment_name):

    path_to_resfile = '../resfiles/resfile_%s'

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


def create_xml(pdb_resfile_dict):
    job_string = '<Job>\n\t<Input>\n\t\t<PDB fillename="{}"/>\n\t</Input>\n' \
                 '\t<TASKOPERATIONS>\n\t\t<ReadResfile name="read_resfile" filename="{}"/>\n' \
                 '\t</TASKOPERATIONS>\n</Job>\n'
    if not os.path.exists('top_50_frags/fixbb'):
        os.makedirs('top_50_frags/fixbb')
    with open('top_50_frags/fixbb/design.xml', 'w') as xml_file:
        xml_file.write('<JobDefinitionFile>\n')
        for pdb, resfile in pdb_resfile_dict.items():
            xml_file.write(job_string.format(pdb, resfile))
        xml_file.write('</JobDefinitionFile>')
    return xml_file


def run_fixbb():

    os.chdir('top_50_frags/fixbb')
    print("running fixbb design")
    os.system(FIXBB_SCRIPTS.format('design.xml'))


def build_peptide(pep_seq):

    # Build extended peptide
    os.system(BUILD_PEPTIDE.format(pep_seq))

    # Change chain ID to 'B'
    renamed_peptide = []
    with open('peptide.pdb', 'r') as pep:
        pdb_lines = pep.readlines()
    for line in pdb_lines:
        if line[21].isalpha():
            new_line = list(line)
            new_line[21] = 'B'
            renamed_peptide.append("".join(new_line))
    os.remove('peptide.pdb')
    with open('peptide.pdb', 'w') as new_peptide:
        for bchain_line in renamed_peptide:
            new_peptide.write(bchain_line)


if __name__ == "__main__":

    with open(sys.argv[1], 'r') as peptide:
        peptide_seq = peptide.readline().strip()

    pep_length = len(peptide_seq)

    make_pick_fragments(peptide_seq)

    create_params_file(FRAGS_FILE.format(str(pep_length)))

    pdb_and_resfiles = extract_fragments(peptide_seq)  # extract fragments, create resfiles
    #  and return a dictionary with fragments names and matching resfiles names

    xml_design = create_xml(pdb_and_resfiles) # create xml for running fixbb with RS

    # run fixbb
    run_fixbb()

    # preprocessing for PIPER
    # run PIPER
    # extract top 250 models

    build_peptide(peptide_seq)  # build extended peptide and change it's chain id to 'B'
    # prepack receptor
    # run refinement
    # clustering
