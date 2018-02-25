import sys
import os

# Rosetta directories

ROSETTA_DIR = '/vol/ek/share/rosetta/'
ROSETTA_RELEASE = 'rosetta_src_2017.45.59812_bundle/'
PATH_TO_ROSETTA = ROSETTA_DIR + ROSETTA_RELEASE
ROSETTA_DATABASE = PATH_TO_ROSETTA + 'main/database/'


# Commands

PRODY = '/vol/ek/share/csbw/tools/env2017/bin/prody fetch %s'  # used to fetch PDB, but any other method can be used
# GET_PDB = 'perl ~/rosetta/Rosetta/tools/perl_tools/getPdb.pl {}'
# RENUMBER_PDB = 'python /vol/ek/Home/alisa/rosetta/Rosetta/tools/renumber_pdb.py --pdb={} --output={}'
# CLEAN_PDB = 'python ~/rosetta/Rosetta/tools/protein_tools/scripts/clean_pdb.py {} {} --nopdbout >> log'

SEQUENTIAL_RENUM = ROSETTA_DIR + 'pdbUtil/sequentialPdbResSeq.pl -pdbfile {} > {}'

# EXCISE_PDB = ROSETTA_DIR + 'pdbUtil/excisePdb_v2.pl {} {} {} {} {}'

EXTRACT_PDB = '/vol/ek/share/labscripts/extract_chains_and_range.pl -p {p} -c {c} -r {s}-{e} -o {o} >> extracting.log'

# EXTRACT_CHAINS = '/vol/ek/share/labscripts/extract_chains_and_range.pl -p {p} -a -o {o} >> extracting.log'

GET_FASTA = ''

FIXBB = PATH_TO_ROSETTA + 'main/source/bin/fixbb.linuxgccrelease' \
                          ' -database ' + ROSETTA_DATABASE + ' -in:file:s %s' \
                          ' -resfile resfiles/resfile_%s -ex1 -ex2 -use_input_sc' \
                          ' -scorefile design_score.sc -ignore_zero_occupancy false' \
                          ' >>design.log'

FIXBB_SCRIPTS = PATH_TO_ROSETTA + 'main/source/bin/rosetta_scripts.default.linuxgccrelease -database '\
                + ROSETTA_DATABASE + ' -parser:protocol {}'

BUILD_PEPTIDE = PATH_TO_ROSETTA + 'main/source/bin/BuildPeptide.linuxgccrelease -in:file:fasta {}' \
                                  ' -database ' + ROSETTA_DATABASE + ' -out:file:o peptide.pdb ' \
                                  '> build_peptide.log'

MAKE_FRAGMENTS = 'perl /vol/ek/share/scripts/global_pep_dock/fragpicker_setup/make_fragments.pl -verbose' \
                 ' -id xxxxx {} 2>log'
FRAG_PICKER = PATH_TO_ROSETTA + 'main/source/bin/fragment_picker.linuxgccrelease' \
                                ' -database ' + ROSETTA_DATABASE + ' @flags >makeFrags.log'

COPY = 'cp {} {}'

# Other constants

FRAGS_FILE = 'frags.100.{}mers'
CUR_DIR = os.getcwd()


def make_pick_fragments(pep_seq):
    frags_dir = 'frag_picker'
    if not os.path.exists(frags_dir):
        os.makedirs(frags_dir)
    # os.system(COPY.format('make_fragments.pl', frags_dir))
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


def delete_frag(fragment):  # delete bad fragments
    os.remove(fragment)
    print("The fragment has been deleted because of an incorrect length or zero occupancy")


def review_frags(outfile, start, end):
    # check correctness of fragments
    with open(outfile) as frag:
        residues = set()
        cur_line = frag.readline()
        while 'ATOM' not in cur_line[0:4]:  # find ATOM lines
            cur_line = frag.readline()

            # if there there no atoms...
            if not cur_line:
                delete_frag(outfile)
                return False

        cur_line = frag.readline()  # line with a first atom

        while cur_line[22:27].strip() != start:
            cur_line = frag.readline()

            # if there is no start residue
            if not cur_line:
                delete_frag(outfile)
                return False

        cur_line = frag.readline()  # first atom of the start residue
        while cur_line[22:27].strip() != end:

            # check occupancy
            if cur_line[54:60].strip() == 0.00:
                delete_frag(outfile)
                print("Zero occupancy atoms")
                return False

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


# def convert_mse_to_met(pdb_full):
#     file_name = os.path.basename(pdb_full)
#     os.rename(pdb_full, 'tmp.pdb')
#     with open(file_name, 'w') as new_pdb:
#         with open('tmp.pdb') as old_pdb:
#             for line in old_pdb:
#                 if 'HETATM' in line:
#                     if 'MSE' in line:
#                         editing_line = list(line)
#                         editing_line[0] = 'A'
#                         editing_line[1] = 'T'
#                         editing_line[2] = 'O'
#                         editing_line[3] = 'M'
#                         editing_line[4] = ' '
#                         editing_line[5] = ' '
#                         editing_line[18] = 'E'
#                         editing_line[19] = 'T'
#                         new_line = "".join(editing_line)
#                         new_pdb.write(new_line)
#                     else:
#                         new_pdb.write(line)
#                 else:
#                     new_pdb.write(line)
#     os.remove('tmp.pdb')


def extract_frags(pep_sequence):

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

            # preprocess pdb:
            # convert_mse_to_met(pdb_full)  # should help in some cases

            os.rename(pdb_full, 'tmp.pdb')
            os.system(SEQUENTIAL_RENUM.format('tmp.pdb', pdb_full))  # Renumber PDB if it doesn't start from 1
            os.remove('tmp.pdb')

            if chain == '_':
                chain = 'A'

            print("Extracting fragment")
            os.system(EXTRACT_PDB.format(p=pdb_full, c=chain, s=start, e=end, o=outfile))

            os.remove(pdb_full)

            if os.path.exists(outfile) and os.path.getsize(outfile) > 0:

                if review_frags(outfile, start, end):
                    cur_dir = os.getcwd()
                    frags_count = len([frag for frag in os.listdir('.') if
                                       os.path.isfile(os.path.join(cur_dir, frag))])
                    create_resfile(pep_sequence, chain, start, sequence, fragment_name)
                    pdb_resfiles_dict[outfile] = 'resfile_' + fragment_name
                    if frags_count >= 51:   # there are also extracting log in the folder
                        print("Got top 50 fragments")
                        break
                else:
                    continue

            else:
                print("Failed to extract fragment")
                continue
        else:
            print("Failed to fetch pdb")  # The PDB could be obsolete
            continue  # if failed to fetch PDB

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
    for frag in os.listdir('.'):
        resfile_name = 'resfile_' + os.path.splitext(os.path.basename(frag))[0]
        os.system(FIXBB.format(frag, '../../resfiles/' + resfile_name))
    # os.system(FIXBB_SCRIPTS.format('design.xml'))


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
    #
    # make_pick_fragments(peptide_seq)
    #
    # create_params_file(FRAGS_FILE.format(str(pep_length)))

    # extract fragments, create resfiles and return a dictionary of fragments names and matching resfiles names
    pdb_and_resfiles = extract_frags(peptide_seq)

    # xml_design = create_xml(pdb_and_resfiles)  # create xml for running fixbb with RS
    #
    # # run fixbb
    # run_fixbb()
    #
    # # preprocessing for PIPER
    # # run PIPER
    # # extract top 250 models
    #
    # build_peptide(peptide_seq)  # build extended peptide and change it's chain id to 'B'
    # # prepack receptor
    # # run refinement
    # # clustering
