#!/usr/bin/python3

import argparse
import sys
import os
import pfpd_const as pfpd


###########################################################################
"""It is not recommended to change the flags, unless you know what you 
are doing. If you do, you may change the flags in the following function"""
###########################################################################


def fragments_flags_and_cfg(psipred='xxxxx.psipred_ss2', checkpoint='xxxxx.checkpoint', n_frags=100):
    # Create psi_L1.cfg file:
    with open('psi_L1.cfg', 'w') as scores:
        scores.write('#score\tname\tpriority\twght\tmin_allowed\textras\n'
                     'SecondarySimilarity\t350\t2.0\t-\tpsipred\n'
                     'ProfileScoreL1\t200\t1.0\t-\n')
    # Write flags files
    with open('flags', 'w') as flags_file:
        flags_file.write('-in:file:vall\t' + pfpd.ROSETTA_TOOLS +
                         'fragment_tools/vall.jul19.2011.gz\n'
                         '-in:file:checkpoint\t{check}\n'
                         '-frags:describe_fragments\tfrags.fsc\n'
                         '-frags:frag_sizes\t{len}\n'
                         '-frags:n_candidates\t2000\n'
                         '-frags:n_frags\t{n_frags}\n'
                         '-out:file:frag_prefix\tfrags\n'
                         '-frags:ss_pred\t{psi} psipred\n'
                         '-frags:scoring:config\tpsi_L1.cfg\n'
                         '-frags:bounded_protocol\ttrue\n'
                         '-mute\tcore.util.prof\n'
                         '-mute\tcore.conformation\n'
                         '-mute\tcore.chemical\n'
                         '-mute\tprotocols.jumping'.format(check=checkpoint, len=pep_length,
                                                           psi=psipred, n_frags=n_frags))

####################################################
"""ATTENTION!!!!! Do not change the code below!!!"""
####################################################


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


def create_custom_pp(ss_pred, psip_name='xxxxx.psipred_custom', ori_psip='xxxxx.psipred_ss2'):
    with open(psip_name, 'w') as psi_new:
        with open(ori_psip, 'r') as psipred:
            psi_new.write(psipred.readline())  # write a header to a new custom psipred
            psipred_lines = psipred.readlines()
            for i in range(len(psipred_lines)):
                new_line = psipred_lines[i].split()
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
    return psip_name


def create_windows(struct):
    wins = []
    if struct == 'b':
        sign = 'E'
        window = 'E' * pfpd.WINDOWS_LENGTH
    else:  # for alpha-helix
        sign = 'H'
        window = 'H' * pfpd.WINDOWS_LENGTH
    if pfpd.WINDOWS_LENGTH < pep_length:
        for start in range(0, pep_length - pfpd.WINDOWS_LENGTH + 1, 2):
            remainder = pep_length - start - pfpd.WINDOWS_LENGTH
            if remainder >= 2:
                custom_pp = 'C' * start + window + 'C' * remainder
            else:
                custom_pp = 'C' * start + 'E' * (pep_length - start)
            wins.append(custom_pp)
    else:
        wins.append(sign * pep_length)
    return wins


def add_frag_to_list():
    with open('frags.1.{}mers'.format(pep_length), 'a') as one_frag:
        with open(pfpd.FRAGS_FILE.format(pep_length), 'r') as frag_file:
            all_frags = frag_file.readlines()
        for frag_line in all_frags[2:]:
            one_frag.write(frag_line)
        os.system('mv frags.1.{}mers '.format(pep_length) + pfpd.FRAGS_FILE.format(pep_length))


def add_alpha_beta_frags(psip_file, check_file):
    """For 12 aa peptide create 6aas beta and alpha windows with 2 aas leaps.
    For 6aa peptides or shorter - full helices or beta-strands. Call frag_picker with custom psipred for all the
    windows. Max number of additional fragments will be 8, they will replace the last 8 fragments. Minimal number - 2"""
    beta_wins = create_windows('b')
    alpha_wins = create_windows('a')
    for win in beta_wins:
        fragments_flags_and_cfg(create_custom_pp(win, 'b_psipred', psip_file), check_file, 1)
        os.system(pfpd.FRAG_PICKER)
        add_frag_to_list()
    for win in alpha_wins:
        fragments_flags_and_cfg(create_custom_pp(win, 'a_psipred', psip_file), check_file, 1)
        os.system(pfpd.FRAG_PICKER)
        add_frag_to_list()


def make_pick_fragments(pep_seq, ss_pred=None):
    """Run fragment picker"""
    if not os.path.exists(frag_picker_dir):
        os.makedirs(frag_picker_dir)
    # Create fasta file:
    with open(os.path.join(frag_picker_dir, 'xxxxx.fasta'), 'w') as fasta_file:
        fasta_file.write('>|' + pep_seq + '\n' + pep_seq + '\n')
    os.chdir(frag_picker_dir)
    print('************Running BLAST/PsiPred***************')

    os.system(pfpd.MAKE_FRAGMENTS.format('xxxxx.fasta'))  # Run make_fragments.pl script

    if sec_struct:
        fragments_flags_and_cfg(create_custom_pp(ss_pred))  # Create custom psipred file and flag file with its name
        psi_p_file, checkpoint_file = create_custom_pp(ss_pred), 'xxxxx.checkpoint'
    elif full_p:
        psi_p_file = cut_file(1)  # cut psipred_ss2
        checkpoint_file = cut_file(0)  # cut checkpoint
        fragments_flags_and_cfg(psi_p_file, checkpoint_file)  # Write flags files
    else:
        fragments_flags_and_cfg()
        psi_p_file, checkpoint_file = 'xxxxx.psipred_ss2', 'xxxxx.checkpoint'
    print("**************Picking fragments**************")
    os.system(pfpd.FRAG_PICKER)  # Run fragment picker
    add_alpha_beta_frags(psi_p_file, checkpoint_file)
    os.system(pfpd.COPY.format(pfpd.FRAGS_FILE.format(pep_length), root))  # Copy fragments file (frags.100.nmers)
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
            if cur_line[:4] == 'ATOM' and cur_line[17:20] == pfpd.THREE_TO_ONE_AA[sequence[i]]:
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
        os.system(pfpd.CLEAN_PDB.format(pdb, chain))  # get clean pdb and it's fasta

        pdb_full = pfpd.PDB.format(pdb.upper(), chain)  # names of clean_pdb output files
        fasta_name = pfpd.FASTA.format(pdb.upper(), chain)

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
                if frags_count >= frags_num + add_frags_num:
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
    if frags_count < frags_num:
        return frags_num - frags_count
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
    add_frags = add_frags[frags_num + bad_frags_shift:]
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
            os.system(pfpd.FIXBB_JD3_TALARIS.format('design.xml'))
        else:
            for pdb, resfile in pdb_and_resfiles.items():
                os.system(pfpd.FIXBB_TALARIS.format(frag=os.path.join(fragments_dir, pdb),
                                                    resfile=os.path.join(resfiles_dir, resfile)))
    else:
        if jd3:
            create_xml(pdb_and_resfiles)
            os.system(pfpd.FIXBB_JD3.format('design.xml'))
        else:
            for pdb, resfile in pdb_and_resfiles.items():
                os.system(pfpd.FIXBB.format(frag=os.path.join(fragments_dir, pdb),
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


def arg_parser():
    parser = argparse.ArgumentParser(description='Please provide a text file with peptide sequence (or FASTA file).\n'
                                                 'If you want to run the protocol with talaris2014, add '
                                                 '"--restore_talaris_behavior" option and make sure you have'
                                                 ' the 2016 version of Rosetta. For running with jd3 - add'
                                                 ' --jd3 option. If the secondary structure is known add --sec_struct ' 
                                                 'following by a file with secondary structure. If you want  '
                                                 'PsiPred of the full protein (or just longer version of your actual '
                                                 'peptide), which is recommended, add --long_pep_psipred\n'
                                                 'option followed by your full/long protein fasta file')

    parser.add_argument('peptide_sequence')
    parser.add_argument('--restore_talaris_behavior', dest='talaris', action='store_true', default=False)
    parser.add_argument('--jd3', dest='job_distributor', action='store_true', default=False)
    parser.add_argument('--sec_struct', dest='ss_pred', default=None)
    parser.add_argument('--long_pep_psipred', dest='full_prot_seq', default=None)  # pdb_name,chain (e.g. 1abc,A)
    parser.add_argument('--n_frag', dest='num_frags', default=pfpd.NUM_OF_FRAGS)

    return parser


def run_protocol(peptide_sequence):

    if sec_struct:
        with open(sec_struct, 'r') as ss_h:
            ss_pred = ss_h.readline().strip()
        for char in ss_pred:
            if char not in pfpd.PSIPRED_OUTPUT:
                print('Wrong secondary structure file format. A valid file should contain only 1 line with C, H or E '
                      'letters')
                sys.exit()
        make_pick_fragments(peptide_sequence, ss_pred)

    elif full_p:
        create_psipred_from_full_protein(os.path.abspath(full_p))
    else:
        make_pick_fragments(peptide_sequence)

    all_frags = create_params_file(pfpd.FRAGS_FILE.format(str(pep_length)))

    # extract fragments, create resfiles and return a dictionary of fragments names and matching resfiles names
    pdb_and_resfiles = process_frags(peptide_seq, all_frags)

    # run fixbb
    run_fixbb(pdb_and_resfiles)


if __name__ == "__main__":

    arguments = arg_parser().parse_args()

    with open(arguments.peptide_sequence, 'r') as peptide:
        peptide_seq = peptide.readlines()
        if peptide_seq[0][0] == '>':
            peptide_seq = peptide_seq[1].strip()
        else:
            peptide_seq = peptide_seq[0].strip()

    pep_length = len(peptide_seq)

    talaris = arguments.talaris
    jd3 = arguments.job_distributor
    sec_struct = arguments.ss_pred
    full_p = arguments.full_prot_seq
    frags_num = int(arguments.num_frags)

    # Define all the directories that will be created:
    root = os.getcwd()
    frag_picker_dir = os.path.join(root, 'frag_picker')
    fragments_dir = os.path.join(root, 'top_50_frags')
    resfiles_dir = os.path.join(root, 'resfiles')
    fixbb_dir = os.path.join(fragments_dir, 'fixbb')
    bad_frags_dir = 'bad_fragments'

    run_protocol(peptide_seq)
