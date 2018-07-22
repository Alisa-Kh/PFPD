#!/usr/bin/python3

import os

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
                                                                            ' -database ' + ROSETTA_DB + \
                ' -out:file:o peptide.pdb > build_peptide.log'

MAKE_FRAGMENTS = 'perl ' + os.path.join(PFPD_SCRIPTS, 'make_fragments.pl') + \
                 ' -verbose -id xxxxx {} 2>log'

# MAKE_FRAGMENTS = 'perl ' + os.path.join(ROSETTA_TOOLS, 'fragment_tools/make_fragments.pl') + \
#                  ' -old_name_format -verbose -id xxxxx {} 2>log'

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
                       "do if [ ! -f {f}.pdb ]; then " + PIPER_DIR + "apply_ftresult.py -i $f ft.000.00 " \
                       + PIPER_DIR + "prms/rot70k.0.0.4.prm %s --out-prefix $f;fi;done"
APPLY_FTRESULTS = 'python ' + PIPER_DIR + 'apply_ftresult.py -i {model} ft.000.00 ' \
                  + PIPER_DIR + 'prms/rot70k.0.0.4.prm {lig} --out-prefix {out}'

# Other commands

COPY = 'cp {} {}'

PREPARE_FPD_IN = "piper_run=`pwd | awk -F'/' '{print $NF}'`\nfor f in `ls [0-9]*.pdb`;" \
                 "do cat %s ${f} | grep ATOM >%s/${piper_run}.${f};gzip %s/${piper_run}.${f};done"

# Other constants

NUM_OF_FRAGS = 50  # default = 50
PIPER_MODELS_EXTRACTED = str(250)  # default = 250
N_ROTS = str(70000)  # default = 70000
WINDOWS_LENGTH = 6  # default = 6

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

# Common functions


def count_pdbs(folder):
    """Count files in a directory"""
    return len([frag for frag in os.listdir(folder) if
                os.path.isfile(os.path.join(folder, frag)) and
                os.path.splitext(os.path.basename(frag))[1] == '.pdb'])


def count_dirs(folder):
    """Count directories in a directory"""
    return len([d for d in os.listdir(folder) if
                os.path.isdir(os.path.join(folder, d))])


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
