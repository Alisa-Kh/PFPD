#!/usr/bin/python

import sys
import os
import pfpd_protocol as protocol
import math

radius = sys.argv[1]
native = sys.argv[2]
silent_input = sys.argv[3]

HEADER = 'DecoyID\t\tClusterN MemberID\tI_sc\tReweighted_sc'
FINAL_DIR = '../../FINAL_RESULTS'
CLUSTERING_DIR = os.getcwd()


def create_pdb_list():
    with open('../score.sc', 'r') as score_file:
        scores = score_file.readlines()
        header = scores[1].split()
        scores = scores[2:]  # SEQUENCE line + header out
        reweighted_column = header.index('reweighted_sc')
        description_column = header.index('description')

        total_decoys = len(scores)
        clustering_pool_num = int(total_decoys / 100)

    with open('pdb_list', 'w') as pdbs:
        i = 0
        for line in sorted(scores, key=lambda sc_line: float(sc_line.split()[reweighted_column])):
            pdbs.write(line.split()[description_column] + '\n')
            i += 1
            if i >= clustering_pool_num:
                break
    return clustering_pool_num


def define_actual_r():
    with open(native, 'r') as n:
        calphas = 0
        pep_calphas = 0
        for line in n:
            if line[13:15] == 'CA':
                calphas += 1
                if line[21] == 'B':
                    pep_calphas += 1
    actual_r = math.sqrt(float(pep_calphas) / float(calphas)) * float(radius)
    print('Actual radius is ' + str(actual_r))
    return actual_r


def run_clustering(actual_r):
    os.system(os.path.join(protocol.ROSETTA_BIN, 'cluster.linuxgccrelease') +
              ' -in:file:silent {decoys} -in:file:silent_struct_type binary -database {DB}'
              ' -cluster:radius {actualR} -in:file:fullatom -tags `cat pdb_list`'
              ' -silent_read_through_errors > clog'.format(decoys='../decoys.silent',
                                                           DB=protocol.ROSETTA_DB, actualR=actual_r))


def results_processing(clustering_pool_num):
    """Create score score lists"""
    print('Now printing results')
    with open('clog', 'r') as log:
        whole_log = log.readlines()
    with open('cluster_list', 'w') as cluster_lst:
        decoys_lst = []
        final_lines_inx = whole_log.index('Timing: \n')
        relevant_lines = whole_log[final_lines_inx - clustering_pool_num + 1:final_lines_inx]
        for line in relevant_lines:
            cluster_lst.write(line.split()[2] + '\t' + line.split()[3] + '\t' + line.split()[4] + '\n')
            decoys_lst.append(line.split()[2] + '\t' + line.split()[3] + '\t' + line.split()[4])
    with open('../score.sc', 'r') as score_file:
        scores = score_file.readlines()
        header = scores[1].split()
        scores = scores[2:]
    interface_sc_idx = header.index('I_sc')
    reweighted_sc_idx = header.index('reweighted_sc')
    description_idx = header.index('description')
    with open('cluster_list_sc', 'w') as cluster_lst_sc:
        cluster_lst_sc.write(HEADER + '\n')
        for decoy_line in decoys_lst:
            for scores_line in scores:
                if scores_line.split()[description_idx] == decoy_line.split()[0]:
                    cluster_lst_sc.write(decoy_line + '\t' + scores_line.split()[interface_sc_idx] + '\t' +
                                         scores_line.split()[reweighted_sc_idx] + '\n')

    os.system('echo "{}" >cluster_list_I_sc_sorted'.format(HEADER))
    os.system('sort -nk 4 cluster_list_sc | sort -u -k2,2 | '
              'sort -nk 4 | head -20 >>cluster_list_I_sc_sorted')

    os.system('echo "{}" >cluster_list_reweighted_sc_sorted'.format(HEADER))
    os.system('sort -nk 5 cluster_list_sc | sort -u -k2,2 | '
              'sort -nk 5 | head -20 >>cluster_list_reweighted_sc_sorted')


def collect_results():
    if not os.path.exists(FINAL_DIR):
        os.makedirs(FINAL_DIR)
    os.system(protocol.COPY.format(os.path.join(CLUSTERING_DIR, 'cluster_list_reweighted_sc_sorted'),
                                   FINAL_DIR))
    with open(os.path.join(FINAL_DIR, 'cluster_list_reweighted_sc_sorted')) as top_reweighted:
        top_reweighted.readline()
        for i in range(10):
            cur_line = top_reweighted.readline().split()
            struct = 'c.{cluster}.{member}.pdb'.format(cluster=cur_line[1], member=cur_line[2])
            os.system(protocol.COPY.format(os.path.join(CLUSTERING_DIR, struct), FINAL_DIR))

if __name__ == "__main__":
    os.rename('../score.sc', '../score_before_rescoring.sc')
    os.rename('../decoys.slilent', '../old_decoys.silent')
    os.rename('../rescore.sc', '../score.sc')
    os.rename('../decoys.rescored.silent', '../decoys.silent')
    clustering_pool = create_pdb_list()
    actual_radius = define_actual_r()
    run_clustering(actual_radius)
    results_processing(clustering_pool)
    collect_results()
