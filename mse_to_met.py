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