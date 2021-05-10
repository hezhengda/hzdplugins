import os 
import sys

for i in sys.argv:
    print('i = {}'.format(i))
    os.system("grep 'Hubbard energy' {}/out_scf_w_wannier".format(i))
