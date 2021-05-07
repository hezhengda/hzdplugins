# This script can help us calculate the Wannier function
import os
import sys
import subprocess
import re

from hzdplugins.aiidaplugins.utilities import pwInputToDict, dictToPwInput
from hzdplugins.aiidaplugins.info import getOptimizedStructure, getXMLFromPW, getAtomicSpeciesList
from hzdplugins.aiidaplugins.utilities import getSubmitFile

# process the input commands
for ind, value in enumerate(sys.argv):
    if value == '--xml':
        filename_xml = sys.argv[ind+1]
    if value == '--inp':
        filename_inp = sys.argv[ind+1]

# get optimized structure from xml files
data_xml = getXMLFromPW(xml_file=filename_xml)
newCell, newAtomicPositions = getOptimizedStructure(data_xml)
atomic_species_list = getAtomicSpeciesList(data_xml)

# prepare the input files
dict_input = pwInputToDict(loc_file=filename_inp, atomic_species_list=atomic_species_list)

# prepare the scf input file: inp_scf

dict_input['CONTROL']['calculation'] = "'scf'"
dict_input['CONTROL']['outdir'] = "'./'"
dict_input['CONTROL']['prefix'] = "'hzd'"
dict_input['CONTROL']['pseudo_dir'] = "'./'"
dict_input['CONTROL']['restart_mode'] = "'from_scratch'"

dict_input['SYSTEM']['hubbard_u'] = {
    'Ni': 1e-7,
    'Fe': 1e-7
}
dict_input['SYSTEM']['starting_magnetization'] = {
    'Ni': 0.1,
    'Fe': 0.1
}
dict_input['SYSTEM']['nbnd'] = 400

dict_input['CELL_PARAMETERS'] = newCell
dict_input['ATOMIC_POSITIONS'] = newAtomicPositions

dictToPwInput(dict_input=dict_input, location='inp_scf',atomic_species_list=atomic_species_list)

print('input file for scf calculation ... finished!')

# prepare the input file for pmw.x: inp_pmw

first_band = 130 # This number can be systematically tested
last_band = 400

f = open('inp_pmw', 'w')
f.write('&inputpp\n')
f.write("  outdir='./'\n")
f.write("  prefix='hzd'\n")
f.write("  first_band={}\n".format(first_band))
f.write("  last_band={}\n".format(last_band))
f.write('/')
f.close()

print('input file for pmw.x calculation ... finished!')

# prepare the input file: inp_scf_w_wannier

dict_inp_scf_w_wannier = pwInputToDict(loc_file='inp_scf', atomic_species_list=atomic_species_list)

dict_inp_scf_w_wannier['CONTROL']['restart_mode'] = "'restart'"
dict_inp_scf_w_wannier['CONTROL']['tprnfor'] = '.false.'
dict_inp_scf_w_wannier['CONTROL']['tstress'] = '.false.'

dict_inp_scf_w_wannier['SYSTEM']['U_projection_type'] = "'file'"
dictToPwInput(dict_input=dict_inp_scf_w_wannier, location='inp_scf_w_wannier',atomic_species_list=atomic_species_list)

print('input file for scf calculation with wannier function ... finished')

# prepare the input file for projwfc.x: inp_ldos

ngauss = 0
degauss = 0.05
emin = -40.0
emax = 40.0
deltaE = 0.01

f = open('inp_ldos', 'w')
f.write('&projwfc\n')
f.write("  outdir='./'\n")
f.write("  prefix='hzd'\n")
f.write("  lsym=.true.\n") 
f.write("  filpdos='hzd'\n")
f.write("  filproj='hzd.proj.dat'\n")
f.write("  ngauss={}\n".format(ngauss))
f.write("  degauss={}\n".format(degauss))
f.write("  Emin={}\n".format(emin))
f.write("  Emax={}\n".format(emax))
f.write("  DeltaE={}\n".format(deltaE))
f.write("/")
f.close()

print('input file for PDOS calculation ... finished!')

# prepare the submit file for inp_scf
getSubmitFile(filename='job_scf', 
              computer='juwels', 
              typeCalculation='qe', 
              inpDict={
                  'job_name':'scf',
                  'scheduler_stdout':'stdout_scf',
                  'scheduler_stderr':'stderr_scf',
                  'resources': {'num_machines':2},
                  'max_wallclock_seconds':7200,
                  'modules': ['Stages/2020', 'Intel/2020.2.254-GCC-9.3.0', 'ParaStatinMPI/5.4.7-1', 'QuantumESPRESSO/6.6'],
                  'cmd': 'srun pw.x -nk 2 < inp_scf > out_scf'
                  })

print('submission script for scf calculation ... finished!')

# prepare the submit file for inp_pmw
getSubmitFile(filename='job_pmw', 
              computer='juwels', 
              typeCalculation='qe', 
              inpDict={
                  'job_name':'pmw',
                  'scheduler_stdout':'stdout_pmw',
                  'scheduler_stderr':'stderr_pmw',
                  'resources': {'num_machines':1},
                  'max_wallclock_seconds':7200,
                  'modules': ['Stages/2020', 'Intel/2020.2.254-GCC-9.3.0', 'ParaStatinMPI/5.4.7-1', 'QuantumESPRESSO/6.6'],
                  'cmd': 'srun pmw.x < inp_pmw > out_pmw'
                  })

print('submission script for pmw.x calculation ... finished')

# prepare the submit file for inp_scf_w_wannier
getSubmitFile(filename='job_scf_w_wannier', 
              computer='juwels', 
              typeCalculation='qe', 
              inpDict={
                  'job_name':'scf_w_wannier',
                  'scheduler_stdout':'sww_out',
                  'scheduler_stderr':'sww_err',
                  'resources': {'num_machines':1},
                  'max_wallclock_seconds':7200,
                  'modules': ['Stages/2020', 'Intel/2020.2.254-GCC-9.3.0', 'ParaStatinMPI/5.4.7-1', 'QuantumESPRESSO/6.6'],
                  'cmd': 'srun pw.x < inp_scf_w_wannier > out_scf_w_wannier'
                  })

print('submission script for scf calculation with wannier function ... finished!')

# prepare the submit file for inp_ldos 
getSubmitFile(filename='job_ldos', 
              computer='juwels', 
              typeCalculation='qe', 
              inpDict={
                  'job_name':'ldos',
                  'scheduler_stdout':'ldos_out',
                  'scheduler_stderr':'ldos_err',
                  'resources': {'num_machines':1},
                  'max_wallclock_seconds':7200,
                  'modules': ['Stages/2020', 'Intel/2020.2.254-GCC-9.3.0', 'ParaStatinMPI/5.4.7-1', 'QuantumESPRESSO/6.6'],
                  'cmd': 'srun projwfc.x < inp_ldos > out_ldos'
})

print('submission script for PDOS calculation ... finished!')

## regular expression patterns for the output
#pattern = re.compile('([0-9]+)')
#
## do the scf simulation
#print('submitting scf job ...')
#out_str_scf = subprocess.run('sbatch job_scf', capture_output=True)
#tmp = pattern.search(out_str_scf)
#id_str_scp = tmp.group()
#
## do the pmw.x calculation
#print('submitting pmw.x job ...')
#out_str_pmw = subprocess.run('sbatch --dependency=afterok:{} job_pmw'.format(id_str_scp))
#tmp = pattern.search(out_str_pmw)
#id_str_pmw = tmp.group()
#
## use the wannier function as projector
#print('submitting scf with wannier function ...')
#out_str_w_wannier = subprocess.run('sbatch --dependency=afterok:{} job_scf_w_wannier'.format(id_str_pmw))
#tmp = pattern.search(out_scf_w_wannier)
#id_str_w_wannier = tmp.group()
#
## calculate the pdos
#print('submitting PDOS job ...')
#out_str_ldos = subprocess.run('sbatch --dependency=afterok:{} job_ldos'.format(id_str_w_wannier))
