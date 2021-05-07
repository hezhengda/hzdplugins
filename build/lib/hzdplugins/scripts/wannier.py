# This script can help us calculate the Wannier function
import os
import sys

from hzdplugins.aiidaplugins.utilities import pwInputToDict, dictToPwInput
from hzdplugins.aiidaplugins.info import getOptimizedStructure, getXMLFromPW, getAtomicSpeciesList

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

# prepare the input file: inp_scf_w_wannier

dict_inp_scf_w_wannier = pwInputToDict(loc_file='inp_scf', atomic_species_list=atomic_species_list)

dict_inp_scf_w_wannier['CONTROL']['restart_mode'] = "'restart'"
dict_inp_scf_w_wannier['CONTROL']['tprnfor'] = '.false.'
dict_inp_scf_w_wannier['CONTROL']['tstress'] = '.false.'

dict_inp_scf_w_wannier['SYSTEM']['U_projection_type'] = "'file'"
dictToPwInput(dict_input=dict_inp_scf_w_wannier, location='inp_scf_w_wannier',atomic_species_list=atomic_species_list)

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

# do the scf simulation

# do the pmw.x calculation

# use the wannier function as projector

# calculate the pdos
