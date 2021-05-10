#################################################
# Author: Zheng-Da He
# Date: 28/04/2021
# Location: Aachen, Germany
#################################################

import dpdata
import os

# find all files with end of xml
allXmlFiles = []
for filename in os.listdir('.'):
    if os.path.isfile(filename) and'xml' in filename:
        allXmlFiles.append(filename)

print(allXmlFiles)

# create the file
for name in allXmlFiles:
    tmp_list = name.split('.')
    os.system('mkdir '+tmp_list[0])
    d_xml = dpdata.LabeledSystem(name, fmt='vasp/xml')
    d_xml.to_deepmd_raw(tmp_list[0]+'/')
    print("{} is done!".format(name))

# combine the file
f_box = open('./box.raw', 'w')
f_coord = open('./coord.raw', 'w')
f_energy = open('./energy.raw', 'w')
f_force = open('./force.raw', 'w')

for name in allXmlFiles:
    tmp_list = name.split('.')
    f_tmp_box = open('{}/box.raw'.format(tmp_list[0]), 'r')
    for line in f_tmp_box.readlines():
        f_box.write(line)
    f_tmp_box.close()
    print('copy box from {} to the box.raw'.format(tmp_list[0]))

    f_tmp_coord = open('{}/coord.raw'.format(tmp_list[0]), 'r')
    for line in f_tmp_coord.readlines():
        f_coord.write(line)
    f_tmp_coord.close()
    print('copy coord from {} to the coord.raw'.format(tmp_list[0]))
   
    f_tmp_energy = open('{}/energy.raw'.format(tmp_list[0]), 'r')
    for line in f_tmp_energy.readlines():
        f_energy.write(line)
    f_tmp_energy.close()
    print('copy energy from {} to the energy.raw'.format(tmp_list[0]))
    
    f_tmp_force = open('{}/force.raw'.format(tmp_list[0]), 'r')
    for line in f_tmp_force.readlines():
        f_force.write(line)
    f_tmp_force.close()
    print('copy force from {} to the force.raw'.format(tmp_list[0]))

# copy the type.raw and type_map.raw
os.system('cp {}/type.raw {}/type_map.raw .'.format(allXmlFiles[0].split('.')[0],allXmlFiles[0].split('.')[0]))

# delete all the tmp folders
for name in allXmlFiles:
    tmp_list = name.split('.')
    os.system('rm -rf '+tmp_list[0])

f_box.close()
f_coord.close()
f_energy.close()
f_force.close()
