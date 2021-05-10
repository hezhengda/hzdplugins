from .constants import slurm_options
from copy import deepcopy

def pwInputToDict(loc_file, atomic_species_list):
    """
    Create the input dictionary from the INP_PWSCF

    :param loc_file: The location of the input file for pw.x
    :type loc_file: Python string
    
    :param atomic_species_list: contains the order of the atomic species, very useful for the hubbard_u, etc.
    :type atomic_species_list: Python dictionary

    :returns: a python dictionary that contains all the relevant information
    """
    
    # reverse the atomic_species_list
    r_asl = {v:k for k,v in atomic_species_list.items()}
        
    # read the file line by line
    f = open(loc_file,'r')

    # some basic keywords
    levelOneParameters = ['&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL', 'ATOMIC_SPECIES', 'K_POINTS', 'CELL_PARAMETERS', 'ATOMIC_POSITIONS']

    # process each line, and add them to the dictionary
    input_dict = {}
    for line in f.readlines():
        line = line.rstrip()
        keyword = line.split(' ')[0]

        # face /, just move to the next line
        if keyword == '/':
            continue

        # if meet keyword, then start a new key 
        if keyword in levelOneParameters:
            
            if '&' in keyword:
                keyword = keyword.split('&')[1]
                
                input_dict[keyword] = {}
                currentKeyword = keyword # current keyword
        
            else: 
                input_dict[keyword] = {}
                currentKeyword = keyword
                
                if keyword in 'CELL_PARAMETERS':
                    input_dict[keyword] = []

                elif keyword in 'ATOMIC_POSITIONS':
                    input_dict[keyword] = []

        # if it is not the keyword
        else:
            
            lineContent = line.split(' ')
            tmp_list = []
            for item in lineContent:
                if item != '':
                    tmp_list.append(item)
            
            if currentKeyword in ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL']:
                
                if 'hubbard_u' in tmp_list[0]:
                                        
                    if 'hubbard_u' not in list(input_dict[currentKeyword].keys()):
                        input_dict[currentKeyword]['hubbard_u'] = {}
                        input_dict[currentKeyword]['hubbard_u'][r_asl[int(tmp_list[0][-2])]] = tmp_list[2]
                    else:
                        input_dict[currentKeyword]['hubbard_u'][r_asl[int(tmp_list[0][-2])]] = tmp_list[2]
                        
                elif 'starting_magnetization' in tmp_list[0]:
                    
                    if 'starting_magnetization' not in list(input_dict[currentKeyword].keys()):
                        input_dict[currentKeyword]['starting_magnetization'] = {}
                        input_dict[currentKeyword]['starting_magnetization'][r_asl[int(tmp_list[0][-2])]] = tmp_list[2]
                    else:
                        input_dict[currentKeyword]['starting_magnetization'][r_asl[int(tmp_list[0][-2])]] = tmp_list[2]
                    
                elif 'starting_ns_eigenvalue' in tmp_list[0]:
                    
                    if 'starting_ns_eigenvalue' not in list(input_dict[currentKeyword].keys()):
                        input_dict[currentKeyword]['starting_ns_eigenvalue'] = []
                        snm_list = [
                            int(tmp_list[0][-6]),
                            int(tmp_list[0][-4]),
                            r_asl[int(tmp_list[0][-2])],
                            tmp_list[2]
                        ]
                        input_dict[currentKeyword]['starting_ns_eigenvalue'].append(snm_list)
                    else:
                        snm_list = [
                            int(tmp_list[0][-6]),
                            int(tmp_list[0][-4]),
                            r_asl[int(tmp_list[0][-2])],
                            tmp_list[2]
                        ]
                        input_dict[currentKeyword]['starting_ns_eigenvalue'].append(snm_list)
                    
                else:
                    input_dict[currentKeyword][tmp_list[0]] = tmp_list[2]

            elif currentKeyword in 'ATOMIC_SPECIES':
                input_dict[currentKeyword][tmp_list[0]] = {
                            'atomic_weight': tmp_list[1],
                            'pseudopotential': tmp_list[2]
                        }
            
            elif currentKeyword in 'K_POINTS':
                input_dict[currentKeyword] = {
                            'kpoints': [int(item) for item in tmp_list[0:3]],
                            'displacement': [int(item) for item in tmp_list[3:]]
                        }
            
            elif currentKeyword in 'CELL_PARAMETERS':
                tmp_list = [float(item) for item in tmp_list] 
                input_dict[currentKeyword].append(tmp_list)

            elif currentKeyword in 'ATOMIC_POSITIONS':
                tmp_list = [tmp_list[0]] + [float(item) for item in tmp_list[1:]]
                input_dict[currentKeyword].append(tmp_list)

    f.close()
    return input_dict

def listToStr(l, separator):
    """
    :param l: the list that we want to print
    :type l: Python list

    :param separator: the separator that we want to add, e.g. ' ' (space), etc.
    :type separator: Python string
    """

    tmp_str = str(l[0])
    for i in l[1:-1]:
        tmp_str += separator + str(i)
    tmp_str += separator + str(l[-1])

    return tmp_str

def dictToPwInput(dict_input, location, atomic_species_list):
    """
    To convert the dictionary to the input file for pw.x
    
    :param dict_input: contains all the information about the structure
    :type dict_input: Python dictionary

    :param location: contains the location of the output file
    :type location: Python string

    :param atomic_species_list: contains the order of the atomic species, very useful for the hubbard_u, etc.
    :type atomic_species_list: Python dictionary

    :returns: a file named in location
    """

    f = open(location, 'w')
   
    for key, value in dict_input.items():
    
        if key in ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL']:
            
            f.write('&{}\n'.format(key))
            
            for param, p_v in value.items():

                if 'hubbard_u' in param:
                    
                    for k, v in p_v.items():
                        f.write('  hubbard_u({}) = {}\n'.format(atomic_species_list[k], v))
                
                elif 'starting_magnetization' in param:

                    for k, v in p_v.items():
                        f.write('  starting_magnetization({}) = {}\n'.format(atomic_species_list[k], v))
                
                elif 'starting_ns_eigenvalue' in param:

                    for l in p_v:
                        f.write('  starting_ns_eigenvalue({},{},{}) = {}\n'.format(l[0], l[1], atomic_species_list[l[2]], l[3]))
                
                else:

                    f.write('  {} = {}\n'.format(param, p_v))
            
            f.write('/\n')

        elif key in ['ATOMIC_SPECIES']:
            
            f.write('{}\n'.format(key))
            
            for element, e_v in value.items():
                f.write('{} {} {}\n'.format(element, e_v['atomic_weight'], e_v['pseudopotential']))

        elif key in ['K_POINTS']:
            
            kpts = value['kpoints']
            disp = value['displacement']
            f.write('K_POINTS automatic\n')
            f.write('{} {} {} {} {} {}\n'.format(
                kpts[0], kpts[1], kpts[2],\
                disp[0], disp[1], disp[2]
            ))

        elif key in ['CELL_PARAMETERS']:
            f.write('CELL_PARAMETERS (angstrom)\n')
            tmp_str = ''
            for l in value:
                f.write('{}\n'.format(listToStr(l, ' ')))

        elif key in ['ATOMIC_POSITIONS']:
            
            f.write('ATOMIC_POSITIONS (angstrom)\n')
            
            for l in value:
                f.write('{}\n'.format(listToStr(l, ' ')))
        
        else:
            pass

    f.close()

def getSubmitFile(filename, computer, typeCalculation, inpDict):
    """
    This function can help us generate the submitting file for the supercomputer.

    :param location: the filename of the submitting script
    :type location: Python string

    :param computer: the name of the computer that we want to run on, the options are: 
                    ['rwth-claix', 'juwels-mac', 'jureca-dc-mac', 'jureca-booster-mac'] (currently)
    :type computer: Pythong string
        
    :param typeCalculation: the type of the simulation that we want to conduct
    :type typeCalculation: Python string (e.g. ['qe'] or others)

    :param inpDict: The dictionary that contains all the information that we need for constructing the input file.
                    All the keywords are: ['job_name', 'scheduler_stdout', 'scheduler_stderr', 
                                           'queue_name', 'account', 'resources', 
                                           'max_wallclock_seconds', 'modules', 'cmd', 'ntasks_per_node']
    :type inpDict: Python dictionary

    :return: a file named filename, no other return
    """

    f = open(filename, 'w')

    tmpDict = deepcopy(slurm_options[computer+'-mac'][typeCalculation])
    for k, v in inpDict.items():
        if k == 'resources':
            for k1, v1 in v.items():
                tmpDict[k][k1] = v1
        else:
            tmpDict[k] = v
    
    f.write('#!/bin/bash\n')

    for k, v in tmpDict.items():

        # choose different things to f.write
        if k == 'job_name':
            f.write('#SBATCH --job-name={}\n'.format(v))

        if k == 'scheduler_stdout':
            f.write('#SBATCH --output={}\n'.format(v))

        if k == 'scheduler_stderr':
            f.write('#SBATCH --error={}\n'.format(v))

        if k == 'queue_name':
            f.write('#SBATCH --partition={}\n'.format(v))

        if k == 'account':
            f.write('#SBATCH --account={}\n'.format(v))

        if k == 'resources':
            f.write('#SBATCH --nodes={}\n'.format(v['num_machines']))

            if 'gpus' in list(v.keys()):
                f.write('#SBATCH --gres:gpu={}\n'.format(v['gpus']))
        
        if k == 'ntasks_per_node':
            f.write('#SBATCH --ntasks-per-node={}\n'.format(v))

        if k == 'max_wallclock_seconds':
            
            hours = v // 3600
            minutes = v % 3600 // 60
            seconds = v % 3600 % 60
            
            f.write('#SBATCH --time={:02d}:{:02d}:{:02d}\n'.format(hours,minutes,seconds))

        if k == 'use':

            for item in v:
                f.write('module use {}\n'.format(item))

        if k == 'modules':

            for item in v:
                f.write('module load {}\n'.format(item))

        if k == 'cmd':

            f.write('{}\n'.format(v))

    f.close()