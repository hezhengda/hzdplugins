def pwInputToDict(loc_file):
    """
    Create the input dictionary from the INP_PWSCF

    :param loc_file: The location of the input file for pw.x
    :type loc_file: Python string

    :returns: a python dictionary that contains all the relevant information
    """
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

def dictToPwInput(dict_input, location):
    """
    To convert the dictionary to the input file for pw.x
    
    :param dict_input: contains all the information about the structure
    :type dict_input: Python dictionary

    :param location: contains the location of the output file
    :type location: Python string

    :returns: a file named in location
    """

    f = open(location, 'w')
   
    for key, value in dict_input.items():
    
    if key in ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL']:
        
        f.write('&{}\n'.format(key))
        
        for param, p_v in value.items():
            f.write('  {}={}\n'.format(param, p_v))
        
        f.write('/\n')

    elif key in ['ATOMIC_SPECIES']:
        
        f.write('{}\n'.format(key))
        
        for element, e_v in value.items():
            f.write('{} {} {}\n'.format(element, e_v['atomic_weight'], e_v['pseudopotential']))

    elif key in ['K_POINTS']:
        
        kpts = value['kpoints']
        disp = value['displacement']
        f.write('K_POINTS automatic')
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
        
        f.write('ATOMIC_POSITIONS (angstrom)')
        
        for l in value:
            f.write('{}\n'.format(listToStr(l, ' ')))
    
    else:
        pass

    f.close()
