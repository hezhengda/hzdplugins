import json
import sys

# read the file line by line
f = open(sys.argv[1],'r')
fJson = open('pwInput.json', 'w')

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

# dump the dictionary
json.dump(input_dict, fJson)
f.close()
fJson.close()
