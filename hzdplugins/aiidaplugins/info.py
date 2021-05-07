# This document contains many small functions that can help us find relevant information about each node
# It would be better if the key of results dictionary is the uuid of each node, because uuid is unique, but pk value
# is not.

import pandas as pd
import numpy as np
from aiida.orm import load_node, QueryBuilder, Node
from copy import deepcopy
from hzdplugins.aiidaplugins.constants import results_keys_set
import json
import qeschema

def getChargeAndMagneticMoments(uuid, traj_index=-1):
    """

    :code:`get_ChargeAndMagneticMoments` function will output the charge and the magnetic moments for each atom if
    we have spin-polarized simulation.

    :param uuid: The uuid of the computational node.
    :type uuid: python string object

    :param traj_index: The index of the structure. -1 means the end structure, but sometimes we need to use -2
                       because the last structure was created by SCF calculation, and for vc-relax simulation,
                       the last step can sometimes be questionable.

    :returns: A table that shows the charge and atomic_magnetic_moments for each species (labeled as :code:`[number][
              atomic_species]`, e.g. 1Ni)

    """

    node = load_node(uuid=uuid)
    trajectory = node.outputs.output_trajectory
    atomic_species = trajectory.get_array('atomic_species_name')
    magnetic_moments = trajectory.get_array('atomic_magnetic_moments')[traj_index]  # the last step
    charges = trajectory.get_array('atomic_charges')[traj_index]

    results = {}
    for i in range(len(atomic_species)):
        name = atomic_species[i] + str(i)  # e.g. Ni1
        results[name] = {}
        results[name]['charge'] = charges[i]
        results[name]['magnetic_moment'] = magnetic_moments[i]

    df = pd.DataFrame.from_dict(results, orient='index')

    return df

def getTotalForces(uuid):
    """

    :code:`get_TotalForces` function will output the total force for each atomic step.

    :param uuid: The uuid of the computational node.
    :type uuid: python string object

    :returns: A matplotlib figure that shows the convergence, and also the last 5 steps of total_forces.

    """

    node = load_node(uuid=uuid)
    trajectory = node.outputs.output_trajectory
    total_force = trajectory.get_array('total_force')

    iteration = []
    tf = []

    for force_id in range(len(total_force)):
        iteration.append(force_id)
        tf.append(total_force[force_id])

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(iteration, tf)
    plt.xlabel('iterations')
    plt.ylabel('total force / eV/A')
    plt.show()

    return tf[-5:-1]

def getStructureAnalysis(structure, bond_length=2.5, atom_index=[], is_Metal=False):
    """

    :code:`get_StructureAnalysis` is a function that can analyze the local structure of each atom in the structure.

    :param structure: The structure that we want to investigate
    :type structure: aiida.orm.StructureData

    :param bond_length: The maximum bond length that we consider as a "neighbour", 2.5 is sufficiently large enough,
                        but if can be adjusted by the user.
    :type bond_length: python float object

    :param atom_index: A list that tell the code which atom that we want to investigate, put the atom_id of the atom
                       in the list.
    :type atom_index: python list object

    :param is_Metal: A boolean variable. If you are only interested in the metal elements, then you put that to True,
                     else False.
    :type is_Metal: python boolean object

    :returns: A dictionary that shows the distance of the central atom with its surrouding atoms. Since metals are
              important, so we mainly focus on Metal atoms. Later maybe I can add a boolean parameters to let the
              user choose.

    """

    structure = structure.get_ase()

    cell = structure.cell

    if len(atom_index) > 0:
        investigate_Atoms = atom_index
    else:
        investigate_Atoms = [i for i in range(len(structure))]  # get all the atoms

    results = {}  # the final dictionary

    # calculate the distance between the atoms
    for atom_id in investigate_Atoms:
        atom = structure[atom_id]
        name_atom = atom.symbol + str(atom_id)
        if is_Metal:
            if isMetal(atom.symbol):
                results[name_atom] = {}
                for id2, atom2 in enumerate(structure):
                    dist_results = checkDistance(cell, atom, atom2, bond_length)
                    for tmp in dist_results:
                        check, length, symbol = tmp
                        if (atom_id != id2) and check and (isMetal(atom.symbol)):
                            name_atom2 = atom2.symbol + str(id2) + symbol
                            results[name_atom][name_atom2] = length
        else:
            results[name_atom] = {}
            for id2, atom2 in enumerate(structure):
                dist_results = checkDistance(cell, atom, atom2, bond_length)
                for tmp in dist_results:
                    check, length, symbol = tmp
                    if (atom_id != id2) and check:
                        name_atom2 = atom2.symbol + str(id2) + symbol
                        results[name_atom][name_atom2] = length

    # sort the bond length in the dictionary, easy for comparison:
    for key, value in results.items():
        tmp = {k: v for k, v in sorted(value.items(), key=lambda item: item[1])}
        results[key] = tmp

    return results

def getStructure(structure):
    """

    :code:`getStructure` can give you the structure by using ase_gui.

    :param structure: The uuid of the node.
    :type structure: aiida.orm.StructureData

    :returns: A ase-gui figure represents the structure, which you can view; An ase structure file which you can
              manipulate later.

    """

    structure = structure.get_ase()

    from ase.visualize import view

    v = view(structure, viewer='ngl')

    # setting the output window for the nglview
    # nglview is a really good tool, and I need to learn more about that.
    v.view.add_ball_and_stick()
    v.view.center()
    v.view.layout.width = '800px'
    v.view.layout.height = '800px'
    v.view.add_label(color='blue', radius=1.0, labelType='text',
                     labelText=[structure[i].symbol + str(i) for i in range(len(structure))], zOffset=2.0,
                     attachment='middle_center')
    v.view.gui_style = 'ngl'

    return v, structure

def viewStructure(structure):

    """

    :param structure: The structure that we want to show, it has to be pymatgen structure
    :type structure: pymatgen.core.structure object

    :returns: a nglview variable to show the structure
    :rtype: nglview object

    """

    from ase.visualize import view
    from pymatgen.io.ase import AseAtomsAdaptor

    structure = AseAtomsAdaptor.get_atoms(structure=structure)

    v = view(structure, viewer='ngl')

    # setting the output window for the nglview
    # nglview is a really good tool, and I need to learn more about that.
    v.view.add_ball_and_stick()
    v.view.center()
    v.view.layout.width = '800px'
    v.view.layout.height = '800px'
    v.view.add_label(color='blue', radius=1.0, labelType='text',
                     labelText=[structure[i].symbol + str(i) for i in range(len(structure))], zOffset=2.0,
                     attachment='middle_center')
    v.view.gui_style = 'ngl'

    return v

def getPdos(uuid, index, is_spin, set_angular_momentum=[0, 1, 2]):
    """

    :code:`getPdos` will give you the PDOS that you want, in order to do the analysis later.

    :param uuid: The uuid of the computational node.
    :type uuid: python string object

    :param index: Shows the atoms that we want to investigate.
    :type index: python list object

    :param is_spin: If is_spin is True, then we need to look for spin-up and spin-down, but if is_spin is false,
                    then we only need to care about the orbital, not the spin.
    :type is_spin: python boolean object

    :param set_angular_momentum: A list. [0, 1, 2] means ['s', 'p', 'd'] orbitals, and if we are considering
                                 f-electrons, then we need to add '3' to this list.
    :type set_angular_momentum: python list object

    :returns: A list and A dictionary. A list is the energy list, which is modified by the Fermi energy. The second
              return is a dictionary. First key is the index of atom in the cell, second key is the angular_momentum,
              third key is the magnetic number, and in there is the pdos of certain magnetic number,
              for each angular_momentum, we will have a 'tot' element which shows the combination of all the magnetic
              number.

              The structure can be represented as:

              .. code-block:: python

                    res = {
                        1: {
                            's': {
                                's':
                                'tot':
                            },
                            'p': {
                                'px':,
                                'py',
                                'pz'
                                'tot':
                            }
                            ...
                        }
                    }
    """

    # for the name of the orbital
    from aiida.tools.data.orbital.realhydrogen import RealhydrogenOrbital as RHO
    import numpy as np

    # preparation
    results = {}
    node = load_node(uuid=uuid)

    keys = node.inputs.parameters.get_dict().keys()

    # check whether the user chooses the correct node.
    if 'projwfc' in keys or 'PROJWFC' in keys:
        pass
    else:
        return ValueError("You should pick a projwfc calculation node.")

    # check whether the calculation is using spin-polarization
    if is_spin:
        projections_up = node.outputs.projections_up
        projections_dw = node.outputs.projections_down
    else:
        pass  # wait till I have more information

    # basic dictionary
    dict_ang_mom = {
        0: 's',
        1: 'p',
        2: 'd',
        3: 'f'
    }

    # get energy
    energy = node.outputs.Dos.get_array('x_array')

    for atom_id in index:

        results[atom_id] = {}

        for angular_momentum in set_angular_momentum:  # 0->s; 1->p; 2->d; 3->f

            am_label = dict_ang_mom[angular_momentum]
            results[atom_id][am_label] = {}

            if is_spin:
                tot_up = np.zeros(len(energy))
                tot_dw = np.zeros(len(energy))
            else:
                pass

            for magnetic_number in range(2 * (angular_momentum + 1) - 1):

                print('Assign the pdos (angular_momentum = {}, magnetic_number = {}) to the atom {}'.format(
                    angular_momentum, magnetic_number, atom_id))

                start_id = atom_id * (2 * angular_momentum - 1) + magnetic_number
                orbital_name = RHO.get_name_from_quantum_numbers(angular_momentum=angular_momentum,
                                                                 magnetic_number=magnetic_number).lower()
                results[atom_id][am_label][orbital_name] = {}

                if is_spin:
                    results[atom_id][am_label][orbital_name]['up'] = projections_up.get_pdos(
                        angular_momentum=angular_momentum, magnetic_number=magnetic_number)[start_id][1]
                    tot_up = tot_up + results[atom_id][am_label][orbital_name]['up']
                    results[atom_id][am_label][orbital_name]['dw'] = projections_dw.get_pdos(
                        angular_momentum=angular_momentum, magnetic_number=magnetic_number)[start_id][1]
                    tot_dw = tot_dw + results[atom_id][am_label][orbital_name]['dw']
                else:
                    pass

            results[atom_id][am_label]['tot'] = {}
            if is_spin:
                results[atom_id][am_label]['tot']['up'] = tot_up
                results[atom_id][am_label]['tot']['dw'] = tot_dw

    return energy, results

def getDos(uuid):
    """

    :code:`getDos` can return the DOS of the system

    :param uuid: The uuid of the projwfc calc node
    :type uuid: python string object

    :returns: - **energy** (`python list object`): list of energy
              - **dos** (`python list object`): dos of energy

    """

    node = load_node(uuid=uuid)

    energy = node.outputs.Dos.get_array('x_array')
    dos = node.outputs.Dos.get_array('y_array_0')

    return energy, dos

def getLastScf(uuid):

    """

    :code:`getLastScf` will return the value of the last scf accuracy and the stdout of :code:`grep 'scf accuracy'
    aiida.out`

    :param uuid: The uuid of the simulation
    :type uuid: python string object

    :returns: * The last force of the simulation
              * The stdout of the cmd

    """

    from hzdplugins.aiidaplugins.io import setCmdOnRemoteComputer
    r, stdout, stderr = setCmdOnRemoteComputer(cmd="grep 'scf accuracy' aiida.out", uuid=uuid)

    list_scf_accuracy = []
    strings = stdout.split('\n')
    for string in strings:
        if string != '':
            scf = string.split(' ')[-2]
            list_scf_accuracy.append(scf)

    if len(list_scf_accuracy) == 0: # no scf deteced
        return -1, stdout
    else:
        return list_scf_accuracy[-1], stdout

def getLastForce(uuid):

    """

    :code:`getLastForce` will return the value of the last total force and the stdout of :code:`grep 'Total force'
    aiida.out`

    :param uuid: The uuid of the simulation
    :type uuid: python string object

    :returns: * The last force of the simulation
              * The stdout of the cmd

    """

    from hzdplugins.aiidaplugins.io import setCmdOnRemoteComputer
    r, stdout, stderr = setCmdOnRemoteComputer(cmd="grep 'Total force' aiida.out", uuid=uuid)

    list_force = []
    strings = stdout.split('\n')
    for string in strings:
        if string != '':
            tmp = string.split(' ')
            list_force.append(tmp[12])

    if len(list_force) == 0: # no force detected
        return -1, stdout
    else:
        return list_force[-1], stdout

# In the below are some functions that may help

def checkDistance(cell, atom1, atom2, bond_length):
    """

    :code:`check_distance` function can help us determine whether the two atoms in the slab structure are close
    enough (distance < bond_length) or not. In here we should notice that all slab structures have periodic boundary
    condition (PBC), which means that not only we need to consider the position of atom2, we also need to consider 6
    different atom positions that is in translational symmetry with atom2.

    :param cell: A 3x3 array. The cell parameters of the slab, which can be easily accessed by :code:`structure.cell`
    :type cell: python list object

    :param atom1: atom1
    :type atom1: ase.Atom object

    :param atom2: atom2
    :type atom2: ase.Atom object

    :param bond_length: The threshold of the bond length that we are interested in, can be set by the user.
    :type bond_length: python float object

    :returns: * If true, then function will return [True, distance, add_string]
              * If false, then function will return [False, -1, '']

    """

    def distance(v1, v2):
        return np.linalg.norm(v1 - v2)

    # in here no matter what x,y,z is, the results are the same, it just easy to write
    x_cell = cell[0]
    y_cell = cell[1]
    z_cell = cell[2]

    results = []

    for ax in [0, 1, -1]:
        for ay in [0, 1, -1]:
            for az in [0, 1, -1]:
                atom2_modifyPosition = atom2.position + ax * x_cell + ay * y_cell + az * z_cell
                if distance(atom1.position, atom2_modifyPosition) <= bond_length:
                    results.append((True, distance(atom1.position, atom2_modifyPosition), 'ax:{}, ay:{}, az:{}'.format(ax, ay,
                                                                                                              az)))
                    # return True, distance(atom1.position, atom2_modifyPosition), 'ax:{}, ay:{}, az:{}'.format(ax, ay,
                    #                                                                                           az)
    if len(results) != 0:
        return results
    else:
        return [(False, -1, '')]

def isMetal(atom_symbol):

    """

    :code:`isMetal` can determine whether an element is metal or not.

    :param atom_symbol: The symbol of the atom.
    :type atom_symbol: python string object

    :returns: If it is metal, the return True, else return False.

    """
    if atom_symbol in ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne',
                       'Si', 'P', 'S', 'Cl', 'Ar', 'Ge', 'As', 'Se',
                       'Br', 'Kr', 'Sb', 'Te', 'I', 'Xe', 'Rn']:
        return False
    else:
        return True

def getGCN(structure, bond_length, atom_list, cn_max):

    """
    This function is used to calculate the GCN (generalized coordination number of certain atom) of certain site.

    :param structure: Structure of the slab
    :type structure: aiida.orm.StructureData

    :param bond_length: the maximum bond_length that we need to investigate
    :type bond_length: python real object

    :param atom_list: a list of atoms which constructs a unique site (could be ensemble)
    :type atom_list: python list object

    :param cn_max: the maximum coordination number of certain atom in the structure
    :type cn_max: python integer

    :param atom_name: the name of the atom (useful in dealing with the dictionary)
    :type atom_name: python string object
    """

    import re

    # get the coordination dictionary
    coord_dict = getStructureAnalysis(structure=structure, atom_index=atom_list, bond_length = bond_length, is_Metal=False)

    # do the analysis

    conn_set = [] # for all the atoms that connected to the site
    for symbol, conn_dict in coord_dict.items():
        for name, distance in conn_dict.items():
            match = re.match(r"([a-z]+)(\d+)([a-z]*)", name, re.I)
            if match:
                items = match.groups()
            index = int(items[1]) # get the index of that atom
            conn_set.append(index)

    sum_cn = 0
    for ind in conn_set:
        coord_dict = getStructureAnalysis(structure=structure, atom_index=[ind], bond_length=bond_length, is_Metal=False)
        for key, value in coord_dict.items():
            sum_cn += len(value)

    return sum_cn/cn_max

def checkTotalSpin(symbol, electron_configuration):
    """
    Get the spin of certain configurations.
    """
    results = 0
    for line in electron_configuration:
        if line[1] == 1:
            results += line[3]*1.0
        else:
            results += line[3]*(-1.0)
    return results

def getXMLFromPW(xml_file):
    """
    Get the dictionary from the xml output.

    :param xml_file: The location of the .xml file that we need
    :returns: return the python dictionary that contains all the information
    """

    pw_document = qeschema.PwDocument()
    pw_document.read(xml_file)
    dict_data = pw_document.to_dict()
    
    return dict_data

def getOptimizedStructure(data_xml):
    """
    Get the last structure (atomic positions and cell) from the xml data (parsed as dictionary)

    :param data_xml: The python dictionary parsed from the xml file generated by pw.x
    :returns: There are two returns, the first is the cell of the optimized structure, the second is the atomic positions.
    """

    convertBohrToAngstrom = 0.5291772109

    cell = data_xml['qes:espresso']['output']['atomic_structure']['cell']
    tmp_cell = []
    for key, value in cell.items():
        np_value = np.array(value) * convertBohrToAngstrom
        tmp_cell.append(np_value.tolist())
    
    atomic_positions = data_xml['qes:espresso']['output']['atomic_structure']['atomic_positions']['atom']
    tmp_atomic_positions = []
    for atom in atomic_positions:
        tmp_list = [atom['@name'], atom['$'][0]*convertBohrToAngstrom, atom['$'][1]*convertBohrToAngstrom, atom['$'][2]*convertBohrToAngstrom]
        tmp_atomic_positions.append(tmp_list)
    
    return tmp_cell, tmp_atomic_positions
    
def getAtomicSpeciesList(data_xml):
    
    tmp_species_list = data_xml['qes:espresso']['output']['atomic_species']['species']

    species_list = {}

    for ind, l in enumerate(tmp_species_list):
        species_list[l['@name']] = ind+1
    
    return species_list
