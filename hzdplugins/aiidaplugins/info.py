# This document contains many small functions that can help us find relevant information about each node
# It would be better if the key of results dictionary is the uuid of each node, because uuid is unique, but pk value
# is not.

import pandas as pd
import numpy as np
from aiida.orm import load_node, QueryBuilder, Node
from copy import deepcopy
from hzdplugins.aiidaplugins.constants import results_keys_set
import json


def showResults(results):
    """

    :code:`showResults` shows results in pandas form.

    :param results: The results dictionary that contains all relevant computational information.
    :type results: python dictionary object

    :returns: A pandas object that we can use to display the dataframe

    """

    df = pd.DataFrame.from_dict(results,
                                orient='index',
                                columns=results_keys_set)
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    return df


def getUnDoneTasks(results):
    """

    :code:`getUnDoneTasks` shows all the unDone computational tasks, which means the exit_status is neither '0' or
    '501', or None.

    :param results: The results dictionary that contains all relevant computational information.
    :type results: python dictionary object

    :returns: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

    """

    subresults = {}

    for key, value in results.items():
        if 'exit_status' in value.keys():
            if value['exit_status'] != '0' and value['exit_status'] != '501':
                subresults[key] = value
        else:
            node = load_node(uuid=key)
            if node.exit_status != '0' and node.exit_status != '501':
                subresults[key] = value

    return subresults


def getUnFinishedTasks(results):
    """

    :code:`getUnFinishedTasks` show all the unFinished computational tasks, which means the is_finished tag is False

    :param results: The results dictionary that contains all relevant computational information.
                    The keys are uuid of each nodes
    :type results: python dictionary object

    :returns: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

    """

    subresults = {}

    for key, value in results.items():
        node = load_node(uuid=key)
        if not node.is_finished:
            subresults[key] = value

    return subresults


def getUnConvergedTasks(results):
    """

    :code:`getUnConvergedTasks` show all the finished by not converged computational tasks, which means the
    `is_finished` tag is True, but `is_finished_ok` is False, although when `exit_status == '501'`, it is still ok.

    :param results: The results dictionary that contains all relevant computational information.
    :type results: python dictionary object

    :returns: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

    """

    subresults = {}

    for key, value in results.items():
        node = load_node(uuid=key)
        if node.is_finished:
            if not node.is_finished_ok:
                if node.exit_status != 501:
                    subresults[key] = value

    return subresults


def assignValue(results):
    """

    :code:`assignValue` will assign the current status of simulaton to results. The function will do two things:
    (1)  clean the current results dictionary, remove any key that does not belong to the new key set. (2) Add the
    values in the new set.

    :param results: The results dictionary that contains all relevant computational information.
    :type results: python dictionary object

    :returns: A dictionary that has modified and assigned values

    """

    results_tmp = deepcopy(results)

    for uuid_node, value in results.items():
        value_tmp = deepcopy(value)
        node = load_node(uuid=uuid_node)
        # clean the results
        for key in value_tmp:
            if not (key in results_keys_set):
                results_tmp[uuid_node].pop(key)
        # assign the value
        results_tmp[uuid_node]['system'] = node.label
        results_tmp[uuid_node]['cluster'] = node.computer.label

        keys = [key.lower() for key in
                node.inputs.parameters.get_dict().keys()]  # get all the keys in lower case, easy for comparison

        if 'control' in keys:  # it is pw.x calculation
            results_tmp[uuid_node]['comp_type'] = node.inputs.parameters.get_dict()['CONTROL']['calculation']
            results_tmp[uuid_node]['xc functional'] = node.inputs.parameters.get_dict()['SYSTEM']['input_dft']
            compcode = 'pw'

        if 'projwfc' in keys:  # it is projwfc calculation
            results_tmp[uuid_node]['comp_type'] = 'PDOS'
            compcode = 'projwfc'

        if 'inputph' in keys:
            results_tmp[uuid_node]['comp_type'] = 'PH'
            compcode = 'ph'

        if node.is_finished_ok or (node.exit_status == 0) or (node.exit_status == 501):
            if compcode == 'pw':
                compType = node.inputs.parameters['CONTROL']['calculation']
                if compType == 'relax' or compType == 'vc-relax':
                    results_tmp[uuid_node]['E/eV'] = node.res.energy
                results_tmp[uuid_node]['is_finished'] = node.is_finished
                results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
                results_tmp[uuid_node]['exit_status'] = str(node.exit_status)
            if compcode == 'projwfc':
                results_tmp[uuid_node]['is_finished'] = node.is_finished
                results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
                results_tmp[uuid_node]['exit_status'] = str(node.exit_status)
            if compcode == 'ph':
                results_tmp[uuid_node]['is_finished'] = node.is_finished
                results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
                results_tmp[uuid_node]['exit_status'] = str(node.exit_status)
        else:
            results_tmp[uuid_node]['E/eV'] = None
            results_tmp[uuid_node]['is_finished'] = node.is_finished
            results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
            if node.is_killed:
                results_tmp[uuid_node]['exit_status'] = 'killed'
            elif node.is_excepted:
                results_tmp[uuid_node]['exit_status'] = 'excepted'
            else:
                results_tmp[uuid_node]['exit_status'] = str(node.exit_status)

    return results_tmp


def pkToUuidConverter(results_pk):
    """

    :code:`pkToUuidConverter` can convert the :code:`pk` based results dictionary to :code:`uuid` based dictionary.

    :param results_pk: The results dictionary that has pk as its key
    :type results_pk: python dictionary object

    :returns: A dictionary that has uuid as its key

    """

    results_tmp = {}  # our empty and temporary dictionary for storing the converted dictionary

    for key, value in results_pk.items():
        # get the node by uuid, pk may not work (different profile)
        uuid_node = value['uuid']
        qb = QueryBuilder()
        qb.append(Node, filters={'uuid': {'==': uuid_node}})
        # node = qb.all()[0][0]  # because uuid is unique, it can only return to one node

        # reconstruct results_tmp, add the checking mechanism
        results_tmp[uuid_node] = {}
        keyset = value.keys()
        if 'system' in keyset:
            results_tmp[uuid_node]['system'] = value['system']
        else:
            results_tmp[uuid_node]['system'] = None

        if 'comp_type' in keyset:
            results_tmp[uuid_node]['comp_type'] = value['comp_type']
        else:
            results_tmp[uuid_node]['comp_type'] = None

        if 'cluster' in keyset:
            results_tmp[uuid_node]['cluster'] = value['cluster']
        else:
            results_tmp[uuid_node]['cluster'] = None

        if 'xc functional' in keyset:
            results_tmp[uuid_node]['xc functional'] = value['xc functional']
        else:
            results_tmp[uuid_node]['xc functional'] = None

        if 'exit_status' in keyset:
            results_tmp[uuid_node]['exit_status'] = value['exit_status']
        else:
            results_tmp[uuid_node]['exit_status'] = None

        if 'is_finished' in keyset:
            results_tmp[uuid_node]['is_finished'] = value['is_finished']
        else:
            results_tmp[uuid_node]['is_finished'] = None

        if 'is_finished_ok' in keyset:
            results_tmp[uuid_node]['is_finished_ok'] = value['is_finished_ok']
        else:
            results_tmp[uuid_node]['is_finished_ok'] = None

        if 'E/eV' in keyset:
            results_tmp[uuid_node]['E/eV'] = value['E/eV']
        else:
            results_tmp[uuid_node]['E/eV'] = None

        if 'remove_remote_folder' in keyset:
            results_tmp[uuid_node]['remove_remote_folder'] = value['remove_remote_folder']
        else:
            results_tmp[uuid_node]['remove_remote_folder'] = None

        # for previous_calc and son_calc, we need to change them to uuid
        if 'previous_calc' in keyset:
            results_tmp[uuid_node]['previous_calc'] = results_pk[value['previous_calc']]['uuid']
        else:
            results_tmp[uuid_node]['previous_calc'] = None

        if ('son_calc' in keyset) and (value['son_calc'] in results_pk.keys()):
            results_tmp[uuid_node]['son_calc'] = results_pk[value['son_calc']]['uuid']
        else:
            results_tmp[uuid_node]['son_calc'] = None

    return results_tmp


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


def getStructureAnalysis(uuid, bond_length=2.5, atom_index=[], is_Metal=False):
    """

    :code:`get_StructureAnalysis` is a function that can analyze the local structure of each atom in the structure.

    :param uuid: The uuid of the node.
    :type uuid: python string object

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

    node = load_node(uuid=uuid)

    if 'CalcJobNode' in node.node_type:  # node is a calcjob node
        if node.is_finished:
            structure = node.outputs.output_structure.get_ase()  # since ase structure is more easy to use
        else:
            structure = node.inputs.structure.get_ase()
    elif 'StructureData' in node.node_type:  # node is a StructureData node
        structure = node.get_ase()
    else:
        raise IOError('You need to input either a CalcJobNode or a StructureData object.')

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
                    check, length, symbol = checkDistance(cell, atom, atom2, bond_length)
                    if (atom_id != id2) and check and (isMetal(atom.symbol)):
                        name_atom2 = atom2.symbol + str(id2) + symbol
                        results[name_atom][name_atom2] = length
        else:
            results[name_atom] = {}
            for id2, atom2 in enumerate(structure):
                check, length, symbol = checkDistance(cell, atom, atom2, bond_length)
                if (atom_id != id2) and check:
                    name_atom2 = atom2.symbol + str(id2) + symbol
                    results[name_atom][name_atom2] = length

    # sort the bond length in the dictionary, easy for comparison:
    for key, value in results.items():
        tmp = {k: v for k, v in sorted(value.items(), key=lambda item: item[1])}
        results[key] = tmp

    return results


def getStructure(uuid):
    """

    :code:`getStructure` can give you the structure by using ase_gui.

    :param uuid: The uuid of the node.
    :type uuid: python string object

    :returns: A ase-gui figure represents the structure, which you can view; An ase structure file which you can
              manipulate later.

    """

    node = load_node(uuid=uuid)

    if 'CalcJobNode' in node.node_type:  # node is a calcjob node
        if node.is_finished:
            structure = node.outputs.output_structure.get_ase()  # since ase structure is more easy to use
        else:
            structure = node.inputs.structure.get_ase()
    elif 'StructureData' in node.node_type:  # node is a StructureData node
        structure = node.get_ase()
    else:
        raise IOError('You need to input either a CalcJobNode or a StructureData object.')

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

    :param set_angular_momentum: A list. [0, 1, 2] means ['s', 'p', 'd'] orbitals, and if we are considering f-electrons, then we need to
                                 add '3' to this list.
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


def saveResults(results, filename):
    """

    :code:`saveResults` can be used for saving results

    :param results: The results dictionary that contains all relevant computational information.
    :type results: python dictionary object

    :param filename: The name of the file that you want to store in.
    :type filename: python string object

    :returns: A dictionary that has modified and assigned values

    """

    # dump results file
    with open(filename, 'w') as json_file:
        json.dump(results, json_file)

    return 'Your results have been successfully saved.'


def readResults(filename):
    """

    :code:`readResults` can be used for reading results from json file

    :param filename: The name of the file that you store all the information in.
    :type filename: python string object

    :returns: A dictionary that has modified and assigned values

    """

    with open(filename, 'r') as json_file:
        results = json.load(json_file)

    return results


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

    :returns: True or False. If the distance is small than bond_length in one of seven conditions, then return True;
              otherwise return False.

    """

    def distance(v1, v2):
        return np.linalg.norm(v1 - v2)

    # in here no matter what x,y,z is, the results are the same, it just easy to write
    x_cell = cell[0]
    y_cell = cell[1]
    z_cell = cell[2]

    for ax in [0, 1, -1]:
        for ay in [0, 1, -1]:
            for az in [0, 1, -1]:
                atom2_modifyPosition = atom2.position + ax * x_cell + ay * y_cell + az * z_cell
                if distance(atom1.position, atom2_modifyPosition) <= bond_length:
                    return True, distance(atom1.position, atom2_modifyPosition), 'ax:{}, ay:{}, az:{}'.format(ax, ay,
                                                                                                              az)

    return False, -1, ''


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
