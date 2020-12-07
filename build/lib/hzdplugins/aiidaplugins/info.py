# This document contains many small functions that can help us find relevant information about each node
# It would be better if the key of results dictionary is the uuid of each node, because uuid is unique, but pk value is not.

import pandas as pd
from aiida.orm import load_node, QueryBuilder, Node
from copy import deepcopy
from hzdplugins.aiidaplugins.constants import results_keys_set
import json

def showResults(results):

    """

    show results in pandas form.

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A panda object that we can use to display the dataframe

    """

    df = pd.DataFrame.from_dict(results,
                                orient='index',
                                columns=results_keys_set)
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    return df

def get_unDoneTasks(results):

    """

    show all the unDone computational tasks, which means the exit_status is neither '0' or '501', or None.

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

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

def get_unFinishedTasks(results):

    """

    show all the unFinished computational tasks, which means the is_finished tag is False

    Parameters:

    results:
        The results dictionary that contains all relevant computational information. The keys are uuid of each nodes

    Return: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

    """

    subresults = {}

    for key, value in results.items():
        node = load_node(uuid=key)
        if node.is_finished == False:
            subresults[key] = value

    return subresults

def get_unConvergedTasks(results):

    """

    show all the finished by not converged computational tasks, which means the `is_finished` tag is True, but `is_finished_ok` is False, although when `exit_status == '501'`, it is still ok.

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

    """

    subresults = {}

    for key, value in results.items():
        node = load_node(uuid=key)
        if node.is_finished == True:
            if (node.is_finished_ok == False):
                if node.exit_status != 501:
                    subresults[key] = value

    return subresults

def assignValue(results):

    """

    Assign the current status of simulaton to results. The function will do two things: (1) clean the current results dictionary, remove any key that does not belong to the new key set. (2) Add the values in the new set.

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A dictionary that has modified and assigned values

    """

    results_tmp = deepcopy(results)

    for uuid_node, value in results.items():
        value_tmp = deepcopy(value)
        node = load_node(uuid=uuid_node)
        # clean the results
        for key in value_tmp:
            if not(key in results_keys_set):
                results_tmp[pk_str].pop(key)
        # assign the value
        results_tmp[uuid_node]['system'] = node.label
        results_tmp[uuid_node]['cluster'] = node.computer.label

        if 'CONTROL' in node.inputs.parameters.get_dict().keys(): # it is pw.x calculation
            results_tmp[uuid_node]['comp_type'] = node.inputs.parameters.get_dict()['CONTROL']['calculation']
            results_tmp[uuid_node]['xc functional'] = node.inputs.parameters.get_dict()['SYSTEM']['input_dft']
            compcode = 'pw'

        if 'projwfc' in node.inputs.parameters.get_dict().keys(): # it is projwfc calculation
            results_tmp[uuid_node]['comp_type'] = 'PDOS'
            compcode = 'projwfc'

        if (node.is_finished_ok) or (node.exit_status == 0) or (node.exit_status == 501):
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
        else:
            results_tmp[uuid_node]['E/eV'] = None
            results_tmp[uuid_node]['is_finished'] = node.is_finished
            results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
            if node.is_killed == True:
                results_tmp[uuid_node]['exit_status'] = 'killed'
            else:
                results_tmp[uuid_node]['exit_status'] = str(node.exit_status)

    return results_tmp

def pkToUuidConverter(results_pk):

    """

    `pkToUuidConverter` can convert the `pk` based results dictionary to `uuid` based dictionary.

    Parameters:

    results_pk:
        The results dictionary that has pk as its key

    Return: A dictionary that has uuid as its key

    """

    results_tmp = {} # our empty and temporary dictionary for storing the converted dictionary

    for key, value in results_pk.items():
        # get the node by uuid, pk may not work (different profile)
        uuid_node = value['uuid']
        qb = QueryBuilder()
        qb.append(Node, filters={'uuid': {'==':uuid_node}})
        node = qb.all()[0][0] # because uuid is unique, it can only return to one node

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

def get_ChargeAndMagneticMoments(uuid):

    """

    :code:`get_ChargeAndMagneticMoments` function will output the charge and the magnetic moments for each atom if we have spin-polarized simulation.

    Parameter:

    uuid:
        The uuid of the computational node.

    Return: A table that shows the charge and atomic_magnetic_moments for each species (labeled as `[number][atomic_species]`, e.g. 1Ni)

    """

    node = load_node(uuid=uuid)
    trajectory = node.outputs.output_trajectory
    atomic_species = trajectory.get_array('atomic_species_name')
    magnetic_moments = trajectory.get_array('atomic_magnetic_moments')[-1] # the last step
    charges = trajectory.get_array('atomic_charges')[-1]

    results = {}
    for i in range(len(atomic_species)):
        name = atomic_species[i] + str(i) # e.g. Ni1
        results[name] = {}
        results[name]['charge'] = charges[i]
        results[name]['magnetic_moment'] = magnetic_moments[i]

    df = pd.DataFrame.from_dict(results, orient='index')

    return df

def get_TotalForces(uuid):

    """

    :code:`get_TotalForces` function will output the total force for each atomic step.

    Parameter:

    uuid:
        The uuid of the computational node.

    Return: A matplotlib figure that shows the convergence, and also the last 5 steps of total_forces.

    """

    node = load_node(uuid=uuid)
    trajectory = node.outputs.output_trajectory
    total_force = trajectory.get_array('total_force')

    iteration = []
    tf = []

    for id in range(len(total_force)):
        iteration.append(id)
        tf.append(total_force[id])

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(iteration, tf)
    plt.xlabel('iterations')
    plt.ylabel('total force / eV/A')
    plt.show()

    return tf[-5:-1]

def saveResults(results, filename):

    """

    `saveResults` can be used for saving results

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    filename:
        The name of the file that you want to store in

    Return: A dictionary that has modified and assigned values

    """

    # dump results file
    with open(filename, 'w') as json_file:
        json.dump(results, json_file)

    return 'Your results have been successfully saved.'

def readResults(filename):

    """

    `readResults` can be used for reading results from json file

    Parameters:

    filename:
        The name of the file that you store all the information in

    Return: A dictionary that has modified and assigned values

    """

    with open(filename, 'r') as json_file:
        results = json.load(json_file)

    return results
