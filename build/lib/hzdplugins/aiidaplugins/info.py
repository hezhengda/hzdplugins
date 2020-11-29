# This document contains many small functions that can help us find relevant information about each node

import pandas as pd
from aiida.orm import load_node

def uuid(pk):

    """

    Give the uuid of certain node

    Parameters:

    pk:
        The pk of the node

    Return: uuid of that node

    """

    node = load_node(pk)
    print(node.uuid)

def showresults(results):

    """

    show results in pandas form.

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A panda object that we can use to display the dataframe

    """

    df = pd.DataFrame.from_dict(results, orient='index')
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    return df

def unDoneTasks(results):

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
            pk = int(key)
            node = load_node(pk)
            if node.exit_status != '0' and node.exit_status != '501':
                subresults[key] = value

    return subresults

def unFinishedTasks(results):

    """

    show all the unFinished computational tasks, which means the is_finished tag is False

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

    """

    subresults = {}

    for key, value in results.items():
        pk = int(key)
        node = load_node(pk)
        if node.is_finished == False:
            subresults[key] = value

    return subresults

def unConvergedTasks(results):

    """

    show all the finished by not converged computational tasks, which means the `is_finished` tag is True, but `is_finished_ok` is False, although when `exit_status == '501'`, it is still ok.

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.

    """

    subresults = {}

    for key, value in results.items():
        pk = int(key)
        node = load_node(pk)
        if node.is_finished == True:
            if (node.is_finished_ok == False):
                if node.exit_status != 501:
                    subresults[key] = value

    return subresults
