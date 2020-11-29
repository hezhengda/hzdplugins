from copy import deepcopy
from aiida.orm import load_node

def qeCleanOneRemoteFolder(results, pk):

    """

    This method is used for clean the content on remote folder generated by Quantum Espresso. Since QE will generate a lot of files (e.g. 'mix', 'hub', 'wfc', 'restart', etc), so the storage on our supercomputer will be overloaded soon, that's not a good thing.

    Parameters:

    results:
        the dictionary that contains all the information about the calculation for the project.

    pk:
        the pk value of certain node (CalcJob)

    Return: 0 if everything works fine, or if something is wrong, then we need to use try-exception to output the error.

    """

    results_tmp = deepcopy(results)

    # check whether pk is in results_tmp or not
    if not(str(pk) in results_tmp.keys()):
        print('pk:{} --- Your input is not in the results dictionary, please check your number again.')
        return results_tmp
    else:
        node = load_node(pk)

    authinfo = node.get_authinfo()
    transport = authinfo.get_transport()
    remote_folder_path = node.get_remote_workdir() # get the path of remote folder, our working directory
    exit_status = node.exit_status
    finished = node.is_finished

    # check whether the calculation is finished (0 or 501) or not
    if (exit_status == 0) or (exit_status == 501) or finished:
        results_tmp[str(pk)]['exit_status'] = str(exit_status)
        results_tmp[str(pk)]['is_finished'] = finished
    else:
        print('pk:{} --- Sorry, your calculation is not finished yet, please wait until it is finished.'.format(pk))
        return results_tmp

    # check whether the remote folder has already been cleaned, so we won't need to do that again
    if results_tmp[str(pk)]['remove_remote_folder'] == True:
        print('pk:{} --- The remote folder has already been cleared.'.format(pk))
        return results_tmp

    # open transport portal
    transport.open()

    # chdir to remote_folder_path
    transport.chdir(remote_folder_path) # now we are in working directory

    # The structure of each working directory is the same:
    # aiida.in aiida.out _aiidasubmit.sh out(folder) pseudo(folder)

    # determine whether out folder is in the path
    list = transport.listdir()
    if 'out' in list:
        transport.chdir(transport.getcwd()+'/out') # move to the out file
        # start deleting the files
        # *mix* *restart* *hub* *wfc*
        list = transport.listdir()
        for item in list:
            if 'hub' in item:
                transport.remove(transport.getcwd()+'/'+item) # delete the *hub* file
            elif 'mix' in item:
                transport.remove(transport.getcwd()+'/'+item) # delete the *mix* file
            elif 'restart' in item:
                transport.remove(transport.getcwd()+'/'+item) # delete the *restart* file
            elif 'wfc' in item:
                transport.remove(transport.getcwd()+'/'+item) # delete the *wfc* file

        # we still have a folder called aiida.save
        if 'aiida.save' in list:
            transport.chdir(transport.getcwd()+'/aiida.save') # move to aiida.save folder
            list = transport.listdir()
            # we should delete wfc file and also charge-density file
            for item in list:
                if 'wfc' in item:
                    transport.remove(transport.getcwd()+'/'+item) # delete the wavefunction file
                if 'charge' in item:
                    transport.remove(transport.getcwd()+'/'+item) # delete the charge-density file
        results_tmp[str(pk)]['remove_remote_folder'] = True
        print('pk:{} -- All the unnecessary files have been deleted.'.format(pk))
        return results_tmp
    else:
        print('pk:{} --- There is no out folder in the working directory. Please check whether the calculation is sucessfully submitted and executed.')
        return results_tmp

def qecleanAllRemoteFolder(results, *args):

    """

    We would like to clean all the remote folder in the results.

    Parameters:

    results:
        A dictionary which key is the pk (string) of each CalcJob node

    args:
        you can add a list of nodes that you want to clear the remote folder. e.g. [2020, 2030, 2040] etc. The function will only deal with the list object, other types of inputs are ignored. But if args is not set, then the function will deal with all the nodes in the results dictionary.

    Return: a modified results dictionary that change all the `remove_remote_folder` to True, which will not be examined in the future.

    """

    results_tmp = deepcopy(results) # first we need to make a copy of the results dictionary

    if len(args) == 0: # because args is (), type is tuple
        for key in results_tmp.keys():
            pk = int(key)
            results_tmp = qeCleanOneRemoteFolder(results_tmp, pk)
        return results_tmp
    else:
        for arg in args:
            if type(arg) == list:
                for item in arg:
                    pk = int(item) # make sure that pk is integer
                    results_tmp = qeCleanOneRemoteFolder(results_tmp, pk)
            else:
                pass
        return results_tmp
