from copy import deepcopy
from aiida.orm import load_node


def qeCleanOneRemoteFolder(results, uuid):
    """

    :code:`qeCleanOneRemoteFolder` This method is used for clean the content on remote folder generated by Quantum
    Espresso. Since QE will generate a lot of files (e.g. 'mix', 'hub', 'wfc', 'restart', etc), so the storage on our
    supercomputer will be overloaded soon, that's not a good thing.

    :param results: the dictionary that contains all the information about the calculation for the project.
    :type results: python dictionary object

    :param uuid: the uuid value of certain node (CalcJob)
    :type uuid: python string object

    :returns: 0 if everything works fine, or if something is wrong, then we need to use try-exception to output the
              error.

    """

    results_tmp = deepcopy(results)

    # check whether uuid is in results_tmp or not
    if not (str(uuid) in results_tmp.keys()):
        print('uuid:{} --- Your input is not in the results dictionary, please check your number again.')
        return results_tmp
    else:
        node = load_node(uuid=uuid)

    authinfo = node.get_authinfo()
    transport = authinfo.get_transport()
    remote_folder_path = node.get_remote_workdir()  # get the path of remote folder, our working directory
    exit_status = node.exit_status
    finished = node.is_finished

    # check whether the calculation is finished (0 or 501) or not
    if (exit_status == 0) or (exit_status == 501) or finished:
        results_tmp[str(uuid)]['exit_status'] = str(exit_status)
        results_tmp[str(uuid)]['is_finished'] = finished
    elif node.is_excepted:
        results_tmp[str(uuid)]['exit_status'] = 'excepted'
        results_tmp[str(uuid)]['is_finished'] = finished
    else:
        print('uuid:{} --- Sorry, your calculation is not finished yet, please wait until it is finished.'.format(uuid))
        return results_tmp

    # check whether the remote folder has already been cleaned, so we won't need to do that again
    if results_tmp[str(uuid)]['remove_remote_folder']:
        print('uuid:{} --- The remote folder has already been cleared.'.format(uuid))
        return results_tmp

    # open transport portal
    transport.open()

    # chdir to remote_folder_path
    transport.chdir(remote_folder_path)  # now we are in working directory

    # The structure of each working directory is the same:
    # aiida.in aiida.out _aiidasubmit.sh out(folder) pseudo(folder)

    # determine whether out folder is in the path
    file_list = transport.listdir()
    if 'out' in file_list:
        transport.chdir(transport.getcwd() + '/out')  # move to the out file
        # start deleting the files
        # *mix* *restart* *hub* *wfc*
        file_list = transport.listdir()
        for item in file_list:
            if 'hub' in item:
                transport.remove(transport.getcwd() + '/' + item)  # delete the *hub* file
            elif 'mix' in item:
                transport.remove(transport.getcwd() + '/' + item)  # delete the *mix* file
            elif 'restart' in item:
                transport.remove(transport.getcwd() + '/' + item)  # delete the *restart* file
            elif 'wfc' in item:
                transport.remove(transport.getcwd() + '/' + item)  # delete the *wfc* file

        # we still have a folder called aiida.save
        if 'aiida.save' in file_list:
            transport.chdir(transport.getcwd() + '/aiida.save')  # move to aiida.save folder
            file_list = transport.listdir()
            # we should delete wfc file and also charge-density file
            for item in file_list:
                if 'wfc' in item:
                    transport.remove(transport.getcwd() + '/' + item)  # delete the wavefunction file
                if 'charge' in item:
                    transport.remove(transport.getcwd() + '/' + item)  # delete the charge-density file
        results_tmp[str(uuid)]['remove_remote_folder'] = True
        print('uuid:{} -- All the unnecessary files have been deleted.'.format(uuid))
        transport.close()
        return results_tmp
    else:
        print(
            'uuid:{} --- There is no out folder in the working directory. Please check whether the calculation is '
            'sucessfully submitted and executed.')
        transport.close()
        return results_tmp


def qecleanAllRemoteFolder(results, *args):
    """

    :code:`qecleanAllRemoteFolder` We would like to clean all the remote folder in the results.

    :param results: A dictionary which key is the uuid (string) of each CalcJob node.
    :type results: python dictionary object

    :param args: you can add a list of nodes that you want to clear the remote folder. e.g. [2020, 2030, 2040] etc. The
                 function will only deal with the list object, other types of inputs are ignored. But if args is not
                 set, then the function will deal with all the nodes in the results dictionary.
    :type args: python list object

    :returns: a modified results dictionary that change all the `remove_remote_folder` to True, which will not be
                  examined in the future.

    """

    results_tmp = deepcopy(results)  # first we need to make a copy of the results dictionary

    if len(args) == 0:  # because args is (), type is tuple
        for key in results_tmp.keys():
            uuid = key
            results_tmp = qeCleanOneRemoteFolder(results_tmp, uuid)
        return results_tmp
    else:
        for arg in args:
            if type(arg) == list:
                for item in arg:
                    uuid = item  # make sure that uuid is integer
                    results_tmp = qeCleanOneRemoteFolder(results_tmp, uuid)
            else:
                pass
        return results_tmp


def qeRetriveAllFiles(uuid, localpath):
    """

    :code:`qeRetriveAllFiles` can help us retrieve all the files from remote_folder_path

    :param uuid: The uuid of the computational node
    :type uuid: python string object

    :param localpath: The absolute path of local folder, which we want to store our information in.
    :type localpath: python string object

    :returns: There is no return.

    """

    import os

    node = load_node(uuid=uuid)

    authinfo = node.get_authinfo()
    transport = authinfo.get_transport()
    remote_folder_path = node.get_remote_workdir()  # get the path of remote folder, our working directory

    transport.open()
    transport.chdir(remote_folder_path)

    os.chdir(localpath)
    localpath_uuid = localpath + '/' + uuid
    os.system('mkdir ' + uuid)

    transport.gettree(remotepath=transport.getcwd(), localpath=localpath_uuid)

    pwd = transport.getcwd()

    transport.close()

    return 'The file from {} have been copied to {}'.format(pwd, localpath_uuid)


def setCmdOnRemoteComputer(cmd, uuid):
    """

    :code:`setCmdOnRemoteComputer` can run command on remote computer from jupyterlab, which is really convenient,
    and also the command is programmable.

    :param cmd: represents the command that needs to be run on the remote server.
    :type cmd: python string object

    :param uuid: represents the uuid of the job. Our cmd will be operated in uuid work directory.
    :type uuid: python string object

    :returns: - **r** (`int`): if 0, means success
              - **stdout** (`list`): the return list of the cmd
              - **stderr** (`str`): the return message of error, if succeeds, then this is empty

    """

    from hzdplugins.aiidaplugins.constants import cmd_shortcut

    # open the transport

    node = load_node(uuid=uuid)

    authinfo = node.get_authinfo()
    transport = authinfo.get_transport()
    remote_folder_path = node.get_remote_workdir()  # get the path of remote folder, our working directory

    transport.open()
    transport.chdir(remote_folder_path)

    if cmd in cmd_shortcut.keys():
        cmd = cmd_shortcut[cmd]

    r, stdout, stderr = transport.exec_command_wait(cmd)

    transport.close()

    return r, stdout, stderr
