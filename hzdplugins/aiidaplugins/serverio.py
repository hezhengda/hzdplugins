from aiida.orm import load_node

def qeCleanOneRemoteFolder(uuid, is_wfc=False, is_save=False):
    """

    :code:`qeCleanOneRemoteFolderNoWfc` This method is used for clean the content on remote folder generated by Quantum
    Espresso. Since QE will generate a lot of files (e.g. 'mix', 'hub', 'wfc', 'restart', etc), so the storage on our
    supercomputer will be overloaded soon, that's not a good thing.

    :param uuid: the uuid value of certain node (CalcJob)
    :type uuid: python string object

    :param is_wfc: (optional, default = False) True if we want to delete all the wavefunction files
    :type is_wfc: python boolean object

    :param is_save: (optional, default = False )True if we want to delete all contents in :code:`aiida.save` folder
    :type is_save: python boolean object

    """

    node = load_node(uuid=uuid)

    authinfo = node.get_authinfo()
    transport = authinfo.get_transport()
    remote_folder_path = node.get_remote_workdir()  # get the path of remote folder, our working directory

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
            if is_wfc:
                if 'wfc' in item:
                    transport.remove(transport.getcwd() + '/' + item)  # delete the *wfc* file

        # we still have a folder called aiida.save
        if is_save:
            if 'aiida.save' in file_list:
                transport.chdir(transport.getcwd() + '/aiida.save')  # move to aiida.save folder
                file_list = transport.listdir()
                # we should delete wfc file and also charge-density file
                for item in file_list:
                    if 'wfc' in item:
                        transport.remove(transport.getcwd() + '/' + item)  # delete the wavefunction file
                    if 'charge' in item:
                        transport.remove(transport.getcwd() + '/' + item)  # delete the charge-density file
                    if 'UPF' in item:
                        transport.remove(transport.getcwd() + '/' + item)  # delete the pseudopotential file
        print('uuid:{} -- All the unnecessary files have been deleted.'.format(uuid))
        transport.close()
    else:
        print(
            'uuid:{} --- There is no out folder in the working directory. Please check whether the calculation is '
            'sucessfully submitted and executed.')
        transport.close()


def qeCleanAllRemoteFolder(uuid_list, is_wfc, is_save):
    """

    :code:`qeCleanAllRemoteFolderWithWfc` We would like to clean all the remote folder in the results.

    :param uuid_list: you can add a list of nodes that you want to clear the remote folder. e.g. [2020, 2030,
                      2040] etc. The function will only deal with the list object, other types of inputs are ignored.
                      But if args is not set, then the function will deal with all the nodes in the results dictionary.
    :type uuid_list: python list object

    :param is_wfc: (optional, default = False) True if we want to delete all the wavefunction files
    :type is_wfc: python boolean object

    :param is_save: (optional, default = False )True if we want to delete all contents in :code:`aiida.save` folder
    :type is_save: python boolean object

    """
    for uuid in uuid_list:
        qeCleanOneRemoteFolderWithWfc(uuid, is_wfc, is_save)

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
