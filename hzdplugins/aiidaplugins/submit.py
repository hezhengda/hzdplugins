# Since the goal of Aiida is to "preservation of provenance", so I should write this plugin with the same goal as
# Aiida, that's why all the relevant information will be provided in aiida.orm types. (I can convert them in the
# code, it is very straightforward)

from copy import deepcopy

from aiida.engine.launch import submit
from aiida.orm import KpointsData, Dict, Code
from aiida.orm import load_node
from aiida.orm.nodes.data.upf import get_pseudos_from_structure

from hzdplugins.aiidaplugins.constants import slurm_options, pwParameter, projwfcParameter, phParameter

def qePwOriginalSubmit(codename, structure, kpoints, pseudo_family, metadata, pseudo_dict={}, add_parameters={},
                       del_parameters={}, cluster_options={}, settings_dict={}):
    """

    :code:`qePwOriginalSubmit` will submit an original computational task to the desired computer by using certain code.

    :param codename: (mandatory) A string represents the code for pw.x that you want to use.
    :type codename: python string object

    :param structure: (mandatory) The structure of your system.
    :type structure: aiida.orm.StructureData object

    :param add_parameters: (optional, default = {}) The desired parameters that you want to state, it can be
                           incomplete, because inside the function there is a default setting for parameters which can
                           be used in most cases, but if you have specific need, you can put that in parameters,
                           the format is similar as pw.x input file.

                           If you want to assign DFT+U and spin-polarization, you need to specify it on your own.

                           In Aiida, there is a very efficient way to specify the :code:`hubbard_u`,
                           :code:`starting_magnetization` and :code:`starting_ns_eigenvalue`. I give some examples
                           in below:

                           .. code-block:: python

                                # hubbard_u
                                'SYSTEM': {
                                    'hubbard_u': {
                                        'Fe': 5.0,
                                        'Fe3': 5.0 # if you have different spins of same atom, then you should use
                                        newStructure function to create the structure
                                    },
                                    'starting_magnetization': {
                                        'Fe': 0.1,
                                        'Fe3': 0.1,
                                    },
                                    'starting_ns_eigenvalue': [
                                        [1, 1, 'Fe', 1.0] # represent: starting_ns_eigenvalue(1, 1, 1)=1.0
                                        # others are the same, if you want to assign to Fe3, just replace Fe with Fe3.
                                    ]
                                }

    :type add_parameters: python dictionary

    :param del_parameters: (optional, default = {}) The tags that we would like to delete, for example if we do not
                           want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure
                           as add_parameters.

                           e.g. :code:`{'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}`
    :type del_parameters: python dictionary object

    :param kpoints: (mandatory) The kpoints that you want to use, if the kpoints has only 1 list, then it is the kpoint
                    mesh, but if two lists are detected, then the first will be k-point mesh, the second one will be the
                    origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]
    :type kpoints: python list object

    :param pseudo_family: (mandatory) The pseudopotential family that you want to use. Make sure that you already have
                          that configured, otherwise an error will occur.
    :type pseudo_family: python string object.

    :param pseudo_dict: (optional, default = {}) which contains the pseudopotential files that we want to use in the
                        simulation. In here it is very important to note that the path of the pseudopotential file
                        has to be in the absolute path.

                        e.g.

                        .. code-block:: python

                            pseudo_dict = {
                                'Fe': UpfData(absolute_path),
                                'Fe3': UpfData(absolute_path)
                            }
    :type pseudo_dict: python dictionary object.

    :param cluster_options: (optional, default = {}) The detailed option for the cluster. Different cluster may have
                            different settings. Only the following 3 keys can have effects: (1) resources (2) account
                            (3) queue_name
    :type cluster_options: python dictionary object

    :param metadata: (mandatory) The dictionary that contains information about metadata. For example: label and
                     description. label and description are mendatory.

                     e.g. :code:`{'label':{}, 'description':{}}`
    :type metadata: python dictionary object

    :param settings_dict: (optional, default = {}) which contains the additional information for the pw.x
                          calculation. e.g. Fixed atom, retrieving more files, parser options, etc. And the
                          command-line options.
    :type settings_dict: python dictionary object

    :returns: uuid of the new CalcJobNode

    """

    code = Code.get_from_string(codename)
    computer = codename.split('@')[1]  # get the name of the cluster
    pw_builder = code.get_builder()

    # pseudopotential
    # check whether pseudo_family and pseudo_dict are set at the same time, if true, then break
    if len(pseudo_family) > 0 and len(pseudo_dict) > 0:
        return ValueError("You cannot set pseudo_family and pseudo_dict at the same time")
    if len(pseudo_family) == 0 and len(pseudo_dict) == 0:
        return ValueError("You need to specify at least one in pseudo_family or pseudo_dict.")

    if len(pseudo_family) != 0:
        pw_builder.pseudos = get_pseudos_from_structure(structure, family_name=pseudo_family)
    if len(pseudo_dict) != 0:
        pw_builder.pseudos = pseudo_dict

    # set kpoints
    kpts = KpointsData()
    if len(kpoints) == 1:
        kpts.set_kpoints_mesh(mesh=kpoints[0])
    else:
        kpts.set_kpoints_mesh(mesh=kpoints[0], offset=kpoints[1])

    # parameters
    parameters_default = Dict(dict=pwParameter)

    # add parameters in add_parameters
    parameters_tmp = deepcopy(parameters_default)

    for key, value in add_parameters.items():
        for key2, value2 in value.items():
            parameters_tmp[key][key2] = value2

    # delete parameters in del_parameters
    for key, value in del_parameters.items():
        tmp = parameters_tmp[key]
        for key2 in value:
            if key2 in tmp.keys():
                tmp.pop(key2)
            else:
                pass

    parameters_default = parameters_tmp

    # set labels and description
    pw_builder.metadata.label = metadata['label']
    pw_builder.metadata.description = metadata['description']

    # set default options for slurm
    pw_builder.metadata.options['resources'] = slurm_options[computer]['qe']['resources']  # in here machine = node
    pw_builder.metadata.options['max_wallclock_seconds'] = slurm_options[computer]['qe'][
        'max_wallclock_seconds']  #in here machine = node
    pw_builder.metadata.options['account'] = slurm_options[computer]['qe']['account']  # in here machine = node
    pw_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['qe']['scheduler_stderr']
    pw_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['qe']['scheduler_stderr']
    pw_builder.metadata.options['queue_name'] = slurm_options[computer]['qe']['queue_name']

    # revised by cluster_options
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            pw_builder.metadata.options['resources'] = cluster_options['resources']
        if 'account' in cluster_options.keys():
            pw_builder.metadata.options['account'] = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            pw_builder.metadata.options['queue_name'] = cluster_options['queue_name']

    # initialize the settings_dict
    if len(settings_dict) == 0:
        settings_dict['cmdline'] = ['-nk', '4']
    else:
        pass  # do nothing

    # get atomic occupations
    if 'lda_plus_u' in parameters_default['SYSTEM']:
        if parameters_default['SYSTEM']['lda_plus_u']:
            settings_dict['parser_options'] = {'parse_atomic_occupations': True}

    # launch the simulation
    pw_builder.structure = structure
    pw_builder.kpoints = kpts
    pw_builder.parameters = parameters_default
    pw_builder.settings = Dict(dict=settings_dict)
    calc = submit(pw_builder)

    return calc.uuid

def qePwContinueSubmit(uuid, pseudo_family, pseudo_dict={}, codename='', parent_folder=True, add_parameters={},
                       del_parameters={}, kpoints=[], cluster_options={}, metadata={}, settings_dict={}):
    """

    `qePwContinueSubmit` will continue a simulation with similar or modified input parameters. All the parameters are
    listed in the kwargs.

    :param uuid: (mandatory) The uuid of previous calculation. We will start our calculation from there. Because uuid
                 is the unique identification number for each CalcJobNode

                    **Notice**: The uuid must be in the results dictionary, if not the program will shout KeyError.
                    And if you are testing, you could use assignValue to quickly create a dictionary that contains
                    the uuid that you want to continue.
    :type uuid: python string object

    :param pseudo_family: (mandatory) The pseudopotential family that you want to use. Make sure that you already have
                          that configured, otherwise an error will occur. This is mendatory.
    :type pseudo_family: python string object

    :param pseudo_dict: (optional, default = {}) Which contains the pseudopotential files that we want to use in the
                        simulation.
    :type pseudo_dict: python dictionary object

    :param codename: (optional, default = '') Represent the code for pw.x that you want to use. If you want to use the
                     same as previous calculation, then you need to use Str('')
    :type codename: python string object

    :param parent_folder: (optional, default = True) If parent_folder is True, then the calculation will start with the
                          output files from previous calculations.
    :type parent_folder: python boolean object

    :param add_parameters: (optional, default = {}) The desired parameters that you want to state, it can be incomplete,
                           because inside the function there is a default setting for parameters which can be used in
                           most cases, but if you have specific need, you can put that in parameters, the format is
                           similar as pw.x input file.

                           If you want to assign DFT+U and spin-polarization, you need to specify it on your own.

                           e.g. :code:`{'CONTROL':{}, 'SYSTEM':{}}`

                           **Notice**: more options in qePwOriginalSubmit function. In qePwContinueSubmit,
                           we assume that the user wants to restart from previous converged wave functions and
                           charge density, so we set ['CONTROL']['restart_mode']='restart', ['ELECTRON'][
                           'startingwfc']='file and ['ELECTRON']['startingpot']='file'.
    :type add_parameters: python dictionary object

    :param del_parameters: (optional, default = {})The tags that we would like to delete, for example if we do not
                           want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure as
                           add_parameters.

                           e.g. :code:`{'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}`
    :type del_parameters: python dictionary object

    :param kpoints: (optional, default = []), if you want to keep the k-points for previous calculation, just use an
                    empty list :code:`[]`. The kpoints that you want to use, if the kpoints has only 1 list,
                    then it is the kpoint mesh, but if two lists are detected, then the first will be k-point mesh,
                    the second one will be the origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]
    :type kpoints: python list object

    :param cluster_options: (optional, default = {}) The detailed option for the cluster. Different cluster may have
                            different settings. Only the following 3 keys can have effects: (1) resources (2)
                            account (3) queue_name. If value is :code:`{}`, then it means we will use previous settings
    :type cluster_options: python dictionary object

    :param metadata: (optional, default = {}) The dictionary that contains information about metadata. For example:
                     label and description.label and description are mendatory. If value is :code:`{}`,
                     then it means we will use previous settings.
    :type metadata: python dictionary object

    :param settings_dict: (optional, default = {}) which contains the additional information for the pw.x calculation.
                          e.g. Fixed atom, retrieving more files, parser options, etc. And the command-line options.
                          If value is :code:`{}`, then it means we will use previous settings.
    :type settings_dict: python dictionary object

    :returns: uuid of the CalcJobNode of the newest calculation.

    """

    node = load_node(uuid=uuid)

    if len(codename) == 0:  # not going to change cluster
        computer = node.computer.label
        restart_builder = node.get_builder_restart()  # get the restart_builder
    else:
        computer = codename.split('@')[1]
        code = Code.get_from_string(codename)
        restart_builder = code.get_builder()

    parameters_tmp = deepcopy(node.inputs.parameters)

    parameters_dict = parameters_tmp.get_dict()
    calc_type = parameters_dict['CONTROL']['calculation']

    # change the parameters (since this is the continuation of the previous calculation)
    parameters_tmp['CONTROL']['restart_mode'] = 'restart'
    parameters_tmp['ELECTRONS']['startingwfc'] = 'file'  # from wave function in aiida.save
    parameters_tmp['ELECTRONS']['startingpot'] = 'file'  # from charge density in aiida.save

    if calc_type == 'relax' or calc_type == 'vc-relax':
        structure = node.outputs.output_structure
    elif calc_type == 'scf' or calc_type == 'nscf':
        structure = node.inputs.structure

    # assign parameters in add_parameters
    for key, value in add_parameters.items():
        for key2, value2 in value.items():
            parameters_tmp[key][key2] = value2

    # delete parameters in del_parameters
    for key, value in del_parameters.items():
        tmp = parameters_tmp[key]
        for key2 in value:
            if key2 in tmp.keys():
                tmp.pop(key2)

    parameters_default = parameters_tmp

    # reset the kpoints
    if len(kpoints) > 0:
        kpts = KpointsData()
        if len(kpoints) == 1:
            kpts.set_kpoints_mesh(mesh=kpoints[0])
        else:
            kpts.set_kpoints_mesh(mesh=kpoints[0], offset=kpoints[1])
    else:
        kpts = node.inputs.kpoints

    # pseudopotential
    # check whether pseudo_family and pseudo_dict are set at the same time, if true, then break
    if len(pseudo_family) > 0 and len(pseudo_dict) > 0:
        return ValueError("You cannot set pseudo_family and pseudo_dict at the same time")
    if len(pseudo_family) == 0 and len(pseudo_dict) == 0:
        return ValueError("You need to specify at least one in pseudo_family or pseudo_dict.")

    if len(pseudo_family) != 0:
        restart_builder.pseudos = get_pseudos_from_structure(structure, family_name=pseudo_family)
    if len(pseudo_dict) != 0:
        restart_builder.pseudos = pseudo_dict

    # set default options for slurm
    restart_builder.metadata.options['resources'] = slurm_options[computer]['qe']['resources']  # in here machine = node
    restart_builder.metadata.options['max_wallclock_seconds'] = slurm_options[computer]['qe'][
        'max_wallclock_seconds']  # in here machine = node
    restart_builder.metadata.options['account'] = slurm_options[computer]['qe']['account']  # in here machine = node
    restart_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['qe']['scheduler_stderr']
    restart_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['qe']['scheduler_stderr']
    restart_builder.metadata.options['queue_name'] = slurm_options[computer]['qe']['queue_name']

    # reset cluster_options:
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            restart_builder.metadata.options['resources'] = cluster_options['resources']
        if 'account' in cluster_options.keys():
            restart_builder.metadata.options['account'] = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            restart_builder.metadata.options['queue_name'] = cluster_options['queue_name']

    # reset metadata
    if len(metadata) > 0:
        if 'label' in metadata.keys():
            restart_builder.metadata.label = metadata['label']
        else:
            restart_builder.metadata.label = node.label

        if 'description' in metadata.keys():
            restart_builder.metadata.description = metadata['description']
        else:
            restart_builder.metadata.description = node.description
    else:
        restart_builder.metadata.label = node.label
        restart_builder.metadata.description = node.description

    # assign the parent_folder
    if parent_folder:
        restart_builder.parent_folder = node.outputs.remote_folder

    # set settings_dict
    if len(settings_dict) > 0:
        pass
    else:
        settings_dict = node.inputs.settings.get_dict()

    # submit the calculation
    restart_builder.structure = structure
    restart_builder.kpoints = kpts
    restart_builder.parameters = parameters_default
    restart_builder.settings = Dict(dict=settings_dict)
    calc = submit(restart_builder)

    return calc.uuid


def projwfcOriginalSubmit(uuid, codename, metadata, add_parameters={}, del_parameters={}, cluster_options={}):
    """

    :code:`projwfcOriginalSubmit` can submit a projwfc simulation to get the PDOS. It must follow a nscf simulation

    :param uuid: (mandatory) The uuid of previous calculation. We will start our calculation from there. Because uuid
                 is the unique identification number for each CalcJobNode
    :type uuid: python string object

    :param codename: (mandatory) Represent the code for pw.x that you want to use. If you want to use the same as
                     previous calculation, then you need to use Str('')
    :type codename: python string object

    :param add_parameters: (optional, default = {}) The desired parameters that you want to state, it can be incomplete,
                           because inside the function there is a default setting for parameters which can be used in
                           most cases, but if you have specific need, you can put that in parameters, the format is
                           similar as pw.x input file.

                           e.g. :code:`{'PROJWFC':{}}`
    :type add_parameters: python dictionary object

    :param del_parameters: (optional, default = {}) The tags that we would like to delete, for example if we do not
                           want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure
                           as add_parameters.

                           e.g. :code:`{'PROJWFC': [key1, key2, key3]}`
    :type del_parameters: python dictionary object

    :param metadata: (optional, default = {}) The dictionary that contains information about metadata. For example:
                     label and description. label and description are mendatory.
    :type metadata: python dictionary object

    :param cluster_options: (optional, default = {}) The detailed option for the cluster. Different cluster may have
                            different settings. Only the following 3 keys can have effects: (1) resources (2)
                            account (3) queue_name
    :type cluster_options: python dictionary object

    :returns: uuid of the CalcJobNode object of the newest calculation.

    """

    node = load_node(uuid=uuid)

    # check whether it is nscf simulation
    if node.inputs.parameters.get_dict()['CONTROL']['calculation'] != 'nscf':
        return ValueError("You need to provide a nscf simulation with higher k-points.")

    computer = codename.split('@')[1]
    code = Code.get_from_string(codename)
    projwfc_builder = code.get_builder()

    # parameters
    projwfc_parameter = Dict(dict=projwfcParameter)

    # add parameters in add_parameters
    for key, value in add_parameters.items():
        for key2, value2 in value.items():
            projwfc_parameter[key][key2] = value2

    # delete parameters in del_parameters
    for key, value in del_parameters.items():
        tmp = projwfc_parameter[key]
        for key2 in value:
            if key2 in tmp.keys():
                tmp.pop(key2)

    # set default options for slurm
    projwfc_builder.metadata.options['resources'] = slurm_options[computer]['projwfc']['resources']  # in here machine =
    # node
    projwfc_builder.metadata.options['max_wallclock_seconds'] = slurm_options[computer]['projwfc'][
        'max_wallclock_seconds']  # in here machine = node
    projwfc_builder.metadata.options['account'] = slurm_options[computer]['projwfc']['account']  # in here machine = node
    projwfc_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['projwfc']['scheduler_stderr']
    projwfc_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['projwfc']['scheduler_stderr']
    projwfc_builder.metadata.options['queue_name'] = slurm_options[computer]['projwfc']['queue_name']

    # reset cluster_options:
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            projwfc_builder.metadata.options['resources'] = cluster_options['resources']
        if 'account' in cluster_options.keys():
            projwfc_builder.metadata.options['account'] = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            projwfc_builder.metadata.options['queue_name'] = cluster_options['queue_name']

    projwfc_builder.parameters = projwfc_parameter
    projwfc_builder.parent_folder = node.outputs.remote_folder
    projwfc_builder.metadata.label = metadata['label']
    projwfc_builder.metadata.description = metadata['description']

    calc = submit(projwfc_builder)

    return calc.uuid


def phOriginalSubmit(uuid, codename, natlist, qpoints=[[0.0, 0.0, 0.0]], add_parameters={}, del_parameters={},
                     metadata={}, cluster_options={}):
    """

    :code:`phOriginalSubmit` can submit a ph.x simulation to get the PDOS. It must follow a scf simulation.

    :param uuid: (mandatory) The uuid of previous calculation. We will start our calculation from there. Because uuid
                 is the unique identification number for each CalcJobNode
    :type uuid: python string object

    :param codename: (mandatory) Represent the code for pw.x that you want to use. If you want to use the same as
                     previous calculation, then you need to use Str('')
    :type codename: python string object

    :param natlist: (mandatory) Assign the atoms which we want to do the vibrational frequency analysis.
    :type natlist: python list object

    :param qpoints: (optional, default = [[0.0, 0.0, 0.0]] It is like k-points, but useful when calculating
                    vibrational frequencies.
    :type qpoints: python list object

    :param add_parameters: (optional, default = {}) The desired parameters that you want to state, it can be incomplete,
                           because inside the function there is a default setting for parameters which can be used in
                           most cases, but if you have specific need, you can put that in parameters, the format is
                           similar as pw.x input file.

                           e.g. :code:`{'PROJWFC':{}}`
    :type add_parameters: python dictionary object

    :param del_parameters: (optional, default = {}) The tags that we would like to delete, for example if we do not
                           want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure
                           as add_parameters.

                           e.g. :code:`{'PROJWFC': [key1, key2, key3]}`
    :type del_parameters: python dictionary object

    :param metadata: (optional, default = {}) The dictionary that contains information about metadata. For example:
                     label and description. label and description are mendatory.
    :type metadata: python dictionary object

    :param cluster_options: (optional, default = {}) The detailed option for the cluster. Different cluster may have
                            different settings. Only the following 3 keys can have effects: (1) resources (2)
                            account (3) queue_name
    :type cluster_options: python dictionary object

    :returns: uuid of the CalcJobNode object of the newest calculation.


    """

    node = load_node(uuid=uuid)

    # check whether it is nscf simulation
    if node.inputs.parameters.get_dict()['CONTROL']['calculation'] != 'nscf':
        return ValueError("You need to provide a nscf simulation with higher k-points.")

    computer = codename.split('@')[1]
    code = Code.get_from_string(codename)
    ph_builder = code.get_builder()

    # parameters
    ph_parameter = Dict(dict=phParameter)

    # add parameters in add_parameters
    for key, value in add_parameters.items():
        for key2, value2 in value.items():
            ph_parameter[key][key2] = value2

    # delete parameters in del_parameters
    for key, value in del_parameters.items():
        tmp = ph_parameter[key]
        for key2 in value:
            if key2 in tmp.keys():
                tmp.pop(key2)

    # set kpoints
    qpts = KpointsData()
    if len(qpoints) == 1:
        qpts.set_kpoints_mesh(mesh=qpoints[0])
    else:
        qpts.set_kpoints_mesh(mesh=qpoints[0], offset=qpoints[1])

    # set default options for slurm
    # set first, then modify
    ph_builder.metadata.options['resources'] = slurm_options[computer]['ph']['resources']  # in here machine =
    # node
    ph_builder.metadata.options['max_wallclock_seconds'] = slurm_options[computer]['projwfc'][
        'max_wallclock_seconds']  # in here machine = node
    ph_builder.metadata.options['account'] = slurm_options[computer]['ph']['account']  # in here machine = node
    ph_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['ph']['scheduler_stderr']
    ph_builder.metadata.options['scheduler_stderr'] = slurm_options[computer]['ph']['scheduler_stderr']
    ph_builder.metadata.options['queue_name'] = slurm_options[computer]['ph']['queue_name']

    # reset cluster_options:
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            ph_builder.metadata.options['resources'] = cluster_options['resources']
        if 'account' in cluster_options.keys():
            ph_builder.metadata.options['account'] = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            ph_builder.metadata.options['queue_name'] = cluster_options['queue_name']

    ph_builder.parameters = Dict(dict=ph_parameter)
    ph_builder.parent_folder = node.outputs.remote_folder
    ph_builder.metadata.label = metadata['label']
    ph_builder.metadata.description = metadata['description']
    ph_builder.qpoints = qpts

    calc = submit(ph_builder)

    return calc.uuid


def ppOriginalSubmit(results, uuid, codename, add_parameters, del_parameters, metadata):
    pass


def nebOriginalSubmit(results, uuid, codename, add_parameters, del_parameters, metadata):
    pass
