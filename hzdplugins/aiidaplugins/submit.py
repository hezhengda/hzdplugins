# Since the goal of Aiida is to "preservation of provenance", so I should write this plugin with the same goal as
# Aiida, that's why all the relevant information will be provided in aiida.orm types. (I can convert them in the
# code, it is very straightforward)

from copy import deepcopy

from aiida.engine.launch import submit
from aiida.orm import KpointsData, Dict, Code
from aiida.orm import load_node
from aiida.orm.nodes.data.upf import get_pseudos_from_structure

from hzdplugins.aiidaplugins.constants import slurm_options


def qePwOriginalSubmit(results, codename, structure, kpoints, pseudo_family, pseudo_dict, metadata, add_parameters,
                       del_parameters, cluster_options, settings_dict):
    """

    :code:`qePwOriginalSubmit` will submit an original computational task to the desired computer by using certain code.

    :param results: A dictionary that has all the relevant information about the simulation, its key is the uuid of
                    the CalcJobNode
    :type results: python dictionary

    :param codename: A string represents the code for pw.x that you want to use.
    :type codename: python string object

    :param structure: The structure of your system.
    :type structure: aiida.orm.StructureData object

    :param add_parameters: The desired parameters that you want to state, it can be incomplete, because inside the
                           function there is a default setting for parameters which can be used in most cases,
                           but if you have specific need, you can put that in parameters, the format is similar as
                           pw.x input file.

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

    :param del_parameters: The tags that we would like to delete, for example if we do not want to use spin-polarized
                           simulation, then 'nspin' needs to be deleted. Same structure as add_parameters.

                           e.g. :code:`{'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}`
    :type del_parameters: python dictionary object

    :param kpoints: The kpoints that you want to use, if the kpoints has only 1 list, then it is the kpoint mesh,
                    but if two lists are detected, then the first will be k-point mesh, the second one will be the
                    origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]
    :type kpoints: python list object

    :param pseudo_family: The pseudopotential family that you want to use. Make sure that you already have that
                          configured, otherwise an error will occur.
    :type pseudo_family: python string object.

    :param pseudo_dict: which contains the pseudopotential files that we want to use in the simulation.
                        In here it is very important to note that the path of the pseudopotential file has to be in
                        the absolute path.

                        e.g.

                        .. code-block:: python

                            pseudo_dict = {
                                'Fe': UpfData(absolute_path),
                                'Fe3': UpfData(absolute_path)
                            }

    :type pseudo_dict: python dictionary object.

    :param cluster_options: The detailed option for the cluster. Different cluster may have different settings. Only
                            the following 3 keys can have effects: (1) resources (2) account (3) queue_name
    :type cluster_options: python dictionary object

    :param metadata: The dictionary that contains information about metadata. For example: label and
                     description. label and description are mendatory.

                     e.g. :code:`{'label':{}, 'description':{}}`
    :type metadata: python dictionary object

    :param settings_dict: which contains the additional information for the pw.x calculation. e.g. Fixed atom,
                          retrieving more files, parser options, etc. And the command-line options.
    :type settings_dict: python dictionary object

    :returns: - **results** (`python dictionary object`): a modified results dictionary with the latest submitted job

              - **uuid** (`python string object`): the uuid of that CalcJob

    """

    results_tmp = deepcopy(results)  # first we need to create a copy for our simulation

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
    parameters_default = Dict(dict={
        'CONTROL': {
            'calculation': 'vc-relax',
            'max_seconds': 86000,
            'restart_mode': 'from_scratch',
            'wf_collect': True,
            'nstep': 50000,
            'tstress': True,
            'tprnfor': True,
            'etot_conv_thr': 1e-06,
            'forc_conv_thr': 0.001,
            'disk_io': 'low',
            'verbosity': 'low'
        },
        'SYSTEM': {
            'ibrav': 0,
            'nosym': False,
            'ecutwfc': 80.0,
            'ecutrho': 640.0,
            'occupations': 'smearing',
            'degauss': 0.002,
            'smearing': 'gaussian',
            'input_dft': 'PBE'
        },
        'ELECTRONS': {
            'electron_maxstep': 1000,
            'conv_thr': 1.0e-6,
            'diagonalization': 'david',
            'mixing_mode': 'plain',
            'mixing_beta': 0.3,
            'mixing_ndim': 10
        }
    })

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
    pw_builder.metadata.options.resources = slurm_options[computer]['qe']['resources']  # in here machine = node
    pw_builder.metadata.options.max_wallclock_seconds = slurm_options[computer]['qe'][
        'max_wallclock_seconds']  # in here machine = node
    pw_builder.metadata.options.account = slurm_options[computer]['qe']['account']  # in here machine = node
    pw_builder.metadata.options.scheduler_stderr = slurm_options[computer]['qe']['scheduler_stderr']
    pw_builder.metadata.options.scheduler_stdout = slurm_options[computer]['qe']['scheduler_stderr']
    pw_builder.metadata.options.queue_name = slurm_options[computer]['qe']['queue_name']

    # revised by cluster_options
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            pw_builder.metadata.options.resources = cluster_options['resources']
        if 'account' in cluster_options.keys():
            pw_builder.metadata.options.resources = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            pw_builder.metadata.options.resources = cluster_options['queue_name']

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

    # results
    results_tmp[str(calc.uuid)] = {}
    results_tmp[str(calc.uuid)]['system'] = metadata['label']
    results_tmp[str(calc.uuid)]['comp_type'] = parameters_default['CONTROL']['calculation']
    results_tmp[str(calc.uuid)]['E/eV'] = None
    results_tmp[str(calc.uuid)]['remove_remote_folder'] = False
    results_tmp[str(calc.uuid)]['cluster'] = codename.split('@')[
        1]  # becase all the codename have same structure "code@computer"
    results_tmp[str(calc.uuid)]['xc functional'] = parameters_default['SYSTEM']['input_dft']
    results_tmp[str(calc.uuid)]['exit_status'] = None
    results_tmp[str(calc.uuid)]['is_finished'] = None
    results_tmp[str(calc.uuid)]['is_finished_ok'] = None
    results_tmp[str(calc.uuid)]['previous_calc'] = 0  # 0 represent original
    results_tmp[str(calc.uuid)]['son_calc'] = None  # currently no son_calc node
    # results_tmp[ads_str(calc.pk)]['description'] = metadata['description']

    return results_tmp, calc.uuid


def qePwContinueSubmit(results, uuid, pseudo_family, pseudo_dict, codename, parent_folder, add_parameters,
                       del_parameters, kpoints, cluster_options, metadata, settings_dict):
    """

    `qePwContinueSubmit` will continue a simulation with similar or modified input parameters. All the parameters are
    listed in the kwargs.

    :param results: A dictionary that has all the relevant information about the simulation, its key is the uuid of the
                    CalcJobNode
    :type results: python dictionary object

    :param uuid: The uuid of previous calculation. We will start our calculation from there. Because uuid is the
                 unique identification number for each CalcJobNode

                    **Notice**: The uuid must be in the results dictionary, if not the program will shout KeyError.
                    And if you are testing, you could use assignValue to quickly create a dictionary that contains
                    the uuid that you want to continue.
    :type uuid: python string object

    :param pseudo_family: The pseudopotential family that you want to use. Make sure that you already have that
                          configured, otherwise an error will occur. This is mendatory.
    :type pseudo_family: python string object

    :param pseudo_dict: Which contains the pseudopotential files that we want to use in the simulation.
    :type pseudo_dict: python dictionary object

    :param codename: Represent the code for pw.x that you want to use. If you want to use the same as
                     previous calculation, then you need to use Str('')
    :type codename: python string object

    :param parent_folder: If parent_folder is True, then the calculation will start with the output files from
                          previous calculations.
    :type parent_folder: python boolean object

    :param add_parameters: The desired parameters that you want to state, it can be incomplete,
                           because inside the function there is a default setting for parameters which can be used in
                           most cases, but if you have specific need, you can put that in parameters, the format is
                           similar as pw.x input file.

                           If you want to assign DFT+U and spin-polarization, you need to specify it on your own.

                           e.g. :code:`{'CONTROL':{}, 'SYSTEM':{}}`

                           **Notice**: more options in qePwOriginalSubmit function.
    :type add_parameters: python dictionary object

    :param del_parameters: The tags that we would like to delete, for example if we do not want to use spin-polarized
                           simulation, then 'nspin' needs to be deleted. Same structure as add_parameters.

                           e.g. :code:`{'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}`
    :type del_parameters: python dictionary object

    :param kpoints: optional, if you want to keep the k-points for previous calculation, just use an empty list
                    :code:`[]`. The kpoints that you want to use, if the kpoints has only 1 list, then it is the
                    kpoint  mesh, but if two lists are detected, then the first will be k-point mesh, the second
                    one will be the origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]
    :type kpoints: python list object

    :param cluster_options: The detailed option for the cluster. Different cluster may have different
                            settings. Only the following 3 keys can have effects: (1) resources (2) account (3)
                            queue_name
    :type cluster_options: python dictionary object

    :param metadata: The dictionary that contains information about metadata. For example: label and description.
                     label and description are mendatory.
    :type metadata: python dictionary object

    :param settings_dict: which contains the additional information for the pw.x calculation. e.g. Fixed atom,
                          retrieving more files, parser options, etc. And the command-line options.
    :type settings_dict: python dictionary object

    :returns: - results (`python dictionary object`): a modified results dictionary with the latest submitted job

              - uuid (`python dictionary object`): the uuid of that CalcJob

    """

    results_tmp = deepcopy(results)

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

    if calc_type == 'relax' or calc_type == 'vc-relax':
        structure = node.outputs.output_structure
    elif calc_type == 'scf' or calc_type == 'nscf':
        structure = node.inputs.structure

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

    # reset cluster_options:
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            restart_builder.metadata.options.resources = cluster_options['resources']
        if 'account' in cluster_options.keys():
            restart_builder.metadata.options.account = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            restart_builder.metadata.options.queue_name = cluster_options['queue_name']
    else:
        restart_builder.metadata.options.update({
            'resources': slurm_options[computer]['qe']['resources'],
            'max_wallclock_seconds': slurm_options[computer]['qe']['max_wallclock_seconds'],
            'account': slurm_options[computer]['qe']['account'],
            'scheduler_stderr': slurm_options[computer]['qe']['scheduler_stderr'],
            'scheduler_stdout': slurm_options[computer]['qe']['scheduler_stdout'],
            'queue_name': slurm_options[computer]['qe']['queue_name']
        })

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

    # results
    results_tmp[str(calc.uuid)] = {}
    results_tmp[str(calc.uuid)]['system'] = restart_builder.metadata.label
    results_tmp[str(calc.uuid)]['comp_type'] = parameters_default['CONTROL']['calculation']
    results_tmp[str(calc.uuid)]['E/eV'] = None
    results_tmp[str(calc.uuid)]['remove_remote_folder'] = False
    results_tmp[str(calc.uuid)]['cluster'] = computer
    results_tmp[str(calc.uuid)]['xc functional'] = parameters_default['SYSTEM']['input_dft']
    results_tmp[str(calc.uuid)]['exit_status'] = None
    results_tmp[str(calc.uuid)]['is_finished'] = None
    results_tmp[str(calc.uuid)]['is_finished_ok'] = None
    results_tmp[str(calc.uuid)]['previous_calc'] = uuid  # pk is the previous calculation
    results_tmp[str(calc.uuid)]['son_calc'] = None  # right now there is no son_node

    # change son_calc with previous simulation
    results_tmp[uuid]['son_calc'] = calc.uuid

    return results_tmp, calc.uuid


def projwfcOriginalSubmit(results, uuid, codename, add_parameters, del_parameters, metadata, cluster_options):
    """

    :code:`projwfcOriginalSubmit` can submit a projwfc simulation to get the PDOS. It must follow a nscf simulation


    :param results: A dictionary that has all the relevant information about the simulation, its key is the uuid of the
                    CalcJobNode
    :type results: python dictionary object

    :param uuid: The uuid of previous calculation. We will start our calculation from there. Because uuid is the
                 unique identification number for each CalcJobNode
    :type uuid: python string object

    :param codename: Represent the code for pw.x that you want to use. If you want to use the same as
                     previous calculation, then you need to use Str('')
    :type codename: python string object

    :param add_parameters: The desired parameters that you want to state, it can be incomplete,
                           because inside the function there is a default setting for parameters which can be used in
                           most cases, but if you have specific need, you can put that in parameters, the format is
                           similar as pw.x input file.

                           e.g. :code:`{'PROJWFC':{}}`
    :type add_parameters: python dictionary object

    :param del_parameters: The tags that we would like to delete, for example if we do not want to use
                           spin-polarized simulation, then 'nspin' needs to be deleted. Same structure as
                           add_parameters.

                           e.g. :code:`{'PROJWFC': [key1, key2, key3]}`
    :type del_parameters: python dictionary object

    :param metadata: The dictionary that contains information about metadata. For example: label and
                     description. label and description are mendatory.
    :type metadata: python dictionary object

    :param cluster_options: The detailed option for the cluster. Different cluster may have different
                            settings. Only the following 3 keys can have effects: (1) resources (2) account (3)
                            queue_name
    :type cluster_options: python dictionary object

    :returns: - **results** (`python dictionary object`): results contains the latest calculation information.

              - **uuid** (`python string object`): uuid of the latest calculation

    """

    results_tmp = deepcopy(results)

    node = load_node(uuid=uuid)

    # check whether it is nscf simulation
    if node.inputs.parameters.get_dict()['CONTROL']['calculation'] != 'nscf':
        return ValueError("You need to provide a nscf simulation with higher k-points.")

    computer = codename.split('@')[1]
    code = Code.get_from_string(codename)
    projwfc_builder = code.get_builder()

    # parameters
    projwfc_parameter = {
        'PROJWFC': {
            'DeltaE': 0.01,
            'ngauss': 0,
            'degauss': 0.015,
            'Emin': -40,
            'Emax': 40
        }
    }

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

    # reset cluster_options:
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            projwfc_builder.metadata.options.resources = cluster_options['resources']
        if 'account' in cluster_options.keys():
            projwfc_builder.metadata.options.account = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            projwfc_builder.metadata.options.queue_name = cluster_options['queue_name']
    else:
        projwfc_builder.metadata.options.update({
            'resources': slurm_options[computer]['projwfc']['resources'],
            'max_wallclock_seconds': slurm_options[computer]['projwfc']['max_wallclock_seconds'],
            'account': slurm_options[computer]['projwfc']['account'],
            'scheduler_stderr': slurm_options[computer]['projwfc']['scheduler_stderr'],
            'scheduler_stdout': slurm_options[computer]['projwfc']['scheduler_stdout'],
            'queue_name': slurm_options[computer]['projwfc']['queue_name']
        })

    projwfc_builder.parameters = Dict(dict=projwfc_parameter)
    projwfc_builder.parent_folder = node.outputs.remote_folder
    projwfc_builder.metadata.label = metadata['label']
    projwfc_builder.metadata.description = metadata['description']

    calc = submit(projwfc_builder)

    # results
    results_tmp[str(calc.uuid)] = {}
    results_tmp[str(calc.uuid)]['system'] = projwfc_builder.metadata.label
    results_tmp[str(calc.uuid)]['comp_type'] = 'projwfc'
    results_tmp[str(calc.uuid)]['E/eV'] = None
    results_tmp[str(calc.uuid)]['remove_remote_folder'] = False
    results_tmp[str(calc.uuid)]['cluster'] = computer
    results_tmp[str(calc.uuid)]['xc functional'] = node.inputs.parameters.get_dict()['SYSTEM']['input_dft']
    results_tmp[str(calc.uuid)]['exit_status'] = None
    results_tmp[str(calc.uuid)]['is_finished'] = None
    results_tmp[str(calc.uuid)]['is_finished_ok'] = None
    results_tmp[str(calc.uuid)]['previous_calc'] = uuid  # pk is the previous calculation
    results_tmp[str(calc.uuid)]['son_calc'] = None  # right now there is no son_node

    # change son_calc with previous simulation
    results_tmp[uuid]['son_calc'] = calc.uuid

    return results_tmp, calc.uuid


def phOriginalSubmit(results, uuid, codename, qpoints, add_parameters, del_parameters, metadata, cluster_options):
    """

    :code:`phOriginalSubmit` can submit a ph.x simulation to get the PDOS. It must follow a scf simulation.

    :param results: A dictionary that has all the relevant information about the simulation, its key is the uuid of the
                    CalcJobNode
    :type results: python dictionary object

    :param uuid: The uuid of previous calculation. We will start our calculation from there. Because uuid is the
                 unique identification number for each CalcJobNode
    :type uuid: python string object

    :param codename: Represent the code for pw.x that you want to use. If you want to use the same as
                     previous calculation, then you need to use :code:`Str('')`
    :type codename: python string object

    :param qpoints: A list that may contain 1 / 2 lists. To show the qpoints we want to use. e.g. qpoints = [[0.0, 0.0,
                    0.0]] or qpoints = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], first one is mesh, second is offset.
    :type qpoints: python list object

    :param add_parameters: The desired parameters that you want to state, it can be incomplete, because inside the
                           function there is a default setting for parameters which can be used in most cases,
                           but if you have specific need, you can put that in parameters, the format is similar as
                           pw.x input file.

                           e.g. {'INPUTPH':{}}
    :type add_parameters: python dictionary object

    :param del_parameters: The tags that we would like to delete, for example if we do not want to use
                           spin-polarized simulation, then 'nspin' needs to be deleted. Same structure as
                           add_parameters.

                           e.g. {'INPUTPH': [key1, key2, key3]}
    :type del_parameters: python dictionary object

    :param metadata: The dictionary that contains information about metadata. For example: label and description.
                     label and description are mendatory.
    :type metadata: python dictionary object

    :param cluster_options: The detailed option for the cluster. Different cluster may have different
                            settings. Only the following 3 keys can have effects: (1) resources (2) account (3)
                            queue_name
    :type cluster_options: python dictionary object

    :returns: - **results** (`python dictionary object`): results contains the latest calculation information.

              - **uuid** (`python string object`): uuid of the latest calculation

    """

    results_tmp = deepcopy(results)

    node = load_node(uuid=uuid)

    # check whether it is nscf simulation
    if node.inputs.parameters.get_dict()['CONTROL']['calculation'] != 'nscf':
        return ValueError("You need to provide a nscf simulation with higher k-points.")

    computer = codename.split('@')[1]
    code = Code.get_from_string(codename)
    ph_builder = code.get_builder()

    # parameters
    ph_parameter = {
        'INPUTPH': {
            'title_line': 'This is a ph.x calculation after the CalcNode uuid = {}'.format(uuid),
            'max_seconds': 86000,
            'tr2_ph': 1.0e-8,
            'ldisp': False,
            'epsil': False,
            'trans': True
        }
    }

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

    # reset cluster_options:
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            ph_builder.metadata.options.resources = cluster_options['resources']
        if 'account' in cluster_options.keys():
            ph_builder.metadata.options.account = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            ph_builder.metadata.options.queue_name = cluster_options['queue_name']
    else:
        ph_builder.metadata.options.update({
            'resources': slurm_options[computer]['projwfc']['resources'],
            'max_wallclock_seconds': slurm_options[computer]['projwfc']['max_wallclock_seconds'],
            'account': slurm_options[computer]['projwfc']['account'],
            'scheduler_stderr': slurm_options[computer]['projwfc']['scheduler_stderr'],
            'scheduler_stdout': slurm_options[computer]['projwfc']['scheduler_stdout'],
            'queue_name': slurm_options[computer]['projwfc']['queue_name']
        })

    ph_builder.parameters = Dict(dict=ph_parameter)
    ph_builder.parent_folder = node.outputs.remote_folder
    ph_builder.metadata.label = metadata['label']
    ph_builder.metadata.description = metadata['description']
    ph_builder.qpoints = qpts

    calc = submit(ph_builder)

    # results
    results_tmp[str(calc.uuid)] = {}
    results_tmp[str(calc.uuid)]['system'] = ph_builder.metadata.label
    results_tmp[str(calc.uuid)]['comp_type'] = 'ph'
    results_tmp[str(calc.uuid)]['E/eV'] = None
    results_tmp[str(calc.uuid)]['remove_remote_folder'] = False
    results_tmp[str(calc.uuid)]['cluster'] = computer
    results_tmp[str(calc.uuid)]['xc functional'] = node.inputs.parameters.get_dict()['SYSTEM']['input_dft']
    results_tmp[str(calc.uuid)]['exit_status'] = None
    results_tmp[str(calc.uuid)]['is_finished'] = None
    results_tmp[str(calc.uuid)]['is_finished_ok'] = None
    results_tmp[str(calc.uuid)]['previous_calc'] = uuid  # pk is the previous calculation
    results_tmp[str(calc.uuid)]['son_calc'] = None  # right now there is no son_node

    # change son_calc with previous simulation
    results_tmp[uuid]['son_calc'] = calc.uuid

    return results_tmp, calc.uuid


def ppOriginalSubmit(results, uuid, codename, add_parameters, del_parameters, metadata):
    pass


def nebOriginalSubmit(results, uuid, codename, add_parameters, del_parameters, metadata):
    pass
