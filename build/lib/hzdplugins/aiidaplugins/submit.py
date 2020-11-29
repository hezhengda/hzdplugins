from aiida.orm import StructureData, KpointsData, Dict
from copy import deepcopy
from aiida.orm import load_node
from aiida.engine.launch import submit
from hzdplugins.aiidaplugins.constants import results_keys_set, slurm_options

def qePwOriginalSubmit(results, codename, structure, add_parameters={}, del_parameters={}, kpoints, pseudo_family, cluster_options={}, metadata):

    """

    `qePwOriginalSubmit` will submit an original computational task to the desired computer by using certain code.

    Parameters:

    results:
        A dictionary that has all the relevant information about the

    codename:
        A string. Represent the code for pw.x that you want to use

    structure:
        A StructureData object. The structure of your system. It needs to be StructureData type.

    add_parameters:
        A dictionary. The desired parameters that you want to state, it can be incomplete, because inside the function there is a default setting for parameters which can be used in most cases, but if you have specific need, you can put that in parameters, the format is similar as pw.x input file.

        If you want to assign DFT+U and spin-polarization, you need to specify it on your own.

        e.g. {'CONTROL':{}, 'SYSTEM':{}}

    del_parameters:
        A dictionary. The tags that we would like to delete, for example if we do not want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure as add_parameters.

        e.g. {'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}

    kpoints:
        A list of lists. The kpoints that you want to use, if the kpoints has only 1 list, then it is the kpoint mesh, but if two lists are detected, then the first will be k-point mesh, the second one will be the origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]

    pseudo_family:
        A string. The pseudopotential family that you want to use. Make sure that you already have that configured, otherwise an error will occur.

    cluster_options:
        A dictionary. The detailed option for the cluster. Different cluster may have different settings. Only the following 3 keys can have effects: (1) resources (2) account (3) queue_name

    metadata:
        A dictionary. The dictionary that contains information about metadata. For example: label and description. label and description are mendatory.

    Return:
        results: a modified results dictionary with the latest submitted job
        pk: the id of that CalcJob

    """

    results_tmp = deepcopy(results) # first we need to create a copy for our simulation

    code = Code.get_from_string(codename)
    computer = codename.split('@')[1] # get the name of the cluster
    pw_builder = code.get_builder()

    # pseudopotential
    pw_builder.pseudos = get_pseudos_from_structure(structure, family_name=pseudo_family)

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
            'mixing_mode': 'local-TF',
            'mixing_beta': 0.2,
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
    pw_builder.metadata.options.resources = slurm_options[computer]['resources'] # in here machine = node
    pw_builder.metadata.options.max_wallclock_seconds = slurm_options[computer]['max_wallclock_seconds'] # in here machine = node
    pw_builder.metadata.options.account = slurm_options[computer]['account'] # in here machine = node
    pw_builder.metadata.options.scheduler_stderr = slurm_options[computer]['scheduler_stderr']
    pw_builder.metadata.options.scheduler_stdout = slurm_options[computer]['scheduler_stderr']
    pw_builder.metadata.options.queue_name = slurm_options[computer]['queue_name']

    # revised by cluster_options
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            pw_builder.metadata.options.resources = cluster_options['resources']
        if 'account' in cluster_options.keys():
            pw_builder.metadata.options.resources = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            pw_builder.metadata.options.resources = cluster_options['queue_name']

    # launch the simulation
    pw_builder.structure = structure
    pw_builder.kpoints = kpts
    pw_builder.parameters = parameters_default
    pw_builder.settings = Dict(dict={'CMDLINE': ['-npools', '4']})
    calc = submit(pw_builder)

    # results
    results_tmp[str(calc.pk)] = {}
    results_tmp[str(calc.pk)]['uuid'] = calc.uuid
    results_tmp[str(calc.pk)]['system'] = metadata['label']
    results_tmp[str(calc.pk)]['comp_type'] = parameters_default['CONTROL']['calculation']
    results_tmp[str(calc.pk)]['E/eV'] = None
    results_tmp[str(calc.pk)]['remove_remote_folder'] = False
    results_tmp[str(calc.pk)]['cluster'] = codename.split('@')[1] # becase all the codename have same structure "code@computer"
    results_tmp[str(calc.pk)]['xc functional'] = parameters_default['SYSTEM']['input_dft']
    results_tmp[str(calc.pk)]['exit_status'] = None
    results_tmp[str(calc.pk)]['is_finished'] = None
    results_tmp[str(calc.pk)]['is_finished_ok'] = None
    results_tmp[str(calc.pk)]['previous_calc'] = 0 # 0 represent original
    # results_tmp[str(calc.pk)]['description'] = metadata['description']

    return results_tmp, calc.pk

def qePwContinueSubmit(results, pk, codename='', add_parameters={}, del_parameters={}, kpoints=[], pseudo_family='', cluster_options={}, metadata={}):

    """

    `qePwContinueSubmit` will continue a simulation with similar or modified input parameters. All the parameters are listed in the kwargs

    Parameters:

    pk:
        The id of previous calculation. We will start our calculation from there.

    codename:
        A string. Represent the code for pw.x that you want to use.

    add_parameters:
        A dictionary. The desired parameters that you want to state, it can be incomplete, because inside the function there is a default setting for parameters which can be used in most cases, but if you have specific need, you can put that in parameters, the format is similar as pw.x input file.

        If you want to assign DFT+U and spin-polarization, you need to specify it on your own.

        e.g. {'CONTROL':{}, 'SYSTEM':{}}

    del_parameters:
        A dictionary. The tags that we would like to delete, for example if we do not want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure as add_parameters.

        e.g. {'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}

    kpoints:
        A list of lists. The kpoints that you want to use, if the kpoints has only 1 list, then it is the kpoint mesh, but if two lists are detected, then the first will be k-point mesh, the second one will be the origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]

    pseudo_family:
        A string. The pseudopotential family that you want to use. Make sure that you already have that configured, otherwise an error will occur.

    cluster_options:
        A dictionary. The detailed option for the cluster. Different cluster may have different settings. Only the following 3 keys can have effects: (1) resources (2) account (3) queue_name

    metadata:
        A dictionary. The dictionary that contains information about metadata. For example: label and description. label and description are mendatory.

    Return:
        results: a modified results dictionary with the latest submitted job
        pk: the id of that CalcJob

    Return:
        results: a modified results dictionary with the latest submitted job
        pk: the id of that CalcJob

    """

    results_tmp = deepcopy(results)

    node = load_node(pk)

    if len(codename) == 0: # not going to change cluster
        computer = node.computer.label
        restart_builder = node.get_builder_restart() # get the restart_builder
    else:
        copmuter = codename.split('@')[1]
        code = Code.get_from_string(codename)
        restart_builder = code.get_builder()

    parameters_tmp = deepcopy(node.inputs.parameters)
    structure = node.outputs.output_structure

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

    # reset the kpoints
    if len(kpoints) > 0:
        kpts = KpointsData()
        if len(kpoints) == 1:
            kpts.set_kpoints_mesh(mesh=kpoints[0])
        else:
            kpts.set_kpoints_mesh(mesh=kpoints[0], offset=kpoints[1])
    else:
        kpts = node.inputs.kpoints

    # reset the pseudo_family
    if pseudo_family != '':
        restart_builder.pseudos = get_pseudos_from_structure(structure, family_name=pseudo_family)

    # reset cluster_options:
    if len(cluster_options) > 0:
        if 'resources' in cluster_options.keys():
            restart_builder.metadata.options.resources = cluster_options['resources']
        if 'account' in cluster_options.keys():
            restart_builder.metadata.options.account = cluster_options['account']
        if 'queue_name' in cluster_options.keys():
            restart_builder.metadata.options.queue_name = cluster_options['queue_name']

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
        restart_builder.metadata.label = node.label = node.label
        restart_builder.metadata.description = node.description = node.description

    # submit the calculation
    restart_builder.structure = structure
    restart_builder.kpoints = kpts
    restart_builder.parameters = parameters_default
    restart_builder.settings = Dict(dict={'CMDLINE': ['-npools', '4']})
    calc = submit(restart_builder)

    # results
    results_tmp[str(calc.pk)] = {}
    results_tmp[str(calc.pk)]['uuid'] = calc.uuid
    results_tmp[str(calc.pk)]['system'] = metadata['label']
    results_tmp[str(calc.pk)]['comp_type'] = parameters_default['CONTROL']['calculation']
    results_tmp[str(calc.pk)]['E/eV'] = None
    results_tmp[str(calc.pk)]['remove_remote_folder'] = False
    results_tmp[str(calc.pk)]['cluster'] = computer
    results_tmp[str(calc.pk)]['xc functional'] = parameters_default['SYSTEM']['input_dft']
    results_tmp[str(calc.pk)]['exit_status'] = None
    results_tmp[str(calc.pk)]['is_finished'] = None
    results_tmp[str(calc.pk)]['is_finished_ok'] = None
    results_tmp[str(calc.pk)]['previous_calc'] = pk # pk is the previous calculation

    return results_tmp, calc.pk
