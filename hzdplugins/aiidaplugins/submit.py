# Since the goal of Aiida is to "preservation of provenance", so I should write this plugin with the same goal as Aiida, that's why all the relevant information will be provided in aiida.orm types. (I can convert them in the code, it is very straightforward)

from aiida.orm import StructureData, KpointsData, Dict, Code
from copy import deepcopy
from aiida.orm import load_node
from aiida.engine import calcfunction
from aiida.engine.launch import submit
from hzdplugins.aiidaplugins.constants import results_keys_set, slurm_options
from aiida.orm.nodes.data.upf import get_pseudos_from_structure

@calcfunction
def qePwOriginalSubmit(results, codename, structure, kpoints, pseudo_family, metadata, add_parameters, del_parameters, cluster_options):

    """

    `qePwOriginalSubmit` will submit an original computational task to the desired computer by using certain code.

    **Notice**: Every input parameters should be in `aiida.orm` types.

    Parameters:

    results:
        A dictionary that has all the relevant information about the simulation, its key is the uuid of the CalcJobNode

    codename:
        An aiida.orm.Str object. A string represents the code for pw.x that you want to use

    structure:
        An aiida.orm.StructureData object. The structure of your system. It needs to be StructureData type.

    add_parameters:
        An aiida.orm.Dict object. A dictionary. The desired parameters that you want to state, it can be incomplete, because inside the function there is a default setting for parameters which can be used in most cases, but if you have specific need, you can put that in parameters, the format is similar as pw.x input file.

        If you want to assign DFT+U and spin-polarization, you need to specify it on your own.

        e.g. {'CONTROL':{}, 'SYSTEM':{}}

    del_parameters:
        An aiida.orm.Dict object. A dictionary. The tags that we would like to delete, for example if we do not want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure as add_parameters.

        e.g. {'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}

    kpoints:
        An aiida.orm.List object. A list of lists. The kpoints that you want to use, if the kpoints has only 1 list, then it is the kpoint mesh, but if two lists are detected, then the first will be k-point mesh, the second one will be the origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]

    pseudo_family:
        An aiida.orm.Str object. A string. The pseudopotential family that you want to use. Make sure that you already have that configured, otherwise an error will occur.

    cluster_options:
        An aiida.orm.Dict object. A dictionary. The detailed option for the cluster. Different cluster may have different settings. Only the following 3 keys can have effects: (1) resources (2) account (3) queue_name

    metadata:
        An aiida.orm.Dict object. A dictionary. The dictionary that contains information about metadata. For example: label and description. label and description are mendatory.

    Return:
        results: a modified results dictionary with the latest submitted job
        pk: the id of that CalcJob

    """

    # clean the input from aiida.orm to regular datatypes
    codename = codename.value
    add_parameters = add_parameters.get_dict()
    del_parameters = del_parameters.get_dict()
    kpoints = kpoints.get_list()
    pseudo_family = pseudo_family.value
    metadata = metadata.get_dict()
    cluster_options = cluster_options.get_dict()
    # the end of the cleaning procedure

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
    results_tmp[str(calc.uuid)] = {}
    results_tmp[str(calc.uuid)]['system'] = metadata['label']
    results_tmp[str(calc.uuid)]['comp_type'] = parameters_default['CONTROL']['calculation']
    results_tmp[str(calc.uuid)]['E/eV'] = None
    results_tmp[str(calc.uuid)]['remove_remote_folder'] = False
    results_tmp[str(calc.uuid)]['cluster'] = codename.split('@')[1] # becase all the codename have same structure "code@computer"
    results_tmp[str(calc.uuid)]['xc functional'] = parameters_default['SYSTEM']['input_dft']
    results_tmp[str(calc.uuid)]['exit_status'] = None
    results_tmp[str(calc.uuid)]['is_finished'] = None
    results_tmp[str(calc.uuid)]['is_finished_ok'] = None
    results_tmp[str(calc.uuid)]['previous_calc'] = 0 # 0 represent original
    results_tmp[str(calc.uuid)]['son_calc'] = None # currently no son_calc node
    # results_tmp[str(calc.pk)]['description'] = metadata['description']

    return results_tmp, calc.uuid

@calcfunction
def qePwContinueSubmit(results, uuid, pseudo_family, codename, add_parameters, del_parameters, kpoints, cluster_options, metadata):

    """

    `qePwContinueSubmit` will continue a simulation with similar or modified input parameters. All the parameters are listed in the kwargs.

    **Notice**: Every input parameters should be in `aiida.orm` types.

    Parameters:

    uuid:
        An aiida.orm.Str object. The uuid of previous calculation. We will start our calculation from there. Because uuid is the unique identification number for each CalcJobNode

    pseudo_family:
        An aiida.orm.Str object. A string. The pseudopotential family that you want to use. Make sure that you already have that configured, otherwise an error will occur. This is mendatory.

    codename:
        An aiida.orm.Str object. A string. Represent the code for pw.x that you want to use. If you want to use the same as previous calculation, then you need to use Str('')

    add_parameters:
        An aiida.orm.Dict object. A dictionary. The desired parameters that you want to state, it can be incomplete, because inside the function there is a default setting for parameters which can be used in most cases, but if you have specific need, you can put that in parameters, the format is similar as pw.x input file.

        If you want to assign DFT+U and spin-polarization, you need to specify it on your own.

        e.g. {'CONTROL':{}, 'SYSTEM':{}}

    del_parameters:
        An aiida.orm.Dict object. A dictionary. The tags that we would like to delete, for example if we do not want to use spin-polarized simulation, then 'nspin' needs to be deleted. Same structure as add_parameters.

        e.g. {'CONTROL': [key1, key2, key3], 'SYSTEM': [key1, key2, key3]}

    kpoints:
        An aiida.orm.List object. A list of lists. The kpoints that you want to use, if the kpoints has only 1 list, then it is the kpoint mesh, but if two lists are detected, then the first will be k-point mesh, the second one will be the origin of k-point mesh.e.g. [[3, 3, 1]] or [[3, 3, 1],[0.5, 0.5, 0.5]]

    cluster_options:
        An aiida.orm.Dict object. A dictionary. The detailed option for the cluster. Different cluster may have different settings. Only the following 3 keys can have effects: (1) resources (2) account (3) queue_name

    metadata:
        An aiida.orm.Dict object. A dictionary. The dictionary that contains information about metadata. For example: label and description. label and description are mendatory.

    Return:
        results: a modified results dictionary with the latest submitted job
        uuid: the uuid of that CalcJob

    """

    # clean the input from aiida.orm to regular datatypes
    uuid = uuid.value
    codename = codename.value
    add_parameters = add_parameters.get_dict()
    del_parameters = del_parameters.get_dict()
    kpoints = kpoints.get_list()
    pseudo_family = pseudo_family.value
    metadata = metadata.get_dict()
    cluster_options = cluster_options.get_dict()
    # the end of the cleaning procedure

    results_tmp = deepcopy(results)

    node = load_node(uuid=uuid)

    if len(codename) == 0: # not going to change cluster
        computer = node.computer.label
        restart_builder = node.get_builder_restart() # get the restart_builder
    else:
        computer = codename.split('@')[1]
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
    else:
        restart_builder.metadata.options.update({
            'resources': slurm_options[computer]['resources'],
            'max_wallclock_seconds': slurm_options[computer]['max_wallclock_seconds'],
            'account': slurm_options[computer]['account'],
            'scheduler_stderr': slurm_options[computer]['scheduler_stderr'],
            'scheduler_stdout': slurm_options[computer]['scheduler_stdout'],
            'queue_name': slurm_options[computer]['queue_name']
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

    # submit the calculation
    restart_builder.structure = structure
    restart_builder.kpoints = kpts
    restart_builder.parameters = parameters_default
    restart_builder.settings = Dict(dict={'CMDLINE': ['-npools', '4']})
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
    results_tmp[str(calc.uuid)]['previous_calc'] = uuid # pk is the previous calculation
    results_tmp[str(calc.uuid)]['son_calc'] = None # right now there is no son_node

    # change son_calc with previous simulation
    results_tmp[uuid]['son_calc'] = calc.uuid

    return results_tmp, calc.uuid
