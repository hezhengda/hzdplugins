# This python package is used to generate inputs for the aiida-quantumespresso workchain.
# Workchain is a great tool for managing computational workflows, but the input of the workchain is really complex
# due to its layered structure, so this module we will create the similar structure as the aiida-quantumespresso and
# then help us to use the WorkChain more efficiently.

# author: Zheng-Da He
# initial time: 2020.12.19
# Location: Aachen, Germany

from aiida import orm

class PwCalculationInputGenerator():

    """
    The :code:`PwCalculationInputGenerator` class is used to construct the input file for PwCalculation.

    Initialization function for :code:`PwCalculationInputGenerator` class.

    :param code: The code we want to use for the pw.x calculation.
    :type code: aiida.orm.Code object

    :param structure: The structure we want to calculate
    :type code: aiida.orm.StructureData object

    :param pseudos: The pseudopotential we want to use

                    e.g. An example for pseudos

                    .. code-block:: python

                        pseudos = {
                            'Pt': UpfData(absolute_path_of_pseudopotential_for_Pt)
                            'O': UpfData(absolute_path_of_pseudopotential_for_O)
                        }
    :type pseudos: python dictionary object

    :param parameters: The parameters for the pw.x calculation
    :type parameters: aiida.orm.Dict object

    :param settings: The computational settings for the pw.x calculation
    :type settings: aiida.orm.Dict object

    :param metadata: The metadata for the pw.x calculation

                     e.g. An example for metadata

                     .. code-block:: python

                        metadata = {
                            'label': 'The label of your system, easy for query later',
                            'description': 'A description of your calculation',
                            'options': {
                                'resources': {'num_machines': x},
                                'max_wallclock_seconds': 86400,
                                'account': 'xxxxx',
                                'scheduler_stderr': 'stderr',
                                'scheduler_stdout': 'stdout',
                                'queue_name': 'xxxxxx'
                            }
                        }
    :type metadata: python dictionary object

    :param kpoints: kpoints for the simulation
    :type kpoints: aiida.orm.KpointsData
    """

    def __init__(self, code, structure, pseudos, parameters, settings, metadata, kpoints):

        self.code = code
        self.structure = structure
        self.pseudos = pseudos
        self.parameters = parameters
        self.settings = settings
        self.metadata = metadata
        self.kpoints = kpoints

    def outputdict(self):

        """

        :returns: A dictionary. Suitable for the :code:`PwCalculation`, can be used directly in the :code:`submit(
                  calculation, **inputs)`

        """
        tmp = {}

        tmp['code'] = self.code
        tmp['structure'] = self.structure
        tmp['pseudos'] = self.pseudos
        tmp['parameters'] = self.parameters
        tmp['settings'] = self.settings
        tmp['metadata'] = self.metadata
        tmp['kpoints'] = self.kpoints

        return tmp

class PwBaseWorkChainInputGenerator():

    """
    The :code:`PwBaseWorkChainInputGenerator` class is used to construct the input file for PwBaseWorkChain.

    Initialization function for :code:`PwBaseWorkChainInputGenerator` class.

    :param pw: The input of PwCalculation
    :type pw: PwCalculationInputGenerator object

    :param kpoints: The kpoints for the simulation
    :type kpoints: aiida.orm.KpointsData

    :param clean_workdir: If true, then we want to delete all the files in the work directory
    :type clean_workdir: python boolean object
    """

    def __init__(self, pw, kpoints, clean_workdir=True):

        self.pw = pw
        self.kpoints = kpoints
        self.clean_workdir = orm.Bool(clean_workdir)

    def outputdict(self):

        """

        :returns: A dictionary. Suitable for the :code:`PwBaseWorkChain`, can be used directly in the :code:`submit(
                  PwBaseWorkChain, **inputs)`

        """

        tmp = {}

        tmp['pw'] = self.pw.outputdict()

        if 'kpoints' in tmp['pw'].keys():
            tmp['pw'].pop('kpoints')

        tmp['kpoints'] = self.kpoints
        tmp['clean_workdir'] = self.clean_workdir

        return tmp

class PwRelaxWorkChainInputGenerator():

    """
    The :code:`PwRelaxWorkChainInputGenerator` class is used to construct the input file for PwBaseWorkChain.

    Initialization function for :code:`PwRelaxWorkChainInputGenerator` class.

    :param base: Input for the PwBaseWorkChain
    :type base: PwBaseWorkChainInputGenerator object

    :param base_final_scf: Input for the last scf simulation
    :type base_final_scf: PwBaseWorkChainInputGenerator object

    :param structure: The structure we want to relax
    :type structure: aiida.orm.StructureData object

    :param relax_type: How do we want to relax the structure. Default value can be summaried in below:

                       .. code-block:: python

                            relax_type_dict = {
                                'none': 'Nothing can move --> SCF simulation',
                                'atoms: 'Only atomic positions can be relaxed, cell is fixed.', ('relax')
                                'volume': 'Only volume can change, cell shape and atoms are fixed',
                                'shape': 'Only shape is optimized, volume and atomic positions are fixed',
                                'cell': 'Only cell is optimized (both shape and volume), atoms are fixed',
                                'atoms_volume': 'Relax atomic positions and volume',
                                'atoms_shape': 'Relax atomic positions and shape',
                                'atoms_cell': 'Relax both atomic positions and cell' ('vc-relax')
                            }
    :type relax_type: python string object

    :param meta_convergence: Don't know why this exists
    :type meta_convergence: python boolean object

    :param max_meta_convergence_iterations: pass
    :type max_meta_convergence_iterations: python int object

    :param volume_convergence: pass
    :type volume_convergence: python float object

    :param clean_workdir: pass
    :type clean_workdir: python boolean object
    """

    def __init__(self, base, base_final_scf, structure, relaxation_scheme, relax_type, meta_convergence=True,
                 max_meta_convergence_iterations=5, volume_convergence=0.01, clean_workdir=True):

        self.base = base
        self.base_final_scf = base_final_scf
        self.structure = structure
        # self.final_scf = orm.Bool(final_scf)
        self.relaxation_scheme = orm.Str(relaxation_scheme)
        self.relax_type = orm.Str(relax_type)
        self.meta_convergence = orm.Bool(meta_convergence)
        self.max_meta_convergence_iterations = orm.Int(max_meta_convergence_iterations)
        self.volume_convergence = orm.Float(volume_convergence)
        self.clean_workdir = orm.Bool(clean_workdir)

    def outputdict(self):

        """

        :returns: A dictionary. Suitable for the :code:`PwRelaxWorkChain`, can be used directly in the :code:`submit(
                  PwRelaxWorkChain, **inputs)`

        """

        tmp = {}

        tmp['base'] = self.base.outputdict()

        if 'clean_workdir' in tmp['base'].keys():
            tmp['base'].pop('clean_workdir')

        if 'structure' in tmp['base']['pw'].keys():
            tmp['base']['pw'].pop('structure')

        if 'parent_folder' in tmp['base']['pw'].keys():
            tmp['base']['pw'].pop('parent_folder')

        tmp['base_final_scf'] = self.base_final_scf.outputdict()

        if 'clean_workdir' in tmp['base_final_scf'].keys():
            tmp['base_final_scf'].pop('clean_workdir')

        if 'structure' in tmp['base_final_scf']['pw'].keys():
            tmp['base_final_scf']['pw'].pop('structure')

        if 'parent_folder' in tmp['base_final_scf']['pw'].keys():
            tmp['base_final_scf']['pw'].pop('parent_folder')

        tmp['structure'] = self.structure
        # tmp['final_scf'] = self.final_scf
        tmp['relaxation_scheme'] = self.relaxation_scheme
        tmp['relax_type'] = self.relax_type
        tmp['meta_convergence'] = self.meta_convergence
        tmp['max_meta_convergence_iterations'] = self.max_meta_convergence_iterations
        tmp['volume_convergence'] = self.volume_convergence
        tmp['clean_workdir'] = self.clean_workdir

        return tmp