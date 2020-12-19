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
    """

    def __init__(self, code, structure, pseudos, parameters, settings, metadata):

        self.code = code
        self.structure = structure
        self.pseudos = pseudos
        self.parameters = parameters
        self.settings = settings
        self.metadata = metadata

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
        tmp['settngs'] = self.settings
        tmp['metadata'] = self.metadata

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

        if 'kpoints' in tmp['pw'].keys():  # if 'kpoints' has been set, then delete it
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

    :param final_scf: If true, then we need to run the final scf simulation
    :type final_scf: aiida.orm.Bool

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
    :type relax_type: aiida.orm.Str object

    :param metadata_convergence: Don't know why this exists
    :type metadata_convergence: aiida.orm.Bool object

    :param meta_convergence_iterations: pass
    :type meta_convergence: aiida.orm.Int

    :param volume_convergence: pass
    :type volume_convergence: aiida.orm.Float object

    :param clean_workdir: pass
    :type clean_workdir: aiida.orm.Bool object
    """

    def __init__(self, base, base_final_scf, structure, relax_type, final_scf=True,
                 metadata_convergence=True,
                 meta_convergence_iterations=5, volume_convergence=0.01, clean_workdir=True):

        self.base = base
        self.base_final_scf = base_final_scf
        self.final_scf = final_scf
        self.structure = structure
        self.relax_type = relax_type
        self.metadata_convergence = orm.Bool(metadata_convergence)
        self.metadata_convergence_iterations = orm.Int(meta_convergence_iterations)
        self.volume_convergence = orm.Float(volume_convergence)
        self.clean_workdir = orm.Bool(clean_workdir)

    def outputdict(self):

        """

        :returns: A dictionary. Suitable for the :code:`PwRelaxWorkChain`, can be used directly in the :code:`submit(
                  PwRelaxWorkChain, **inputs)`

        """

        tmp = {}

        tmp['base'] = self.base.outputdict()

        # remove 'clean_workdir', 'pw.structure' and 'pw.parent_folder'
        if 'clean_workdir' in tmp['base'].keys():
            tmp['base'].pop('clean_workdir')

        if 'structure' in tmp['base']['pw'].keys():
            tmp['base']['pw'].pop('structure')

        if 'parent_folder' in tmp['base']['pw'].keys():
            tmp['base']['pw'].pop('parent_folder')

        tmp['base_final_scf'] = self.base_final_scf.outputdict()

        # remove 'clean_workdir', 'pw.structure' and 'pw.parent_folder'
        if 'clean_workdir' in tmp['base'].keys():
            tmp['base'].pop('clean_workdir')

        if 'structure' in tmp['base']['pw'].keys():
            tmp['base']['pw'].pop('structure')

        if 'parent_folder' in tmp['base']['pw'].keys():
            tmp['base']['pw'].pop('parent_folder')

        tmp['final_scf'] = self.final_scf
        tmp['structure'] = self.structure
        tmp['relax_type'] = self.relax_type
        tmp['metadata_convergence'] = self.metadata_convergence
        tmp['metadata_convergence_iterations'] = self.metadata_convergence_iterations
        tmp['volume_convergence'] = self.volume_convergence
        tmp['clean_workdir'] = self.clean_workdir

        return tmp