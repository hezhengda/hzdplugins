# all the keys in the results dictionary
results_keys_set = ['system', 'comp_type', 'cluster', 'xc functional', 'exit_status', 'is_finished',
                    'is_finished_ok', 'E/eV', 'remove_remote_folder', 'previous_calc', 'son_calc']

# dictionary for the slurm system on different clusters
slurm_options = {}

# The following are the settings for the macOS settings

# 'claix-mac'
slurm_options['rwth-claix-mac'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on rwth_claix, each node has 48 cores.
        'max_wallclock_seconds': 86400,
        'account': 'jara0037',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'c18m'
    },
    'projwfc': {
        'resources': {'num_machines': 1},  # on rwth_claix, each node has 48 cores.
        'max_wallclock_seconds': 86400,
        'account': 'jara0037',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'c18m'
    }
}

# 'juwels-mac'
slurm_options['juwels-mac'] = {
    'qe': {
        'resources': {'num_machines': 4},  # on juwels, each node has 48 cores.
        'max_wallclock_seconds': 86400,
        'account': 'fzj-mac',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'batch'
    }
}

# 'jusuf-mac'
slurm_options['jusuf-mac'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on jusuf, each node has 128 cores
        'max_wallclock_seconds': 86400,
        'account': 'jiek61',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'batch'
    }
}

# 'jureca-dc-mac'
slurm_options['jureca-dc-mac'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on jureca_booster, each node has 128 cores
        'max_wallclock_seconds': 86400,
        'account': 'jiek61',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'dc-cpu'
    }
}

# 'jureca-booster-mac'
slurm_options['jureca-booster-mac'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on jureca_booster, each node has 68 cores
        'max_wallclock_seconds': 86400,
        'account': 'jiek61',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'booster'
    }
}

# The following are the settings for the ubuntu settings

# 'claix'
slurm_options['rwth-claix'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on rwth_claix, each node has 48 cores.
        'max_wallclock_seconds': 86400,
        'account': 'jara0037',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'c18m'
    }
}

# 'juwels'
slurm_options['juwels'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on juwels, each node has 48 cores.
        'max_wallclock_seconds': 86400,
        'account': 'fzj-mac',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'batch'
    }
}

# 'jusuf'
slurm_options['jusuf'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on jusuf, each node has 128 cores
        'max_wallclock_seconds': 86400,
        'account': 'jiek61',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'batch'
    }
}

# 'jureca_booster'
slurm_options['jureca-booster'] = {
    'qe': {
        'resources': {'num_machines': 2},  # on jureca_booster, each node has 68 cores
        'max_wallclock_seconds': 86400,
        'account': 'jiek61',
        'scheduler_stderr': 'stderr',
        'scheduler_stdout': 'stdout',
        'queue_name': 'booster'
    }
}

# the color scheme for plotting on the Matplotlib
color_dictionary = {
    'Ac': [112, 171, 250],
    'Ag': [192, 192, 192],
    'Al': [129, 178, 214],
    'Am': [84, 92, 242],
    'Ar': [207, 254, 196],
    'As': [116, 208, 87],
    'At': [117, 79, 69],
    'Au': [255, 209, 35],
    'B': [31, 162, 15],
    'Ba': [0, 201, 0],
    'Be': [94, 215, 123],
    'Bh': [224, 0, 56],
    'Bi': [158, 79, 181],
    'Bk': [138, 79, 227],
    'Br': [126, 49, 2],
    'C': [76, 76, 76],
    'Ca': [90, 150, 189],
    'Cd': [255, 217, 143],
    'Ce': [255, 255, 199],
    'Cf': [161, 54, 212],
    'Cl': [49, 252, 2],
    'Cm': [120, 92, 227],
    'Co': [0, 0, 175],
    'Cr': [0, 0, 158],
    'Cs': [87, 23, 143],
    'Cu': [34, 71, 220],
    'Db': [209, 0, 79],
    'Dy': [31, 255, 199],
    'Er': [0, 230, 117],
    'Es': [179, 31, 212],
    'Eu': [97, 255, 199],
    'F': [176, 185, 230],
    'Fe': [181, 113, 0],
    'Fm': [179, 31, 186],
    'Fr': [66, 0, 102],
    'Ga': [158, 227, 115],
    'Gd': [69, 255, 199],
    'Ge': [126, 110, 166],
    'H': [255, 204, 204],
    'He': [252, 232, 206],
    'Hf': [77, 194, 255],
    'Hg': [184, 184, 208],
    'Ho': [0, 255, 156],
    'Hs': [230, 0, 46],
    'I': [148, 0, 148],
    'In': [166, 117, 115],
    'Ir': [23, 84, 135],
    'K': [161, 33, 246],
    'Kr': [250, 193, 243],
    'La': [90, 196, 73],
    'Li': [134, 223, 115],
    'Lr': [199, 0, 102],
    'Lu': [0, 171, 36],
    'Md': [179, 13, 166],
    'Mg': [251, 123, 21],
    'Mn': [167, 8, 157],
    'Mo': [84, 181, 181],
    'Mt': [235, 0, 38],
    'N': [176, 185, 230],
    'Na': [249, 220, 60],
    'Nb': [115, 194, 201],
    'Nd': [199, 255, 199],
    'Ne': [254, 55, 181],
    'Ni': [183, 187, 189],
    'No': [189, 13, 135],
    'Np': [0, 128, 255],
    'O': [254, 3, 0],
    'Os': [38, 102, 150],
    'P': [192, 156, 194],
    'Pa': [0, 161, 255],
    'Pb': [87, 89, 97],
    'Pd': [0, 105, 133],
    'Pm': [163, 255, 199],
    'Po': [171, 92, 0],
    'Pr': [217, 255, 199],
    'Pt': [208, 208, 224],
    'Pu': [0, 107, 255],
    'Ra': [0, 125, 0],
    'Rb': [112, 46, 176],
    'Re': [38, 125, 171],
    'Rf': [204, 0, 89],
    'Rh': [10, 125, 140],
    'Rn': [66, 130, 150],
    'Ru': [36, 143, 143],
    'S': [255, 250, 0],
    'Sb': [158, 99, 181],
    'Sc': [181, 99, 171],
    'Se': [154, 239, 15],
    'Sg': [217, 0, 69],
    'Si': [27, 59, 250],
    'Sm': [143, 255, 199],
    'Sn': [154, 142, 185],
    'Sr': [0, 255, 0],
    'Ta': [77, 166, 255],
    'Tb': [48, 255, 199],
    'Tc': [59, 158, 158],
    'Te': [212, 122, 0],
    'Th': [0, 186, 255],
    'Ti': [120, 202, 255],
    'Tl': [166, 84, 77],
    'Tm': [0, 212, 82],
    'U': [0, 143, 255],
    'V': [229, 25, 0],
    'W': [33, 148, 214],
    'Xe': [66, 158, 176],
    'Y': [148, 255, 255],
    'Yb': [0, 191, 56],
    'Zn': [143, 143, 129],
    'Zr': [0, 255, 0],
}

cmd_shortcut = {
    'gf': "grep 'scf accuracy' aiida.out",
    'gforce': "grep 'Total force' aiida.out"
}


def molToMolecule(filename):
    """

    :code:`molToMolecule` can help us read the .mol file and return a pymatgen Molecule object.

    :param filename: The filename (relative path) that we want.
    :type filename: python string object

    :returns: A pymatgen Molecule object

    """

    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.io import read

    mol = read(filename, format='mol')

    return AseAtomsAdaptor.get_molecule(mol)


# all the adsorbates that I have been interested in so far:
import pkg_resources


def path(filename):
    return pkg_resources.resource_filename('hzdplugins', 'aiidaplugins/adsorbates/' + filename)


adsorbates = {
    'OH': {'mol': molToMolecule(filename=path('OH.mol')), 'ads_site': [0]},
    'OOH': {'mol': molToMolecule(filename=path('OOH.mol')), 'ads_site': [2]},

    'CH3': {'mol': molToMolecule(filename=path('CH3.mol')), 'ads_site': [0]},
    'CH2': {'mol': molToMolecule(filename=path('CH2.mol')), 'ads_site': [0]},
    'CH': {'mol': molToMolecule(filename=path('CH.mol')), 'ads_site': [0]},

    'NH2': {'mol': molToMolecule(filename=path('NH2.mol')), 'ads_site': [0]},
    'NH': {'mol': molToMolecule(filename=path('NH.mol')), 'ads_site': [0]},

    'CH2OH': {'mol': molToMolecule(filename=path('CH2OH.mol')), 'ads_site': [0]},
    'CH3O': {'mol': molToMolecule(filename=path('CH3O.mol')), 'ads_site': [1]},
    'CHOH': {'mol': molToMolecule(filename=path('CHOH.mol')), 'ads_site': [0]},
    'COH': {'mol': molToMolecule(filename=path('COH.mol')), 'ads_site': [0]},
    'CHO': {'mol': molToMolecule(filename=path('CHO.mol')), 'ads_site': [0]},
    'CO': {'mol': molToMolecule(filename=path('CO.mol')), 'ads_site': [0]},

    'HCOO': {'mol': molToMolecule(filename=path('HCOO.mol')), 'ads_site': [2]},
    'COOH': {'mol': molToMolecule(filename=path('COOH.mol')), 'ads_site': [0]}
}

pwParameter = {
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
        'input_dft': 'PBESOL'
    },
    'ELECTRONS': {
        'electron_maxstep': 200,
        'conv_thr': 1.0e-6,
        'diagonalization': 'david',
        'mixing_mode': 'plain',
        'mixing_beta': 0.3,
        'mixing_ndim': 10
    }
}

projwfcParameter = {
    'PROJWFC': {
        'DeltaE': 0.01,
        'ngauss': 0,
        'degauss': 0.015,
        'Emin': -40,
        'Emax': 40
    }
}

phParameter = {
    'INPUTPH': {
        'title_line': 'This is a ph.x calculation',
        'max_seconds': 86000,
        'tr2_ph': 1.0e-8,
        'ldisp': False,
        'epsil': False,
        'trans': True,
    }
}
