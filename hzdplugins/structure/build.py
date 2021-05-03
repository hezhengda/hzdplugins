# This module will be used in generating structures for Aiida simulation, it utilizes pymatgen, ase and catkit packages
# from other group, all the output structures should have Aiida StructureData type.
# In order to maintain provenance, we need to use all the types provided by the aiida.orm, and all functions should be
# decorated by @calcfunction decorator in order to be captured by the aiida program

from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from aiida.orm import StructureData, List, Dict
from ase.io import read
from matplotlib import patches
from matplotlib.path import Path
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, get_rot
from pymatgen.core.surface import SlabGenerator

from hzdplugins.aiidaplugins.constants import color_dictionary, adsorbates
from hzdplugins.aiidaplugins.info import getStructureAnalysis

def getValue(var):
    """

    :code:`getValue` can help us to initialize the variables

    :param var: variable that we want to get value from.
    :type var: aiida.orm Types

    :returns: the value of the variable.
    :rtype: corresponding python Types

    """

    from aiida.orm import Int, Float, Bool, Str, Dict, List

    if var is None:
        return None

    if isinstance(var, Int):
        return var.value

    if isinstance(var, Float):
        return var.value

    if isinstance(var, Bool):
        return var.value

    if isinstance(var, Str):
        return var.value

    if isinstance(var, Dict):
        return var.get_dict()

    if isinstance(var, List):
        return var.get_list()

def bulkFromFile(filename, supercell):
    """

    :code:`bulkFromFile` function can help us create a Bulk from the file.

    :param filename: Usually when we create the structural file, we do it from the strcutural file such as .xyz or
                     .cif, etc.
    :type filename: python string object

    :param supercell: A list that contains the dimension of supercell we would like.
    :type supercell: python list object

    :returns: A StructureData file that can be used directly in Aiida.
    :rtype: aiida.orm.StructurData object

    """

    # transfer from aiida.orm types to the common python types

    if len(filename) == 0:
        raise (IOError("You didn't provide an input file."))

    if len(filename.split('.')[1]) == 0:
        raise (IOError(
            "You didn't provide a correct file_type, please try to name your file .xyz, .cif or other data structure "
            "types."))

    file_type = filename.split('.')[1]

    structure = read(filename, format=file_type)

    return StructureData(ase=structure * supercell)

def bulkFromString(bulkStr, crystal_structure, a, cubic, supercell, b=None, c=None, alpha=None, covera=None, u=None,
                   orthorhombic=None):
    """

    :code:`bulkFromFile` function can help us create a Bulk from the string.

    :param bulkStr: The string for the material.
    :type bulkStr: python string object

    :param crystal_structure: The crystal structure. It need to be in one of those: sc, fcc,
                              bcc,  tetragonal, bct, hcp, rhombohedral, orthorhombic, mlc, diamond, zincblende,
                              rocksalt, cesiumchloride, fluorite or wurtzite.
    :type crystal_structure: python string object

    :param a(/b/c): Lattice constants.
    :type a(/b/c): python float object

    :param supercell: The supercell that we want to get.
    :type supercell: python list object

    :returns: The structure of the bulk.
    :rtype: aiida.orm.StructureData

    """

    from ase.build import bulk

    bulk_ase = bulk(name=bulkStr, crystalstructure=crystal_structure, a=a, b=b, c=c, alpha=alpha, covera=covera, u=u,
                    orthorhombic=orthorhombic, cubic=cubic)

    return StructureData(ase=bulk_ase * supercell)


# In research, not only we need to deal with solids, we also need to deal with surfaces, and add adsorbates on it,
# so in my module, it will be important to create any structures that I want easily and efficiently.

def millerSurfaces(bulk, miller_index, layers, vacuum, get_orthogonal=False, bonds=None):
    """

    :code:`millerSurfaces` can help us generate a list of slab structures that we could use in future studies

    :param bulk: The bulk structure that we are going to use for creating the surface slabs.
    :type bulk: aiida.orm.StructureData object

    :param miller_index: The miller index that we want to get.
    :type miller_index: python list object

    :param layers: Set how many layers you want for the surface slab.
    :type layers: python int object

    :param vacuum: Set how many layers you want for the vacuum layer.
    :type vacuum: python list object

    :param get_orthogonal: Whether we want to get orthogonal slab or not
    :type get_orthogonal: python boolean object

    :param bonds: a dictionary which define the bond length of two atoms. e.g. {('P', 'O'):3}
    :type bonds: python dictionary object

    :returns: A list of slabs that we generated, notice that all the slabs are orthogonal, because we use
              :code:`slab.get_orthogonal_c_slab()` for all the slabs
    :rtype: pymatgen.core.structure.slab object

    """

    # do the pre-process of aiida.orm objects.
    bulk_pmgstructure = bulk.get_pymatgen_structure()

    sg = SlabGenerator(initial_structure=bulk_pmgstructure,
                       miller_index=miller_index,
                       min_slab_size=layers,
                       min_vacuum_size=vacuum,
                       center_slab=True,
                       in_unit_planes=True,
                       primitive=False,
                       # max_normal_search=max(miller_index),
                       reorient_lattice=True)

    listOfStructures = sg.get_slabs(bonds=bonds)
    if get_orthogonal:
        listOfStructures = [slab.get_orthogonal_c_slab() for slab in listOfStructures]

    # check whether all the atoms are in the unit cell
    for slab in listOfStructures:
        for ind in range(len(slab.sites)):
            slab.sites[ind] = slab.sites[ind].to_unit_cell()

    # results = []
    # for structure in listOfStructures:
    #     structure_node = StructureData(pymatgen_structure=structure)
    #     structure_node.store()
    #     results.append(structure_node.uuid)  # since uuid is the unique identifier for each node across the platform
    #
    # listGenerator = List()
    # listGenerator.set_list(results)

    return listOfStructures

def adsorptionSites(slab, **kwargs):
    """

    :code:`AdsorptionSites` can help us visualize the position and tag of each adsorption sites, then we can determine
    where we want to put the adsorbates.

    :param slab: This is our slab, and we want to find how we can put the adsorbates on the slab.
    :type slab: aiida.orm.StructureData

    :param kwargs: * distance: the distance between adsorption site and the surface
                   * symm_reduce: the symmetry reduce (default = 0.01)
                   * near_reduce: the near reduce (default = 0.01)

    :returns: Dictionary contains the dictionary of adsorption sites.
    :rtype: aiida.orm.Dict object

    """

    # the inspiration for this function was from pymatgen.analysis.adsorption.AdsorbateSiteFinder.plot_slab()
    # function, which is really intuitive way of showing the structure and the adsorption sites.
    # since this function does not create any useful data, so it doesn't be decorated with @calfunction decorator

    # get the structure and adsorption sites
    slab = slab.get_pymatgen_structure()
    # end of the conversion

    asf = AdsorbateSiteFinder(slab, selective_dynamics=False)

    if 'distance' in kwargs.keys():
        distance = kwargs['distance']
    else:
        distance = 1.2

    if 'symm_reduce' in kwargs.keys():
        symm_reduce = kwargs['symm_reduce']
    else:
        symm_reduce = 0.01

    if 'near_reduce' in kwargs.keys():
        near_reduce = kwargs['near_reduce']
    else:
        near_reduce = 0.01

    adsorption_sites = asf.find_adsorption_sites(distance=distance, symm_reduce=symm_reduce, near_reduce=near_reduce)

    dictGenerator = Dict()
    dictGenerator.set_dict(adsorption_sites)

    return dictGenerator


def visualizeSlab(slab, plot_adsSite=False, adsorption_sites=None, adssitetype=['ontop', 'bridge', 'hollow'], **kwargs):
    """

    :code:`visualizeSlab` will show the slab

    :param slab: The slab that we want to visualize
    :type slab: aiida.orm.StructureData object

    :param plot_adsSite: If true, then add adsorption sites, if false, then adsorption site not show.
    :type plot_adsSite: python boolean object

    :param adsorption_sites: Shows the adsorption sites.
    :type adsorption_sites: aiida.orm.Dict object

    :param adssitetype: determine which adsorption site you want to plot, the default is all types of sites (ontop,
                        bridge, hollow), but you can specify on your own. Since sometime the program tends to give us
                        more sites, so it is not easy to see, so we can define what kind of site we want to investigate.
    :type adssitetype: python list object

    :param kwargs: Settings for the plot:
                   * repeat: Int
                   * decay: Float
                   * scale: Float
                   * window: Float
    :returns: A graph which represents the slab surface (and the position of adsorption sites).
    :rtype: matplotlib.pyplot object.

    """

    # start plotting the figure, this part of the code was largely adopted from pymatgen github repository

    slab = slab.get_pymatgen_structure()
    if plot_adsSite:
        adsorption_sites = adsorption_sites.get_dict()

    # parameters for the plotting
    if 'repeat' in kwargs.keys():
        repeat = kwargs['repeat']
    else:
        repeat = 2  # create supercell 3x3x1

    if 'decay' in kwargs.keys():
        decay = kwargs['decay']
    else:
        decay = 0.2

    if 'scale' in kwargs.keys():
        scale = kwargs['scale']
    else:
        scale = 0.8

    if 'window' in kwargs.keys():
        window = kwargs['window']
    else:
        window = 1.5

    draw_unit_cell = True

    fig, ax = plt.subplots(figsize=(10, 10))

    orig_slab = deepcopy(slab)
    orig_cell = deepcopy(slab.lattice.matrix)

    slab.make_supercell([repeat, repeat, 1])

    # sort the coordinates by the z component
    coords = np.array(sorted(slab.cart_coords, key=lambda x: x[2]))
    sites = sorted(slab.sites, key=lambda x: x.coords[2])

    alphas = 1 - decay * (np.max(coords[:, 2]) - coords[:, 2])
    alphas = alphas.clip(min=0)
    corner = [0, 0, slab.lattice.get_fractional_coords(coords[-1])[-1]]
    corner = slab.lattice.get_cartesian_coords(corner)[:2]
    verts = orig_cell[:2, :2]
    lattsum = verts[0] + verts[1]

    # Draw circles at sites and stack them accordingly
    for n, coord in enumerate(coords):
        r = sites[n].specie.atomic_radius * scale
        ax.add_patch(
            patches.Circle(
                coord[:2] - lattsum * (repeat // 2), r, color="w", zorder=2 * n
            )
        )
        color = np.array(color_dictionary[sites[n].species_string]) / 255
        ax.add_patch(
            patches.Circle(
                coord[:2] - lattsum * (repeat // 2),
                r,
                facecolor=color,
                alpha=alphas[n],
                edgecolor="k",
                lw=0.3,
                zorder=2 * n + 1,
            )
        )

        if plot_adsSite:
            # Adsorption sites
            # top site
            if 'ontop' in adssitetype:
                ads_sites = adsorption_sites['ontop']
                sop = get_rot(orig_slab)
                ads_sites = [sop.operate(ads_site)[:2].tolist() for ads_site in ads_sites]
                ax.plot(
                    *zip(*ads_sites),
                    color="k",
                    marker="o",
                    markersize=20,
                    mew=1,
                    linestyle="",
                    zorder=10000
                )
                for site_id, ads_site in enumerate(ads_sites):
                    ax.text(ads_site[0], ads_site[1],
                            str(site_id),
                            color='yellow',
                            fontsize=16,
                            ha='center', va='center',
                            zorder=20000)
            # bridge site
            if 'bridge' in adssitetype:
                ads_sites = adsorption_sites['bridge']
                sop = get_rot(orig_slab)
                ads_sites = [sop.operate(ads_site)[:2].tolist() for ads_site in ads_sites]
                ax.plot(
                    *zip(*ads_sites),
                    color="k",
                    marker="s",
                    markersize=20,
                    mew=1,
                    linestyle="",
                    zorder=10000
                )
                for site_id, ads_site in enumerate(ads_sites):
                    ax.text(ads_site[0], ads_site[1],
                            str(site_id),
                            color='yellow',
                            fontsize=16,
                            ha='center', va='center',
                            zorder=20000)
            # hollow site
            if 'hollow' in adssitetype:
                ads_sites = adsorption_sites['hollow']
                sop = get_rot(orig_slab)
                ads_sites = [sop.operate(ads_site)[:2].tolist() for ads_site in ads_sites]
                ax.plot(
                    *zip(*ads_sites),
                    color="k",
                    marker="^",
                    markersize=20,
                    mew=1,
                    linestyle="",
                    zorder=10000
                )
                for site_id, ads_site in enumerate(ads_sites):
                    ax.text(ads_site[0], ads_site[1],
                            str(site_id),
                            color='yellow',
                            fontsize=16,
                            ha='center', va='center',
                            zorder=20000)
        else:
            pass

    # Draw unit cell
    if draw_unit_cell:
        verts = np.insert(verts, 1, lattsum, axis=0).tolist()
        verts += [[0.0, 0.0]]
        verts = [[0.0, 0.0]] + verts
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        verts = [(np.array(vert) + corner).tolist() for vert in verts]
        path = Path(verts, codes)
        patch = patches.PathPatch(
            path, facecolor="none", lw=2, alpha=0.5, zorder=2 * len(coords) + 2
        )
        ax.add_patch(patch)

    ax.set_aspect("equal")
    center = corner + lattsum / 2.0
    extent = np.max(lattsum)
    lim_array = [center - extent * window, center + extent * window]
    x_lim = [ele[0] for ele in lim_array]
    y_lim = [ele[1] for ele in lim_array]
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)

    plt.show()

def addAdsorbates(slab, adsSiteDictionary):
    """

    :code:`addAdsorbates` will take care of these things:

    * Add each adsorbate to a specific site, and with optimal adsorption distance and angle
    * Put the bottom two layers of the slab fixed, other atoms can be relaxed

    :param slab: The slab that we want to adsorb some adsorbates on.
    :type slab: aiida.orm.StructureData object

    :param adsSiteDictionary: In which contains the adsorbates (for common adsorbates, I can assign the
                              adsorption atom, geometry of the adsorbates and its adsorption angle in
                              :code:`constants.py`) and adsorption sites, also I can add new adsorbates to the
                              adsorbates library. (Do a very simple relax simulation of a molecule in the big cell,
                              and just output the structure and assign the adsorption site. For mono-site adsorption,
                              it is really easy, for bi-site adsorption, it is also easy since we can get the
                              molecule to align; for tri-site or more-site adsorption, well, we can just put the
                              molecule closer to the surface and let the program determine how does this molecule
                              interact with the surface, for our usual research, I don't think this is necessary.)

                              .. code-block:: python

                                    adsSiteDictionary = {
                                        'O': [site1, site2, site3] # each site is a 3x1 list [a, b, c]
                                        'F': [site1, site2, site3]
                                    }
    :type adsSiteDictionary: python dictionary object

    :returns: With the adsorbates and modified constrains on all the atoms, ready for the submit functions.
    :rtype: aiida.orm.StructureData object

    """

    # clean the input parameters from aiida.orm types to usual types
    slab = slab.get_pymatgen_structure()
    slab_tmp = deepcopy(slab)
    # asf = AdsorbateSiteFinder(slab_tmp, selective_dynamics = True)
    # end of cleaning process

    for ads, siteList in adsSiteDictionary.items():
        adsorbate = getMoleculeByName(ads)
        if len(ads) == 1:
            ads_site = [0]
        else:
            ads_site = adsorbates[ads]['ads_site']
        if len(ads_site) == 1:
            for site in siteList:
                slab_tmp = hzd_add_adsMono(slab_tmp, molecule=adsorbate, ads_coord=[site], ads_site=ads_site)
        # elif len(ads_site) == 2:
        #     # treate it as mono-dent
        #     for site in siteList:
        #         for adsite in ads_site:
        #             slab_tmp = hzd_add_adsMono(slab_tmp, molecule = adsorbate, ads_coord = site, ads_site = ads_site)
        #     slab_tmp = hzd_add_adsBi(slab_tmp, molecule = adsorbate, ads_coord = siteList, ads_site = ads_site)

    return StructureData(pymatgen_structure=slab_tmp)

def newStructure(structure, changeDict):
    """

    :code:`newStructure` will create a new structure by replacing the changeList atom to new type of atoms,
    very useful if there are different chemical states of the same atom in the structure.

    :param structure: previous structure with conventional atomic symbols
    :type structure: aiida.orm.StructureData object

    :param changeDict: A dictionary where the i-th atom's label will be replaced. e.g. 'Fe' to 'Fe2'
    :type changeDict: python dictionary

    :returns: return an object which stores the change of the label
    :rtype: aiida.orm.StructureData object

    """

    # create the kind_names list
    kind_names = []
    structure_ase = structure.get_ase()

    for atom in structure_ase:
        kind_names.append(atom.symbol)

    for key, value in changeDict.items():
        for ind, atomSymbol in enumerate(kind_names):
            if atomSymbol == key:
                kind_names[ind] = value

    # since I recreate kind_names from structure, then they must be the same.
    # if len(structure.sites) != len(kind_names):
    #     raise ValueError('The number of new kind names must be equal to the number of sites in the structure.')

    new_structure = StructureData(
        cell=structure.cell,
        pbc=structure.pbc,
    )

    for site, kind_name in zip(structure.sites, kind_names):
        kind = structure.get_kind(site.kind_name)

        new_structure.append_atom(
            name=kind_name,
            symbols=kind.symbols,
            weights=kind.weights,
            position=site.position
        )

    return new_structure

# The functions below need to be used with care, because they are not part of the open API.
def getMoleculeByName(ads_str):
    """

    :code:`getMoleculeByName` can return the Molecule object of a specific adsorbates, could be very useful for the
    simulations.

    :param ads_str: an index of adsorbates in the database
    :type ads_str: python string

    :returns: An atom or a molecule
    :rtype: pymatgen molecule object

    """

    from pymatgen.core import Element, Molecule
    from hzdplugins.aiidaplugins.constants import adsorbates

    if len(ads_str) == 1:  # means it is an atom
        if Element(ads_str):  # means ads_str is an element
            return Molecule(species=ads_str, coords=[[0, 0, 0]])
        else:
            raise ValueError("You have entered a wrong string for an element, please try again.")
    else:
        if ads_str in adsorbates.keys():
            return adsorbates[ads_str]['mol']
        else:
            raise ValueError("Sorry, the adsorbates you are asking for haven't been added in the package yet.")

def setFixedCoords(slab, height):
    """

    :code:`setFixedCoords` is a function that can return a list that shows which atom we want to freeze during the
    simulation

    :param slab: The slab that we want to fixed atoms.
    :type slab: aiida.orm.StructureData object

    :param height: shows the boundary of the fixed atom. if the z-component of the position of the atom is below
                   height, then we freeze it.
    :type height: float

    :returns: A Nx3 list which holds all the information about the fixed atoms, can be directly passed to the
              setting_dict.
    :rtype: python list

    """

    slab_ase = slab.get_ase()
    fixed_coords = []

    for atom in slab_ase:
        if atom.position[2] > height:
            fixed_coords.append([False, False, False])
        else:
            fixed_coords.append([True, True, True])

    return fixed_coords

def hzd_add_adsMono(slab, molecule, ads_coord, ads_site):
    """

    :code:`hzd_add_ads` is a help function that can help me add adsorbates.

    :param slab: A slab that we are interested in.
    :type slab: Pymatgen Structure object.

    :param molecule: The adsorbate that we want to add
    :type slab: Pymatgen Molecule object.

    :param ads_coord: The position of adsorption site
    :type ads_coord: python list

    :param ads_site: The adsorption site of the molecule. If there is only 1 site in the molecule, that means it is
                     mono-ads, if there are 2 sites in the molecule, then it can be mono-ads or bi-ads (needs to
                     treat differently).

    :returns: A modified slab with added adsorbates.
    :rtype: pymatgen structure object.

    """

    slab_tmp = deepcopy(slab)
    molecule_tmp = deepcopy(molecule)

    # check whether ads_coord and ads_site are compatible:
    if len(ads_coord) == len(ads_site):
        pass
    else:
        raise ValueError(
            "Sorry, the adsorption site on the slab and on the molecule are different, they must have the same length.")

    if len(ads_site) == 1:  # momo-ads
        coord_ads_site = molecule[ads_site[0]].coords
        vec = ads_coord[0] - coord_ads_site
        for site in molecule_tmp:
            site.coords = site.coords + vec
        # molecule_tmp = hzd_rotate(molecule_tmp, ads_site[0]) # make sure that the ads_site[0] are in the lowest point.
        molecule_tmp.add_site_property(
            "surface_properties", ['adsorbate'] * molecule_tmp.num_sites
        )
        molecule_tmp.add_site_property(
            'selective dynamics', [[True, True, True]] * molecule_tmp.num_sites
        )
        for site in molecule_tmp:
            slab_tmp.append(species=site.specie, coords=site.coords, coords_are_cartesian=True,
                            properties=site.properties)

    return slab_tmp

def setSpinStructure(symbol, hl_spin, num_electrons):
    """
    :code:`setSpinStructure` can help us generate the :code:`starting_ns_eigenvalue` quickly, otherwise it will take too long. Notice in here we only assume that all the symbol are d-groups (just for simplicity for now.)

    :param symbol: The symbol of our atom, e.g. 'Fe'
    :type symbol: python string object
    :param hl_spin: whether we want low_spin ('ls') or high_spin ('hs')
    :type hl_spin: python string object
    :param num_electrons: the number of electrons that we want to assign
    :type num_electrons: python int object

    :raises ValueError: If your :code:`num_electrons` is larger than 10, then it is not possible, because the maximum amount of d electrons is 10.
    :return: A list of lists which contains the information about the spin configuration.
    :rtype: python list object
    """
    results = []

    if num_electrons > 10:
        raise ValueError('The amount of electrons in d orbital should be smaller than 10.')

    tmp = num_electrons

    # assign electrons
    if hl_spin == 'hs': # high spin
        m = 1
        s = 1
        for i in range(tmp):
            tmplist = [m, s, symbol, 1.0]
            results.append(tmplist)
            m += 1
            if m > 5:
                m = 1
                s = 2
    elif hl_spin == 'ls':
        m = 1
        s = 1
        for i in range(tmp):
            tmplist = [m, s, symbol, 1.0]
            s += 1
            if s == 3:
                m += 1
                s = 1
            results.append(tmplist)

    for m in range(5):
        for s in range(2):
            tmplist = [m+1, s+1, symbol, 1.0]
            if tmplist in results:
                pass
            else:
                results.append([m+1, s+1, symbol, 0.0])

    return results

def delAtoms(structure, atom_list):
    """
    :code:`delAtoms` can delete any atoms you want, and return the aiida.orm.StructureData object.

    :param structure: The structure that we want to deal with
    :type structure: aiida.orm.StructureData
    :param atom_list: The list of atoms that we want to delete
    :type atom_list: python list object
    :return: A new structure where the atom in the atom_list has been deleted.
    :rtype: aiida.orm.StructureData
    """

    from aiida.orm import StructureData
    # get ase structure
    str_ase = structure.get_ase()

    # del atoms
    del str_ase[atom_list]

    return StructureData(ase=str_ase)

def atomQuery(structure, query_dict, is_and=True):
    """
    :code:`atomQuery` can help us query the list of atoms that we want to manipulate later.

    :param structure: The structure that we want to investigate.
    :type structure: aiida.orm.StructureData
    :param query_dict: The query dictionary, now support these keys: (1) 'symbol' (2) 'x', 'y', 'z' (3) 'connection'

                       .. code-block:: python

                            query_dict = {
                                'symbol': ['Ni'],
                                'z': {'>', 30.0},
                                'connection': 25
                            }

                       which means we want to query Ni atoms, that z coordinates is larger than 30.0, also have connection to 25th atom.
    :type query_dict: python dictionary object
    :param is_and: whether we want to and the condition in query_dict, defaults to True
    :type is_and: python bool object, optional
    :return: index of atoms
    :rtype: python list object
    """

    results = {} # we want to return a list of symbols
    tmp_ase = structure.get_ase()

    for item, value in query_dict.items():

        # to query certain atom symbol
        if item == 'symbol':
            results['symbol'] = []
            for ind, atom in enumerate(tmp_ase):
                if atom.symbol in value:
                    results['symbol'].append(ind)

        # to query certain position
        if item in ['x', 'y', 'z']:
            results['coords'] = []
            for ind, atom in enumerate(tmp_ase):
                if list(value.keys())[0] == '>':
                    if item == 'x':
                        if atom.position[0] > list(value.values())[0]:
                            results['coords'].append(ind)
                    elif item == 'y':
                        if atom.position[1] > list(value.values())[0]:
                            results['coords'].append(ind)
                    elif item == 'z':
                        if atom.position[2] > list(value.values())[0]:
                            results['coords'].append(ind)

                if list(value.keys())[0] == '<':
                    if item == 'x':
                        if atom.position[0] < list(value.values())[0]:
                            results['coords'].append(ind)
                    elif item == 'y':
                        if atom.position[1] < list(value.values())[0]:
                            results['coords'].append(ind)
                    elif item == 'z':
                        if atom.position[2] < list(value.values())[0]:
                            results['coords'].append(ind)

        # to query connection
        if item == 'connection':
            results['connection'] = []
            connected_atoms = getStructureAnalysis(structure, atom_index=[value])
            for ind, atom in enumerate(tmp_ase):
                atom_str = atom.symbol + str(ind) + 'ax'
                is_in = False
                for name in list(list(connected_atoms.values())[0].keys()):
                    if atom_str in name:
                        is_in = True
                        break
                if is_in:
                    results['connection'].append(ind)

    # analyze all the results
    tmp = list(range(len(structure.sites))) #create the index of all atoms
    if is_and:
        for key, value in results.items():
            tmp = list(set(tmp) & set(value))

    return tmp

def expandLayerDistance(struct, begin_end_layer, setDistance):
    
    """
    expandLayerDistance can expand the distance between layers in the layered structure
    
    :param struct: The structure
    :type struct: aiida.orm.StructureData
    
    :param begin_end_layer: The beginning and the ending of each layer
    :type struct: python list
    
    e.g. begin_end_layer = [[a1, b1], [a2, b2], [a3, b3]]
    
    :param setDistance: The distance between layers that we want to set
    :type setDistance: python float
    
    """
    
    tmp_ase = struct.get_ase()

    originalDistance = begin_end_layer[1][0] - begin_end_layer[0][1]
    diff = setDistance - originalDistance
    c_new = tmp_ase.cell[2][2] + len(begin_end_layer) * diff

    tmp_ase.cell[2] = [0, 0, c_new]

    atomLayers = [] # used to store the atom in different layers

    # divide the system into different layers
    for begin_end in begin_end_layer:
        tmp_layer = []
        for atom in tmp_ase:
            if atom.position[2] > begin_end[0] and atom.position[2] < begin_end[1]:
                tmp_layer.append(atom)
        atomLayers.append(tmp_layer)

    # add the distance to each layer
    for ind, atoms in enumerate(atomLayers):
        print(ind)
        for atom in atoms:
            atom.position += [0, 0, ind*diff]

    struct = StructureData(ase=tmp_ase)
    
    return struct

def buildMoleculeFromSMILE(smileStr):
    """
    Create molecular structure by using 

    :param smileStr: The string of the SMILES
    :type smileStr: Python string

    :return: A molecule structure
    :rtype: ase.Atoms object
    """

    from openbabel import openbabel
    from ase.io import read, write 
    import numpy as np
    import os

    f = open('babel.xyz', 'w')
    gen3d = openbabel.OBOp.FindType('gen3D')
    mol = openbabel.OBMol()

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('smi', 'xyz')
    obConversion.ReadString(mol, smileStr)

    gen3d.Do(mol, '--best')
    outMDL = obConversion.WriteString(mol)
    f.write(outMDL)
    f.close()

    atoms = read('babel.xyz')
    os.system('rm babel.xyz')

    return atoms