# This module will be used in generating structures for Aiida simulation, it utilizes pymatgen, ase and catkit packages from other group, all the output structures should have Aiida StructureData type.
# In order to maintain provenance, we need to use all the types provided by the aiida.orm, and all functions should be decorated by @calcfunction decorator in order to be captured by the aiida program

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.path import Path
from copy import deepcopy

from hzdplugins.aiidaplugins.constants import color_dictionary, adsorbates
from hzdplugins.structure.math import rotation_matrix_euler

from aiida.orm import StructureData, Str, List, Bool, Dict
from aiida.engine import calcfunction

from ase.io import read

from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, get_rot

def getValue(var):

    """

    :code:`getValue` can help us to initialize the variables

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

@calcfunction
def bulkFromFile(filename, supercell):

    """

    `bulkFromFile` function can help us create a Bulk from the file.

    Parameters:

    filename:
        An aiida.orm.Str object. Usually when we create the structural file, we do it from the strcutural file such as .xyz or .cif, etc.

    supercell:
        An aiida.orm.List object. A list that contains the dimension of supercell we would like.

    Return: A StructureData file that can be used directly in Aiida.

    """

    # transfer from aiida.orm types to the common python types
    filename = getValue(filename)
    supercell = getValue(supercell)

    if len(filename)==0:
        raise(IOError("You didn't provide an input file."))

    if len(filename.split('.')[1])==0:
        raise(IOError("You didn't provide a correct type, please try to name your file .xyz, .cif or other data structure types."))

    type = filename.split('.')[1]

    structure = read(filename, format=type)

    return StructureData(ase=structure*supercell)

@calcfunction
def bulkFromString(bulkStr, crystal_structure, a, cubic, supercell, b=None, c=None, alpha=None, covera=None, u=None, orthorhombic=None):

    """

    :code:`bulkFromFile` function can help us create a Bulk from the string.

    Parameters:

    bulkStr:
        An aiida.orm.Str object. The string for the material

    crystal_structure:
        An aiida.orm.Str object. The crystal structure. It need to be in one of those: sc, fcc, bcc, tetragonal, bct, hcp, rhombohedral, orthorhombic, mlc, diamond, zincblende, rocksalt, cesiumchloride, fluorite or wurtzite.

    a, b, c:
        An aiida.orm.Float object. Lattice constants.

    supercell:
        An aiida.orm.List object. The supercell that we want to get.

    """

    from ase.build import bulk

    # convert all aiida
    bulkStr = getValue(bulkStr)
    crystal_structure = getValue(crystal_structure)
    a = getValue(a)
    b = getValue(b)
    c = getValue(c)
    alpha = getValue(alpha)
    covera = getValue(covera)
    u = getValue(u)
    orthorhombic = getValue(orthorhombic)
    cubic = getValue(cubic)
    supercell = getValue(supercell)

    bulk_ase = bulk(name=bulkStr, crystalstructure=crystal_structure, a=a, b=b, c=c, alpha=alpha, covera=covera, u=u, orthorhombic=orthorhombic, cubic=cubic)

    return StructureData(ase=bulk_ase*supercell)

# In research, not only we need to deal with solids, we also need to deal with surfaces, and add adsorbates on it, so in my module, it will be important to create any structures that I want easily and efficiently.

@calcfunction
def millerSurfaces(bulk, miller_index, layers, vacuum, **kwargs):

    """

    `millerSurfaces` can help us generate a list of slab structures that we could use in future studies

    Parameters:

    bulk:
        An aiida.orm.StructureData object. The bulk structure that we are going to use for creating the surface slabs.

    miller_index:
        An aiida.orm.List object. The miller index that we want to get.

    layers:
        An aiida.orm.Int object. Set how many layers you want for the surface slab.

    vacuum:
        An aiida.orm.List object. Set how many layers you want for the vacuum layer.

    Return: A series of uuids of StructureData object, which can be used for further investigations.
        The reason why we use uuid instead of StructureData is that StructureData cannot be put in a list with List object, it will shout out error. Once we have the uuid number, we can get the structure by just using `structure = load_node(uuid=uuid_structure)`, it is very straightforward.

    """

    # do the pre-process of aiida.orm objects.
    bulk_pmgstructure = bulk.get_pymatgen_structure()
    miller_index_list = miller_index.get_list()
    layers_int = layers.value
    vacuum_int = vacuum.value

    sg = SlabGenerator(initial_structure = bulk_pmgstructure,
                       miller_index = miller_index_list,
                       min_slab_size = layers_int,
                       min_vacuum_size = vacuum_int,
                       center_slab = True,
                       in_unit_planes = True,
                       primitive = False,
                       max_normal_search = max(miller_index_list),
                       reorient_lattice = True)

    listOfStructures = sg.get_slabs()

    results = []
    for structure in listOfStructures:
        structure_node = StructureData(pymatgen_structure = structure)
        structure_node.store()
        results.append(structure_node.uuid) # since uuid is the unique identifier for each node across the platform

    listGenerator = List()
    listGenerator.set_list(results)

    return listGenerator

@calcfunction
def adsorptionSites(slab, **kwargs):

    """

    `AdsorptionSites` can help us visualize the position and tag of each adsorption sites, then we can determine where we want to put the adsorbates.

    Parameters:

    slab:
        An aiida.orm.StructureData object. This is our slab, and we want to find how we can put the adsorbates on the slab.

    kwargs:
        * distance: the distance between adsorption site and the surface
        * symm_reduce: the symmetry reduce (default = 0.01)
        * near_reduce: the near reduce (default = 0.01)

    Return:

    Adsorption site dictionary:
        An aiida.orm.Dict object. Contains the information of all the positions of adsorption sites, which corresponse to the matplotlib figure.

    """

    # the inspiration for this function was from pymatgen.analysis.adsorption.AdsorbateSiteFinder.plot_slab() function, which is really intuitive way of showing the structure and the adsorption sites.
    # since this function does not create any useful data, so it doesn't be decorated with @calfunction decorator

    # get the structure and adsorption sites
    slab = slab.get_pymatgen_structure()
    # end of the conversion

    asf = AdsorbateSiteFinder(slab, selective_dynamics = False)

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

    adsorption_sites = asf.find_adsorption_sites(distance = distance, symm_reduce = symm_reduce, near_reduce = near_reduce)

    dictGenerator = Dict()
    dictGenerator.set_dict(adsorption_sites)

    return dictGenerator

def visualizeSlab(slab, plot_adsSite=False, adsorption_sites=None, **kwargs):

    """

    :code:`visualizeSlab` will show the slab

    Parameters:

    slab:
        An aiida.orm.StructureData object.

    plot_adsSite:
        An aiida.orm.Bool object. If true, then add adsorption sites, if false, then adsorption site not show

    adsorption_sites:
        An aiida.orm.Dict object. Shows the adsorption sites.

    kwargs:
        Settings for the plot:
        * repeat: Int
        * decay: Float
        * scale: Float
        * window: Float

    """

    # start plotting the figure, this part of the code was largely adopted from pymatgen github repository

    slab = slab.get_pymatgen_structure()
    if plot_adsSite:
        adsorption_sites = adsorption_sites.get_dict()

    # parameters for the plotting
    if 'repeat' in kwargs.keys():
        repeat = kwargs['repeat']
    else:
        repeat = 2 # create supercell 3x3x1

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
        color = np.array(color_dictionary[sites[n].species_string])/255
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
            for id, ads_site in enumerate(ads_sites):
                ax.text(ads_site[0], ads_site[1],
                        str(id),
                        color='yellow',
                        fontsize=16,
                        ha='center',va='center',
                        zorder=20000)
            # bridge site
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
            for id, ads_site in enumerate(ads_sites):
                ax.text(ads_site[0], ads_site[1],
                        str(id),
                        color='yellow',
                        fontsize=16,
                        ha='center',va='center',
                        zorder=20000)
            # hollow site
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
            for id, ads_site in enumerate(ads_sites):
                ax.text(ads_site[0], ads_site[1],
                        str(id),
                        color='yellow',
                        fontsize=16,
                        ha='center',va='center',
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
            path, facecolor="none", lw=2, alpha=0.5, zorder=2 * n + 2
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

@calcfunction
def addAdsorbates(slab, adsSiteDictionary):

    """

    :code:`addAdsorbates` will take care of these things:

    * Add each adsorbate to a specific site, and with optimal adsorption distance and angle
    * Put the bottom two layers of the slab fixed, other atoms can be relaxed

    Parameters:

    slab:
        An aiida.orm.StructureData object. The slab that we want to adsorb some adsorbates on.

    adsSiteDictionary:
        An aiida.orm.Dict object. In which contains the adsorbates (for common adsorbates, I can assign the adsorption atom, geometry of the adsorbates and its adsorption angle in :code:`constants.py`) and adsorption sites, also I can add new adsorbates to the adsorbates library. (Do a very simple relax simulation of a molecule in the big cell, and just output the structure and assign the adsorption site. For mono-site adsorption, it is really easy, for bi-site adsorption, it is also easy since we can get the molecule to align; for tri-site or more-site adsorption, well, we can just put the molecule closer to the surface and let the program determine how does this molecule interact with the surface, for our usual research, I don't think this is necessary.)

        ```
        adsSiteDictionary = {
            'O': [site1, site2, site3] # each site is a 3x1 list [a, b, c]
            'F': [site1, site2, site3]
        }
        ```

    Return:
        An aiida.orm.StructureData object. With the adsorbates and modified constrains on all the atoms, ready for the submit functions.

    """

    # clean the input parameters from aiida.orm types to usual types
    slab = slab.get_pymatgen_structure()
    slab_tmp = deepcopy(slab)
    # asf = AdsorbateSiteFinder(slab_tmp, selective_dynamics = True)
    adsSiteDictionary = adsSiteDictionary.get_dict()
    # end of cleaning process

    for ads, siteList in adsSiteDictionary.items():
        adsorbate = getMoleculeByName(ads)
        if len(ads) == 1:
            ads_site = [0]
        else:
            ads_site = adsorbates[ads]['ads_site']
        if len(ads_site) == 1:
            for site in siteList:
                slab_tmp = hzd_add_adsMono(slab_tmp, molecule = adsorbate, ads_coord = [site], ads_site = ads_site)
        # elif len(ads_site) == 2:
        #     # treate it as mono-dent
        #     for site in siteList:
        #         for adsite in ads_site:
        #             slab_tmp = hzd_add_adsMono(slab_tmp, molecule = adsorbate, ads_coord = site, ads_site = ads_site)
        #     slab_tmp = hzd_add_adsBi(slab_tmp, molecule = adsorbate, ads_coord = siteList, ads_site = ads_site)

    return StructureData(pymatgen_structure=slab_tmp)

@calcfunction
def newStructure(structure, changeDict):

    """

    :code:`newStructure` will create a new structure by replacing the changeList atom to new type of atoms, very useful if there are different chemical states of the same atom in the structure.

    Parameters:

    structure:
        A aiida.orm.StructureData object.

    changeDict:
        A dictionary where the i-th atom's label will be replaced. e.g. 'Fe' to 'Fe2'

    Return:
        A new aiida.orm.StructureData file which stores the difference

    """

    # create the kind_names list
    kind_names = []
    structure_ase = structure.get_ase()
    for atom in structure_ase:
        kind_names.append(atom.symbol)

    for id, symbol in changeDict.items():
        kind_names[id] = symbol

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
def getMoleculeByName(str):

    """

    :code:`getMoleculeByName` can return the Molecule object of a specific adsorbates, could be very useful for the simulations.

    Parameters:

    str:
        A string. Which we will find whether it matches the adsorbates in the database

    Return:
        A Molecule Object, which can be used for adding the adsorbates.

    """

    from pymatgen.core import Element, Molecule
    from hzdplugins.aiidaplugins.constants import adsorbates

    if len(str) == 1: # means it is an atom
        if Element(str): # means str is an element
            return Molecule(species = (str), coords = [[0, 0, 0]])
        else:
            raise ValueError("You have entered a wrong string for an element, please try again.")
    else:
        if str in adsorbates.keys():
            return adsorbates[str]['mol']
        else:
            raise ValueError("Sorry, the adsorbates you are asking for haven't been added in the package yet.")

def setFixedCoords(slab, height):

    """

    :code:`setFixedCoords` is a function that can return a list that shows which atom we want to freeze during the simulation

    Parameters:

    slab:
        An aiida.orm.StructureData object.

    height:
        A float. which shows the boundary of the fixed atom. if the z-component of the position of the atom is below height, then we freeze it.

    Return:
        A Nx3 list which holds all the information about the fixed atoms, can be directly passed to the setting_dict.

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

    Parameters:

    slab:
        Pymatgen Structure object. A slab that we are interested in.

    molecule:
        The adsorbate that we want to add

    ads_coord:
        The position of adsorption site

    ads_site:
        The adsorption site of the molecule. If there is only 1 site in the molecule, that means it is mono-ads, if there are 2 sites in the molecule, then it can be mono-ads or bi-ads (needs to treat differently).

    Return:
        A modified slab. That all the adsorbates are in [True, True, True] (selective dynamics)

    """

    slab_tmp = deepcopy(slab)
    molecule_tmp = deepcopy(molecule)

    # check whether ads_coord and ads_site are compatible:
    if len(ads_coord) == len(ads_site):
        pass
    else:
        raise ValueError("Sorry, the adsorption site on the slab and on the molecule are different, they must have the same length.")

    if len(ads_site) == 1: # momo-ads
        coord_ads_site = molecule[ads_site[0]].coords
        vec = ads_coord[0] - coord_ads_site
        for site in molecule_tmp:
            site.coords = site.coords + vec
        # molecule_tmp = hzd_rotate(molecule_tmp, ads_site[0]) # make sure that the ads_site[0] are in the lowest point.
        molecule_tmp.add_site_property(
            "surface_properties", ['adsorbate']*molecule_tmp.num_sites
        )
        molecule_tmp.add_site_property(
            'selective dynamics', [[True, True, True]]*molecule_tmp.num_sites
        )
        for site in molecule_tmp:
            slab_tmp.append(species = site.specie, coords = site.coords, coords_are_cartesian=True, properties = site.properties)

    return slab_tmp
