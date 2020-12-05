# This module will be used in generating structures for Aiida simulation, it utilizes pymatgen, ase and catkit packages from other group, all the output structures should have Aiida StructureData type.
# In order to maintain provenance, we need to use all the types provided by the aiida.orm, and all functions should be decorated by @calcfunction decorator in order to be captured by the aiida program

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.path import Path
from copy import deepcopy

from hzdplugins.aiidaplugins.constants import color_dictionary

from aiida.orm import StructureData, Str, List, Bool, Dict
from aiida.engine import calcfunction

from ase.io import read

from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, get_rot

@calcfunction
def bulk(filename, supercell):

    """

    `bulk` function can help us create a Bulk from the file.

    Parameters:

    filename:
        An aiida.orm.Str object. Usually when we create the structural file, we do it from the strcutural file such as .xyz or .cif, etc.

    supercell:
        An aiida.orm.List object. A list that contains the dimension of supercell we would like.

    Return: A StructureData file that can be used directly in Aiida.

    """

    # transfer from aiida.orm types to the common python types
    filename = filename.value
    supercell_list = supercell.get_list()

    if len(filename)==0:
        raise(IOError("You didn't provide an input file."))

    if len(filename.split('.')[1])==0:
        raise(IOError("You didn't provide a correct type, please try to name your file .xyz, .cif or other data structure types."))

    type = filename.split('.')[1]

    structure = read(filename, format=type)

    return StructureData(ase=structure*supercell_list)

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
def adsorptionSites(slab, visualize, **kwargs):

    """

    `AdsorptionSites` can help us visualize the position and tag of each adsorption sites, then we can determine where we want to put the adsorbates.

    Parameters:

    slab:
        An aiida.orm.StructureData object. This is our slab, and we want to find how we can put the adsorbates on the slab.

    visualize:
        An aiida.orm.Bool object. A Boolean variable. If it is True, that means we want to output the figure of adsorption sites. If it is false, then we only want to return the list of adsorption sites.

    **kwargs:
        There are four variables which can be set for tuning the output figure
        * repeat (Int): how many unit cell do we want to show (how large the supercell we want to show)
        * decay (Float in [0,1]): the decay of alpha-value among different layers
        * scale (Float): radius scaling for sites (larger scale, larger circles)
        * window (Float): window for setting the axes limit

    Return:

    Matplotlib figure:
        A matplotlib figure, which represents the distribution of the adsorption sites and its corresponding id.

    Adsorption site dictionary:
        An aiida.orm.Dict object. Contains the information of all the positions of adsorption sites, which corresponse to the matplotlib figure.

    """

    # the inspiration for this function was from pymatgen.analysis.adsorption.AdsorbateSiteFinder.plot_slab() function, which is really intuitive way of showing the structure and the adsorption sites.
    # since this function does not create any useful data, so it doesn't be decorated with @calfunction decorator

    # get the structure and adsorption sites
    slab = slab.get_pymatgen_structure()
    visualize = visualize.value
    # end of the conversion

    asf = AdsorbateSiteFinder(slab, selective_dynamics = False)
    adsorption_sites = asf.find_adsorption_sites(distance = 1.2)

    dictGenerator = Dict()
    dictGenerator.set_dict(adsorption_sites)

    if visualize == False:
        pass
    else:
        # start plotting the figure, this part of the code was largely adopted from pymatgen github repository

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

    return dictGenerator
