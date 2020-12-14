Tutorial
========

In this tutorial, we will introduce how to use my plugins for accelerating the computational research.

My plugins :code:`hzdplugins` have the following module structures:

* aiidaplugins
    * :code:`submit` module: for submitting the computational jobs
    * :code:`io` module: for managing the contents on remote server
    * :code:`info` module: for getting information from structures

* structure
    * :code:`build` module: for building the structure that we want to calculate.

So in the following tutorial, I will show a very simple script that can calculate the OH adsorption on Pt(111) surface.

First you need to setup your Aiida environment, for how to do that I have uploaded 2 articles:

* For Linux: https://zhuanlan.zhihu.com/p/295175866
* For MacOS: https://zhuanlan.zhihu.com/p/300490416

.. code-block::

    # import the important functions
    %aiida # start the basic aiida environment
    from hzdplugins.structure.build import bulkFromString, millerSurfaces, adsorptionSites, visualizeSlab, addAdsorbates, setFixedCoords
    from hzdplugins.aiidaplugins.info import getStructure
    from hzdplugins.aiidaplugins.submit import qePwOriginalSubmit
    from aiida.orm import Int, List, Dict, Str, Bool, Float, Data

    # creating the Pt_bulk
    Pt_bulk = bulkFromString(bulkStr = Str('Pt'), crystal_structure=Str('fcc'), a = Float(3.98), cubic=Bool(True), supercell = List(list=[1, 1, 1]))

    # creating the slab
    slab_list = millerSurfaces(bulk = Pt_bulk, miller_index = List(list=[1, 1, 1]), layers = Int(4), vacuum = Int(8))
    slab_list = slab_list.get_list()
    surface = load_node(slab_list[0])

    # find adsorption sites
    adsorption_sites = adsorptionSites(slab = surface)
    adsorption_sites = adsorption_sites.get_dict()

    # create adsorption configuration
    adsSiteDictionary = {
        'OH': [adsorption_sites['ontop'][0]]
    }

    # create slab with the surface and adsorbates
    surfSlabWithAds = addAdsorbates(slab = surface,
                                    adsSiteDictionary=Dict(dict=adsSiteDictionary))

    # assign fixed_coords (or selective dynamics in VASP term)
    fixed_coords = setFixedCoords(slab=surfSlabWithAds, height=14)

    # settings for the pw.x calculation
    settings_dict = {
    'fixed_coords': fixed_coords,
    'cmdline': ['-nk', '4'],
    'additional_retrieve_list': ['aiida.out']
    }

    # start the pw.x original calculation (no previous calculations are needed)
    results = {}
    results, uuid_new = qePwOriginalSubmit(
        results = results,
        codename = 'pw.x-6.6-juwels@juwels-mac',
        structure = surfSlabWithAds,
        kpoints = [[4, 4, 1]],
        pseudo_family = 'sssp-precision',
        pseudo_dict = {},
        metadata = {
            'label': 'Pt111-OH',
            'description': 'test whether my package is good.'
        },
        add_parameters = {
            'SYSTEM': {
                'input_dft': 'PBESOL'
            },
            'ELECTRONS': {
                'mixing_mode': 'local-TF',
                'mixing_beta': 0.3,
                'mixing_ndim': 10
            }
        },
        del_parameters = {},
        cluster_options = {},
        settings_dict = settings_dict
    )

    # then you can check the results by using
    # node = load_node(uuid=uuid_new)