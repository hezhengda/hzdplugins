Log file for hzdplugins
=======================

* v0.0.40: 2020.12.28
    * Add :code:`setSpinStructure` function in :code:`build.py`

* v0.0.39: 2020.12.19
    * remove all the :code:`results` parameter in the `submit.py`, and add :code:`group` parameter to make sure that we use the group functionality in Aiida, not something we create ourselves. Since the :code:`results` parameter is no longer needed, we also need to remove the :code:`saveResults`, :code:`readResults`, :code:`showResults`, :code:`getUnDoneTasks`, :code:`getUnFinishedTasks`, :code:`getUnConvergedTasks`, :code:`assignValue`, :code:`pkToUuidConverter` from `info.py` to `deprecated.py`
    * add `inputgenerator.py` in `hzdplugins.aiidaplugins` module in order to help me generate input data for the workchains.
        * For the outputdict(), do we need to specify which parameter we need to exclude? (Use a dictionay to store the parameter list) --> For the future.

* v0.0.38: 2020.12.12
    * add :code:`projwfcOriginalSubmit`, :code:`phOriginalSubmit` in `submit.py`
    * add :code:`getPdos` in :code:`info.py` to get the pdos for certain atom
    * add :code:`getDos` in :code:`info.py`
    * Since rotating the adsorbate by algorithm is very time-consuming and not efficient, so I just create the geometric structure of usual adsorbates in the :code:`hzdplugins.aiidaplugins.adsorbates` folder. All formats are :code:`mol`, which is very easy to be read by ase.

* v0.0.34-37: 2020.12.12
    * add some new functionalities in :code:`info.py` which can show the charge & magnetic moments, and also the total forces for each atomic step. Write :code:`get_TotalForces` and :code:`get_ChargeAndMagneticMoments`
    * add the :code:`info.py:getStructure` function which can help visualize the struture.
    * correct the :code:`info.py:checkDistance` function
    * modify the :code:`io.py` to deal with uuid rather than pk.
    * add :code:`math.py` in :code:`hzdplugins/structure` to store some mathematical functions
    * correct functions in `hzdplugins/aiidaplugins/submit.py`
    * add :code:`newStructure` function in `build.py` to deal with the system which we need to specify different chemical states or spin states. Very useful for complex cases.

* v0.0.33: 2020.12.05
    * give some options for the plotting in :code:`build.py` :code:`adsorptionSites` function.

* v0.0.32: 2020.12.05
    * delete :code:`print` in `build.py`, fixed.

* v0.0.31: 2020.12.05
    * Mising :code:`deepcopy` in `build.py`, fixed.

* v0.0.30: 2020.12.05
    * Missing :code:`Dict` in `build.py`, fixed.

* v0.0.29: 2020.12.05
    * Add a new module :code:`structure`, which can help me create bulk structures and slab surfaces more easily.
    * Make all the functions in :code:`structure` module and `submit.py` in aiida's standard, which means all the input variables are from :code:`aiida.orm` types and all the functions are decorated by :code:`calcfunction` in order to preserve the provenance.
    * Add a powerful function in `structure.py` called :code:`adsorptionSites`, in this function, if we set :code:`visualize=False`, then it only gives us the position of the adsorption sites; but if :code:`visualize=True`, then we will get a nice picture that shows all possible adsorption sites, and their ids, so we can visualize it and make plans for later calculations.

* v0.0.26-28: 2020.12.04
    * Assign :code:`uuid` instead of :code:`pk` for all the nodes, because :code:`uuid` is unique, but :code:`pk` is not, so this may cause some difficulties in the implementation. Also add a converter from :code:`pk` based results dictionary to :code:`uuid` based results dictionary.
    * Correct a few mistakes

* v0.0.25: 2020.12.04
    * modified the :code:`constants` file because I transfer the aiida from ubuntu to mac, so all the computer names have changed. e.g. from :code:`rwth-claix` to :code:`rwth-claix-mac`.

* v0.0.24: 2020.12.01
    * modified the :code:`assignValue` method in `info.py`, now it can process both 'pw.x' and 'projwfc' calculations.

* v0.0.23: 2020.11.30
    * Correct `constants.py`

* v0.0.22: 2020.11.30
    * Remember that all the pk have to be change to str, otherwise it will look ugly (e.g. 211.0 not 211)
    * change :code:`assignValue` function in `info.py`: Now the killed process can also be identified.

* v0.0.21: 2020.11.30
    * Add :code:`son_calc` to the constants in order to connect the "son" with the "father" CalcNode, because now not only we can search back, we can also search forward at the same time.

* v0.0.20: 2020.11.30
    * It turns out, the error is in constant.py, where :code:`max_wallclock_seconds` is mistaken to :code:`max_walllock_seconds` ...

* v0.0.19: 2020.11.30
    * add assignment for restart_builder.metadata.options when cluster_options is empty.

* v0.0.18: 2020.11.30
    * change :code:`assignValue` in `info.py`

* v0.0.17: 2020.11.29
    * add :code:`saveResults`, :code:`readResults` in `info.py` for save and read current results via json file.

* v0.0.15: 2020.11.29
    * change :code:`computer` in :code:`qePwContinueSubmit` in `submit.py`

* v0.0.14: 2020.11.29
    * change :code:`results_tmp[str(calc.pk)]['system'] = node.label` in :code:`qePwContinueSubmit` in :code:`submit.py`

* v0.0.13: 2020.11.29
    * resolve a mistake that Code is not imported.

* v0.0.12: 2020.11.29
    * resubmit the package, because some mistakes in :code:`qePwOriginalSubmit` haven't been resolved, but now they are removed.

* v0.0.11: 2020.11.29
    * I've changed the keys in results, now it has the following keys: ['system', 'uuid', 'comp_type', 'cluster', 'xc functional', 'exit_status', 'is_finished', 'is_finished_ok', 'E/eV', 'remove_remote_folder',  'previous_calc']
    * Add functions :code:`unDoneTasks`, :code:`unFinishedTasks`, :code:`unConvergedTasks`, which can be used in selecting the tasks that still needs attention.
    * put all the important information in `constants.py`
    * change the :code:`qePwOriginalSubmit` and :code:`qePwContinueSubmit` with the usage of `constants.py`, now the input becomes simpler.

* v0.0.5: 2020.11.28
    * I've learn that if you want to make a python module, then you need to add `__init__.py` file in the folder.

* v0.0.4: 2020.11.28
    * change the structure of the folder

* v0.0.3: 2020.11.28
    * add qePwOriginalSubmit and qePwContinueSubmit methods.
