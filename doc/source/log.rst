Log file for hzdplugins
=======================

* v0.0.34-37: 2020.12.12
    * add some new functionalities in `info.py` which can show the charge & magnetic moments, and also the total forces for each atomic step. Write `get_TotalForces` and `get_ChargeAndMagneticMoments`
    * add the `info.py:getStructure` function which can help visualize the struture.
    * correct the `info.py:checkDistance` function
    * modify the `io.py` to deal with uuid rather than pk.
    * add `math.py` in `hzdplugins/structure` to store some mathematical functions
    * correct functions in `hzdplugins/aiidaplugins/submit.py`
    * add `newStructure` function in `build.py` to deal with the system which we need to specify different chemical states or spin states. Very useful for complex cases.

* v0.0.33: 2020.12.05
    * give some options for the plotting in `build.py` `adsorptionSites` function.

* v0.0.32: 2020.12.05
    * delete `print` in `build.py`, fixed.

* v0.0.31: 2020.12.05
    * Mising `deepcopy` in `build.py`, fixed.

* v0.0.30: 2020.12.05
    * Missing `Dict` in `build.py`, fixed.

* v0.0.29: 2020.12.05
    * Add a new module `structure`, which can help me create bulk structures and slab surfaces more easily.
    * Make all the functions in `structure` module and `submit.py` in aiida's standard, which means all the input variables are from `aiida.orm` types and all the functions are decorated by `calcfunction` in order to preserve the provenance.
    * Add a powerful function in `structure.py` called `adsorptionSites`, in this function, if we set `visualize=False`, then it only gives us the position of the adsorption sites; but if `visualize=True`, then we will get a nice picture that shows all possible adsorption sites, and their ids, so we can visualize it and make plans for later calculations.

* v0.0.26-28: 2020.12.04
    * Assign `uuid` instead of `pk` for all the nodes, because `uuid` is unique, but `pk` is not, so this may cause some difficulties in the implementation. Also add a converter from `pk` based results dictionary to `uuid` based results dictionary.
    * Correct a few mistakes

* v0.0.25: 2020.12.04
    * modified the `constants` file because I transfer the aiida from ubuntu to mac, so all the computer names have changed. e.g. from `rwth-claix` to `rwth-claix-mac`.

* v0.0.24: 2020.12.01
    * modified the `assignValue` method in `info.py`, now it can process both 'pw.x' and 'projwfc' calculations.

* v0.0.23: 2020.11.30
    * Correct `constants.py`

* v0.0.22: 2020.11.30
    * Remember that all the pk have to be change to str, otherwise it will look ugly (e.g. 211.0 not 211)
    * change `assignValue` function in `info.py`: Now the killed process can also be identified.

* v0.0.21: 2020.11.30
    * Add `son_calc` to the constants in order to connect the "son" with the "father" CalcNode, because now not only we can search back, we can also search forward at the same time.

* v0.0.20: 2020.11.30
    * It turns out, the error is in constant.py, where `max_wallclock_seconds` is mistaken to `max_walllock_seconds` ...

* v0.0.19: 2020.11.30
    * add assignment for restart_builder.metadata.options when cluster_options is empty.

* v0.0.18: 2020.11.30
    * change `assignValue` in `info.py`

* v0.0.17: 2020.11.29
    * add `saveResults`, `readResults` in `info.py` for save and read current results via json file.

* v0.0.15: 2020.11.29
    * change `computer` in `qePwContinueSubmit` in `submit.py`

* v0.0.14: 2020.11.29
    * change `results_tmp[str(calc.pk)]['system'] = node.label` in `qePwContinueSubmit` in `submit.py`

* v0.0.13: 2020.11.29
    * resolve a mistake that Code is not imported.

* v0.0.12: 2020.11.29
    * resubmit the package, because some mistakes in `qePwOriginalSubmit` haven't been resolved, but now they are removed.

* v0.0.11: 2020.11.29
    * I've changed the keys in results, now it has the following keys: ['system', 'uuid', 'comp_type', 'cluster', 'xc functional', 'exit_status', 'is_finished', 'is_finished_ok', 'E/eV', 'remove_remote_folder',  'previous_calc']
    * Add functions `unDoneTasks`, `unFinishedTasks`, `unConvergedTasks`, which can be used in selecting the tasks that still needs attention.
    * put all the important information in `constants.py`
    * change the `qePwOriginalSubmit` and `qePwContinueSubmit` with the usage of `constants.py`, now the input becomes simpler.

* v0.0.5: 2020.11.28
    * I've learn that if you want to make a python module, then you need to add `__init__.py` file in the folder.

* v0.0.4: 2020.11.28
    * change the structure of the folder

* v0.0.3: 2020.11.28
    * add qePwOriginalSubmit and qePwContinueSubmit methods.
