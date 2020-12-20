# This file contains all the functions that is no longer needed by the package, so everything should in comment.

# From info.py

# def saveResults(results, filename):
#     """
#
#     :code:`saveResults` can be used for saving results
#
#     :param results: The results dictionary that contains all relevant computational information.
#     :type results: python dictionary object
#
#     :param filename: The name of the file that you want to store in.
#     :type filename: python string object
#
#     :returns: A dictionary that has modified and assigned values
#
#     """
#
#     # dump results file
#     with open(filename, 'w') as json_file:
#         json.dump(results, json_file)
#
#     return 'Your results have been successfully saved.'

# def readResults(filename):
#     """
#
#     :code:`readResults` can be used for reading results from json file
#
#     :param filename: The name of the file that you store all the information in.
#     :type filename: python string object
#
#     :returns: A dictionary that has modified and assigned values
#
#     """
#
#     with open(filename, 'r') as json_file:
#         results = json.load(json_file)
#
#     return results

# def showResults(results):
#     """
#
#     :code:`showResults` shows results in pandas form.
#
#     :param results: The results dictionary that contains all relevant computational information.
#     :type results: python dictionary object
#
#     :returns: A pandas object that we can use to display the dataframe
#
#     """
#
#     df = pd.DataFrame.from_dict(results,
#                                 orient='index',
#                                 columns=results_keys_set)
#     pd.set_option("display.max_rows", None, "display.max_columns", None)
#     return df

# def getUnDoneTasks(results):
#     """
#
#     :code:`getUnDoneTasks` shows all the unDone computational tasks, which means the exit_status is neither '0' or
#     '501', or None.
#
#     :param results: The results dictionary that contains all relevant computational information.
#     :type results: python dictionary object
#
#     :returns: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.
#
#     """
#
#     subresults = {}
#
#     for key, value in results.items():
#         if 'exit_status' in value.keys():
#             if value['exit_status'] != '0' and value['exit_status'] != '501':
#                 subresults[key] = value
#         else:
#             node = load_node(uuid=key)
#             if node.exit_status != '0' and node.exit_status != '501':
#                 subresults[key] = value
#
#     return subresults

# def getUnFinishedTasks(results):
#     """
#
#     :code:`getUnFinishedTasks` show all the unFinished computational tasks, which means the is_finished tag is False
#
#     :param results: The results dictionary that contains all relevant computational information.
#                     The keys are uuid of each nodes
#     :type results: python dictionary object
#
#     :returns: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.
#
#     """
#
#     subresults = {}
#
#     for key, value in results.items():
#         node = load_node(uuid=key)
#         if not node.is_finished:
#             subresults[key] = value
#
#     return subresults

# def getUnConvergedTasks(results):
#     """
#
#     :code:`getUnConvergedTasks` show all the finished by not converged computational tasks, which means the
#     `is_finished` tag is True, but `is_finished_ok` is False, although when `exit_status == '501'`, it is still ok.
#
#     :param results: The results dictionary that contains all relevant computational information.
#     :type results: python dictionary object
#
#     :returns: A dictionary that contain all the unfinished and unconverged tasks, in total 'unDown' tasks.
#
#     """
#
#     subresults = {}
#
#     for key, value in results.items():
#         node = load_node(uuid=key)
#         if node.is_finished:
#             if not node.is_finished_ok:
#                 if node.exit_status != 501:
#                     subresults[key] = value
#
#     return subresults

# def assignValue(results):
#     """
#
#     :code:`assignValue` will assign the current status of simulaton to results. The function will do two things:
#     (1)  clean the current results dictionary, remove any key that does not belong to the new key set. (2) Add the
#     values in the new set.
#
#     :param results: The results dictionary that contains all relevant computational information.
#     :type results: python dictionary object
#
#     :returns: A dictionary that has modified and assigned values
#
#     """
#
#     results_tmp = deepcopy(results)
#
#     for uuid_node, value in results.items():
#         value_tmp = deepcopy(value)
#         node = load_node(uuid=uuid_node)
#         # clean the results
#         for key in value_tmp:
#             if not (key in results_keys_set):
#                 results_tmp[uuid_node].pop(key)
#         # assign the value
#         results_tmp[uuid_node]['system'] = node.label
#         results_tmp[uuid_node]['cluster'] = node.computer.label
#
#         keys = [key.lower() for key in
#                 node.inputs.parameters.get_dict().keys()]  # get all the keys in lower case, easy for comparison
#
#         if 'control' in keys:  # it is pw.x calculation
#             results_tmp[uuid_node]['comp_type'] = node.inputs.parameters.get_dict()['CONTROL']['calculation']
#             results_tmp[uuid_node]['xc functional'] = node.inputs.parameters.get_dict()['SYSTEM']['input_dft']
#             compcode = 'pw'
#
#         if 'projwfc' in keys:  # it is projwfc calculation
#             results_tmp[uuid_node]['comp_type'] = 'PDOS'
#             compcode = 'projwfc'
#
#         if 'inputph' in keys:
#             results_tmp[uuid_node]['comp_type'] = 'PH'
#             compcode = 'ph'
#
#         if node.is_finished_ok or (node.exit_status == 0) or (node.exit_status == 501):
#             if compcode == 'pw':
#                 compType = node.inputs.parameters['CONTROL']['calculation']
#                 if compType == 'relax' or compType == 'vc-relax':
#                     results_tmp[uuid_node]['E/eV'] = node.res.energy
#                 results_tmp[uuid_node]['is_finished'] = node.is_finished
#                 results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
#                 results_tmp[uuid_node]['exit_status'] = str(node.exit_status)
#             if compcode == 'projwfc':
#                 results_tmp[uuid_node]['is_finished'] = node.is_finished
#                 results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
#                 results_tmp[uuid_node]['exit_status'] = str(node.exit_status)
#             if compcode == 'ph':
#                 results_tmp[uuid_node]['is_finished'] = node.is_finished
#                 results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
#                 results_tmp[uuid_node]['exit_status'] = str(node.exit_status)
#         else:
#             results_tmp[uuid_node]['E/eV'] = None
#             results_tmp[uuid_node]['is_finished'] = node.is_finished
#             results_tmp[uuid_node]['is_finished_ok'] = node.is_finished_ok
#             if node.is_killed:
#                 results_tmp[uuid_node]['exit_status'] = 'killed'
#             elif node.is_excepted:
#                 results_tmp[uuid_node]['exit_status'] = 'excepted'
#             else:
#                 results_tmp[uuid_node]['exit_status'] = str(node.exit_status)
#
#     return results_tmp

# def pkToUuidConverter(results_pk):
#     """
#
#     :code:`pkToUuidConverter` can convert the :code:`pk` based results dictionary to :code:`uuid` based dictionary.
#
#     :param results_pk: The results dictionary that has pk as its key
#     :type results_pk: python dictionary object
#
#     :returns: A dictionary that has uuid as its key
#
#     """
#
#     results_tmp = {}  # our empty and temporary dictionary for storing the converted dictionary
#
#     for key, value in results_pk.items():
#         # get the node by uuid, pk may not work (different profile)
#         uuid_node = value['uuid']
#         qb = QueryBuilder()
#         qb.append(Node, filters={'uuid': {'==': uuid_node}})
#         # node = qb.all()[0][0]  # because uuid is unique, it can only return to one node
#
#         # reconstruct results_tmp, add the checking mechanism
#         results_tmp[uuid_node] = {}
#         keyset = value.keys()
#         if 'system' in keyset:
#             results_tmp[uuid_node]['system'] = value['system']
#         else:
#             results_tmp[uuid_node]['system'] = None
#
#         if 'comp_type' in keyset:
#             results_tmp[uuid_node]['comp_type'] = value['comp_type']
#         else:
#             results_tmp[uuid_node]['comp_type'] = None
#
#         if 'cluster' in keyset:
#             results_tmp[uuid_node]['cluster'] = value['cluster']
#         else:
#             results_tmp[uuid_node]['cluster'] = None
#
#         if 'xc functional' in keyset:
#             results_tmp[uuid_node]['xc functional'] = value['xc functional']
#         else:
#             results_tmp[uuid_node]['xc functional'] = None
#
#         if 'exit_status' in keyset:
#             results_tmp[uuid_node]['exit_status'] = value['exit_status']
#         else:
#             results_tmp[uuid_node]['exit_status'] = None
#
#         if 'is_finished' in keyset:
#             results_tmp[uuid_node]['is_finished'] = value['is_finished']
#         else:
#             results_tmp[uuid_node]['is_finished'] = None
#
#         if 'is_finished_ok' in keyset:
#             results_tmp[uuid_node]['is_finished_ok'] = value['is_finished_ok']
#         else:
#             results_tmp[uuid_node]['is_finished_ok'] = None
#
#         if 'E/eV' in keyset:
#             results_tmp[uuid_node]['E/eV'] = value['E/eV']
#         else:
#             results_tmp[uuid_node]['E/eV'] = None
#
#         if 'remove_remote_folder' in keyset:
#             results_tmp[uuid_node]['remove_remote_folder'] = value['remove_remote_folder']
#         else:
#             results_tmp[uuid_node]['remove_remote_folder'] = None
#
#         # for previous_calc and son_calc, we need to change them to uuid
#         if 'previous_calc' in keyset:
#             results_tmp[uuid_node]['previous_calc'] = results_pk[value['previous_calc']]['uuid']
#         else:
#             results_tmp[uuid_node]['previous_calc'] = None
#
#         if ('son_calc' in keyset) and (value['son_calc'] in results_pk.keys()):
#             results_tmp[uuid_node]['son_calc'] = results_pk[value['son_calc']]['uuid']
#         else:
#             results_tmp[uuid_node]['son_calc'] = None
#
#     return results_tmp

# submit.py

# # results
# results_tmp[str(calc.uuid)] = {}
# results_tmp[str(calc.uuid)]['system'] = restart_builder.metadata.label
# results_tmp[str(calc.uuid)]['comp_type'] = parameters_default['CONTROL']['calculation']
# results_tmp[str(calc.uuid)]['E/eV'] = None
# results_tmp[str(calc.uuid)]['remove_remote_folder'] = False
# results_tmp[str(calc.uuid)]['cluster'] = computer
# results_tmp[str(calc.uuid)]['xc functional'] = parameters_default['SYSTEM']['input_dft']
# results_tmp[str(calc.uuid)]['exit_status'] = None
# results_tmp[str(calc.uuid)]['is_finished'] = None
# results_tmp[str(calc.uuid)]['is_finished_ok'] = None
# results_tmp[str(calc.uuid)]['previous_calc'] = uuid  # pk is the previous calculation
# results_tmp[str(calc.uuid)]['son_calc'] = None  # right now there is no son_node
#
# # change son_calc with previous simulation
# results_tmp[uuid]['son_calc'] = calc.uuid

