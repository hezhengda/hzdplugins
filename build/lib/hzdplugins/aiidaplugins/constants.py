# all the keys in the results dictionary
results_keys_set = ['system', 'uuid', 'comp_type', 'cluster', 'xc functional', 'exit_status', 'is_finished', 'is_finished_ok', 'E/eV', 'remove_remote_folder',  'previous_calc']

# dictionary for the slurm system on different clusters
slurm_options = {}
# 'claix'
slurm_options['rwth-claix'] = {
    'resources': {'num_machines', 2}, # on rwth_claix, each node has 48 cores.
    'max_wallclock_seconds': 86400,
    'account': 'jara0037',
    'scheduler_stderr': 'stderr',
    'scheduler_stdout': 'stdout',
    'queue_name': 'c18m'
}

# 'juwels'
slurm_options['juwels'] = {
    'resources': {'num_machines', 2}, # on juwels, each node has 48 cores.
    'max_wallclock_seconds': 86400,
    'account': 'fzj-mac',
    'scheduler_stderr': 'stderr',
    'scheduler_stdout': 'stdout',
    'queue_name': 'batch'
}

# 'jusuf'
slurm_options['jusuf'] = {
    'resources': {'num_machines': 1}, # on jusuf, each node has 128 cores
    'max_walllock_seconds': 86400,
    'account': 'jiek61',
    'scheduler_stderr': 'stderr',
    'scheduler_stdout': 'stdout',
    'queue_name': 'batch'
}

# 'jureca_booster'
slurm_options['jureca-booster'] = {
    'resources': {'num_machines': 2}, # on jureca_booster, each node has 68 cores
    'max_walllock_seconds': 86400,
    'account': 'jiek61',
    'scheduler_stderr': 'stderr',
    'scheduler_stdout': 'stdout',
    'queue_name': 'booster'
}
