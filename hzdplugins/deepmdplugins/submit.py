import json
class deepmdInputGenerator():
    """
    deepmdInputGenerator will help us to generate the input file for the deepmd-kit
    """

    def __init__(self, model={}, learning_rate={}, loss={}, training={}):
        """
        In the input file of deepmd-kit, there are four major components: (1) model (2) learning rate (3) loss function (4) training

        :param model: represent the model of the system, including the descriptor, defaults to {}
        :type model: python dictionary, optional
        :param learning_rate: the learning rate of neural network, defaults to {}
        :type learning_rate: python dictionary, optional
        :param loss: the parameters for the loss functions, defaults to {}
        :type loss: python dictionary, optional
        :param training: the parameters for the training, defaults to {}
        :type training: python dictionary, optional
        """
        self.model = model
        self.learning_rate = learning_rate
        self.loss = loss 
        self.training = training

    def setDescriptor_se_a(self, sel, rcut=6.0, rcut_smth=0.5, neuron=[10, 20, 40], axis_neuron=4, activation_function='tanh', resnet_dt=False, type_one_side=False, precision='float64', trainable=True, seed=1, exclude_types=[], set_davg_zero=False):
        """The `se_a` is (I think) the most popular descriptor, so that I should only used it in here.

        :param sel: a list of integers. The length of the list should be the same as the number of atom types in the system. `sel[i]` gives the selected number of typi-i neighbors. `sel[i]` is recommended to be larger than the maximally possible number of type-i neighbors in the cut-ff radius. 
        :type sel: python list 
        :param rcut: the cut-off radius, defaults to 6.0
        :type rcut: python float, optional
        :param rcut_smth: the 1/r term is smoothed from `rcut` to `rcut_smth`, defaults to 0.5
        :type rcut_smth: python float, optional
        :param neuron: number of neurons in each hidden layers of the embedding net. When two layers are of the same size or one layer is twice as large as the previous layer, a skip connection is built. defaults to [10, 20, 40]. This list could influence the accuracy of our model.
        :type neuron: python list, optional
        :param axis_neuron: size of the submatrix of G (embedding matrix), defaults to 4
        :type axis_neuron: python int, optional
        :param activation_function: the activation function in the embedding net, supported activation functions are: `relu`, `relu6`, `softplus`, `sigmoid`, `tanh`, `gelu`, defaults to 'tanh'
        :type activation_function: python str, optional
        :param seed: random seed for parameter initialization. Usually it is set to be 1.
        :type seed: python int.
        :param resnet_dt: whether to use a `timestep` in the skip connection, defaults to False
        :type resnet_dt: bool, optional
        :param type_one_side: try to build N_types embedding nets. Otherwise, building N_types^2 embedding nets, defaults to False
        :type type_one_side: bool, optional
        :param precision: the precision of the float numbers, defaults to 'float64'
        :type precision: str, optional
        :param trainable: if the parameters in the embedding net is trainable, defaults to True
        :type trainable: bool, optional
        :param exclude_types: the excluded types, defaults to []
        :type exclude_types: list, optional
        :param set_davg_zero: set the normalization average to zero. This option should be set when `atom_ener` in the energy fitting is used, defaults to False
        :type set_davg_zero: bool, optional
        :return: the descriptor (feature vector) for each atom 
        :rtype: a python list? 
        """
        tmp_descriptor = {}
        tmp_descriptor['type'] = 'se_a'
        tmp_descriptor['sel'] = sel 
        tmp_descriptor['rcut'] = rcut
        tmp_descriptor['rcut_smth'] = rcut_smth
        tmp_descriptor['neuron'] = neuron 
        tmp_descriptor['axis_neuron'] = axis_neuron
        tmp_descriptor['activation_function'] = activation_function
        tmp_descriptor['resnet_dt'] = resnet_dt
        tmp_descriptor['type_one_side'] = type_one_side
        tmp_descriptor['precision'] = precision
        tmp_descriptor['trainable'] = trainable
        tmp_descriptor['seed'] = seed 
        tmp_descriptor['exclude_types'] = exclude_types
        tmp_descriptor['set_davg_zero'] = set_davg_zero
        return tmp_descriptor

    def setFittingNet(self, fntype = 'ener', numb_fparam = 0, numb_aparm = 0, neuron = [120, 120, 120], activation_function = 'tanh', precision = 'float64', resnet_dt = True, trainable = True, rcond = 0.001, seed = 1, atom_ener = []):
        """The fitting of physical properties

        :param fntype: the types of the fitting, usually we like to fit the energy (potential energy surface), but deep potential can also fit `dipole`, `polar` and also `global_polar`, defaults to 'ener'.
        :type fntype: python str, optional
        :param numb_fparam: the dimension of the frame parameter. If set to > 0, file `fparam.npy` should be included to provided the input fparams, defaults to 0
        :type numb_fparam: int, optional
        :param numb_aparm: the dimension of atomic parameter. If set to > 0, the `aparam.npy` should be included to provided the input aparams, defaults to 0
        :type numb_aparm: int, optional
        :param neuron: the number of neurons in each hidden layers of the fitting net, defaults to [120, 120, 120]
        :type neuron: list, optional
        :param activation_function: the activation function used in the neural network, defaults to 'tanh'
        :type activation_function: str, optional
        :param precision: the precision for float numbers, defaults to 'float64'
        :type precision: str, optional
        :param resnet_dt: whether to use a `timestep` in the skip connection, defaults to True
        :type resnet_dt: bool, optional
        :param trainable: whether the parameter in the fitting net is trainable, defaults to True
        :type trainable: bool, optional
        :param rcond: the condition number used to determine the initial energy shift for each type of atoms, defaults to 0.001
        :type rcond: float, optional
        :param seed: random seed for parameter initialization of the fitting net, defaults to 1
        :type seed: int, optional
        :param atom_ener: specify the atomic energy in vacuum for each type, defaults to []
        :type atom_ener: list, optional
        :return: a dictionary that can be prepared in the mode stages 
        :rtype: python dictionary 
        """
        tmp_fitting_net = {}
        tmp_fitting_net['type'] =  fntype
        tmp_fitting_net['numb_fparam'] = numb_fparam
        tmp_fitting_net['numb_aparam'] = numb_aparm
        tmp_fitting_net['neuron'] = neuron 
        tmp_fitting_net['activation_function'] = activation_function
        tmp_fitting_net['precision'] = precision
        tmp_fitting_net['resnet_dt'] = resnet_dt
        tmp_fitting_net['trainable'] = trainable
        tmp_fitting_net['rcond'] = rcond 
        tmp_fitting_net['seed'] = seed 
        tmp_fitting_net['atom_ener'] = atom_ener
        return tmp_fitting_net
        
    def setModel(self, type_map, descriptor, fitting_net):
        """This function can help us set the parameters for the model part

        :param type_map: a list of elements investigated in here
        :type type_map: python list
        :param descriptor: the parameters for the descriptor, this can be generated by using self.setDescriptor_se_a() method
        :type descriptor: python dictionary
        :param fitting_net: the parameters for the fitting new, this can be generated by using self.setFittingNet() method
        :type fitting_net: python dictionary
        """
        tmp_model = {}
        tmp_model['type_map'] = type_map
        tmp_model['descriptor'] = descriptor
        tmp_model['fitting_net'] = fitting_net
        self.model = tmp_model
    
    def setLearningRate(self, lr_type='exp', decay_steps=5000, start_lr=0.001, stop_lr=1e-08):
        tmp_learningrate = {}
        tmp_learningrate['type'] = lr_type
        tmp_learningrate['decay_steps'] = decay_steps
        tmp_learningrate['start_lr'] = start_lr
        tmp_learningrate['stop_lr'] = stop_lr
        self.learning_rate = tmp_learningrate

    def setLoss(self, start_pref_e=0.02, limit_pref_e=1.0, start_pref_f=1000, limit_pref_f=1.0, start_pref_v=0.0, limit_pref_v=0.0):
        tmp_loss = {}
        tmp_loss['start_pref_e'] = start_pref_e 
        tmp_loss['limit_pref_e'] = start_pref_e
        tmp_loss['start_pref_f'] = start_pref_f
        tmp_loss['limit_pref_f'] = limit_pref_f
        tmp_loss['start_pref_v'] = start_pref_v
        tmp_loss['limit_pref_v'] = limit_pref_v
        self.loss = tmp_loss

    def setTrainingData(self, systems=['./'], set_prefix='set', batch_size='auto', auto_prob='prob_sys_size', sys_probs=None):
        """This function is used to set the training_data section in the training section

        :param systems: The data systems for training. This key can be provided with a list that specifies the systems, or be provided with a string by which the prefix of all systems are given and the list of the systems is automatically generated. 
        :type systems: python list 
        :param set_prefix: the prefix of sets in the systems, defaults to 'set'
        :type set_prefix: str, optional
        :param batch_size: automatically determines the batch size so that the batch_size times the number of atoms in the system is no less than 32., defaults to 'auto'
        :type batch_size: str, optional
        :param auto_prob: defaults to 'prob_sys_size', the probability of a system is proportional to the number of batches in the system
        :type auto_prob: str, optional
        :param sys_probs: defaults to None
        :type sys_probs: None, optional
        :return: the dictionary of `training_data` 
        :rtype: python dictionary 
        """
        tmp_training_data = {}
        tmp_training_data['systems'] = systems 
        tmp_training_data['set_prefix'] = set_prefix 
        tmp_training_data['batch_size'] = batch_size 
        tmp_training_data['auto_prob'] = auto_prob
        tmp_training_data['sys_probs'] = sys_probs
        return tmp_training_data
        
    def setValidationData(self, systems=['./'], set_prefix='set', batch_size='auto', auto_prob='prob_sys_size', sys_probs=None, numb_batch=1):
        """This function is used to set the validation_data section in the validation section

        :param systems: The data systems for validation. This key can be provided with a list that specifies the systems, or be provided with a string by which the prefix of all systems are given and the list of the systems is automatically generated. 
        :type systems: python list 
        :param set_prefix: the prefix of sets in the systems, defaults to 'set'
        :type set_prefix: str, optional
        :param batch_size: automatically determines the batch size so that the batch_size times the number of atoms in the system is no less than 32., defaults to 'auto'
        :type batch_size: str, optional
        :param auto_prob: defaults to 'prob_sys_size', the probability of a system is proportional to the number of batches in the system
        :type auto_prob: str, optional
        :param sys_probs: defaults to None
        :type sys_probs: None, optional
        :param numb_batch: an integer that specifies the number of systems to be sampled for each validation period.
        :type numb_batch: python int
        :return: the dictionary of `validation_data` 
        :rtype: python dictionary 
        """
        tmp_validation_data = {}
        tmp_validation_data['systems'] = systems 
        tmp_validation_data['set_prefix'] = set_prefix 
        tmp_validation_data['batch_size'] = batch_size 
        tmp_validation_data['auto_prob'] = auto_prob
        tmp_validation_data['sys_probs'] = sys_probs
        tmp_validation_data['numb_batch'] = numb_batch
        return tmp_validation_data

    def setTraining(self, training_data, validation_data, stop_batch=1000000, seed=1, \
                    numb_test=1, disp_file='lcurve.out', disp_freq='500', save_freq='1000', save_ckpt='model.ckpt', \
                    disp_training=True, time_training=True, profiling=True, profiling_file='timeline.json'):
        """This function can help us generate the parameters for the training section

        :param training_data: parameters about the training_data, this can be generated from self.setTrainingData() method 
        :type training_data: python dictionary 
        :param validation_data: parameters about the validation_data, this can be generated from self.setValidationData() method 
        :type validation_data: python dictionary
        :param stop_batch: the number of batches that we want to train the model, defaults to 1000000
        :type stop_batch: int, optional
        :param seed: the seed for generating the random number, defaults to 1
        :type seed: int, optional
        :param disp_file: the output file for the error, defaults to 'lcurve.out'
        :type disp_file: str, optional
        :param disp_freq: output every `disp_freq` steps, defaults to '500'
        :type disp_freq: str, optional
        :param save_freq: steps for saving the training data, defaults to '1000'
        :type save_freq: str, optional
        :param save_ckpt: the file for storing the training stage, defaults to 'model.ckpt'
        :type save_ckpt: str, optional
        :param time_training: output the training time for each step, defaults to True
        :type time_training: bool, optional
        :param profiling: output the profile (cpu / gpu usages for training), defaults to True
        :type profiling: bool, optional
        :param profiling_file: filename of the profiling data, defaults to 'timeline.json'
        :type profiling_file: str, optional
        """
        tmp_training = {}
        tmp_training['training_data'] = training_data 
        tmp_training['validation_data'] = validation_data 
        tmp_training['stop_batch'] = stop_batch
        tmp_training['seed'] = seed 
        tmp_training['numb_test'] = numb_test
        tmp_training['disp_file'] = disp_file
        tmp_training['disp_freq'] = disp_freq
        tmp_training['save_freq'] = save_freq
        tmp_training['save_ckpt'] = save_ckpt
        tmp_training['disp_training'] = disp_training
        tmp_training['time_training'] = time_training
        tmp_training['profiling'] = profiling
        tmp_training['profiling_file'] = profiling_file
        self.training = tmp_training
 
    def output(self):
        """
        This function can output the current stages of the input file for deepmd-kit, very useful for checking which parameter is still missing
        """
        tmp_dict = {}
        tmp_dict['model'] = self.model
        tmp_dict['learning_rate'] = self.learning_rate
        tmp_dict['loss'] = self.loss
        tmp_dict['training'] = self.training 
        return tmp_dict

    def export(self, filename):
        """`export` function can help us generate the json file

        :param filename: the name of the json file 
        :type filename: python string 
        """

        # check
        if len(list(self.model.keys())) == 0:
            print('There are no paramters in model')
        elif len(list(self.learning_rate.keys())):
            print('There are no parameters in learning_rate')
        elif len(list(self.loss.keys())):
            print('There are no parameters in loss') 
        elif len(list(self.training.keys())):
            print('There are no parameters in training')
        else:
            print('everything is fine, generating the input json file now') 
        
        tmp_dict = {
            'model': self.model,
            'learning_rate': self.learning_rate,
            'loss': self.loss,
            'training': self.training
        }
        f = open(filename,'w+')
        json.dump(tmp_dict, f)

class dpgenInputGenerator():
    """
    dpgenInputGenerator will help us to generate the input json file for the dpgen package.
    """

class lammpsInputGenerator():
    """
    `lammpsInputGenerator` will help us to generate the input file for LAMMPS (if we want to combine the deepmd with lammps)

    There are four parts of the input file for LAMMPS:
    * Initilization
    * System definition
    * Simulation settings
    * Run a simulation
    """

    def __init__(self, sys_name, data_file, dp_graph_file, temperature, timestep, run_num):
        self.dict = {}
        self.dict['boundary'] = 'p p p'
        self.dict['units'] = 'metal'
        self.dict['atom_style'] = 'atomic'
        self.dict['neighbor'] = '2.0 bin'
        self.dict['neigh_modify'] = 'every 10 delay 0 check no'
        self.dict['read_data'] = data_file 
        self.dict['pair_style'] = 'deepmd {}'.format(dp_graph_file)
        self.dict['pair_coeff'] = ''
        self.dict['velocity'] = 'all create {} 23456789'.format(temperature)
        self.dict['minimize'] = '1.0e-8 1.0e-6 100000 100000'
        self.dict['fix'] = '1 all nvt temp {} {} 0.5'.format(temperature, temperature)
        self.dict['timestep'] = str(timestep)
        self.dict['thermo_style'] = 'custom step pe ke etotal temp press vol'
        self.dict['thermo'] = '100'
        self.dict['dump'] = 'hzd_dump all custom 1000 geo.xyz type x y z'
        self.dict['run'] = str(run_num) 
        self.dict['write_data'] = '{}.dat'.format(sys_name)
        self.dict['write_restart'] = '{}.rest'.format(sys_name)

    def output(self):
        for k, v in self.dict.items():
            print(k, v) 

    def export(self, filename):
        f = open(filename, 'w')
        for k, v in self.dict.items():
            f.writelines('{} {}\n'.format(k, v))
        f.close()