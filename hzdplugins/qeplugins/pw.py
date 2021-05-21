class pwInputGenerator():

    def __init__(self, control={}, system={}, electrons={}, ions={}, cell={}):
        self.control = control 
        self.system = system 
        self.electrons = electrons
        self.ions = ions 
        self.cell = cell
        self.dict = {}
        self.dict['CONTROL'] = self.control
        self.dict['SYSTEM'] = self.system
        self.dict['ELECTRONS'] = self.electrons
        self.dict['IONS'] = self.ions
        self.dict['CELL'] = self.cell
    
    def setControl(self, calculation = 'vc-relax', restart_mode = 'from_scratch', wf_collect = True, nstep = 50000, tstress = True, tprnfor = True, disk_io = 'low'):
        self.control['calculation'] = calculation 
        self.control['restart_mode'] = restart_mode
        self.control['wf_collect'] = wf_collect
        self.control['nstep'] = nstep
        self.control['tstress'] = tstress
        self.control['tprnfor'] = tprnfor
        self.control['disk_io'] = disk_io
        self.dict['CONTROL'] = self.control
    
    def setSystem(self, ibrav = 0, ecutwfc = 50.0, ecutrho = 500.0, occupations = 'smearing', degauss = 0.002, smearing = 'gaussian', input_dft = 'pbesol', lda_plus_u = False, nspin = None, hubbard_u = {}, starting_magnetization = {}, starting_ns_eigenvalue = [], tot_magnetization = None, nbnd = None):
        self.system['ibrav'] = ibrav 
        self.system['ecutwfc'] = ecutwfc
        self.system['ecutrho'] = ecutrho
        self.system['occupations'] = occupations
        self.system['smearing'] = smearing
        self.system['degauss'] = degauss
        self.system['input_dft'] = input_dft

        # if we want to use DFT+U method
        if lda_plus_u:
            self.system['lda_plus_u'] = lda_plus_u
            self.system['hubbard_u'] = hubbard_u

        # if we want to include spin
        if nspin == 2:
            self.system['nspin'] = 2
            self.system['starting_magnetization'] = starting_magnetization
            if len(starting_ns_eigenvalue) > 0:
                self.system['starting_ns_eigenvalue'] = starting_ns_eigenvalue
            if not tot_magnetization is None:
                self.system['tot_magnetization'] = tot_magnetization
            if not nbnd is None:
                self.system['nbnd'] = nbnd
        
        self.dict['SYSTEM'] = self.system
    
    def setElectrons(self, conv_thr = 1e-06, diagonalization = 'david', electron_maxstep = 1000, mixing_fixed_ns = None, mixing_beta = 2e-01, mixing_mode = 'plain', mixing_ndim = 10):
       self.electrons['conv_thr'] = conv_thr
       self.electrons['diagonalization'] = diagonalization
       self.electrons['electron_maxstep'] = electron_maxstep
       if not mixing_fixed_ns is None:
          self.electrons['mixing_fixed_ns'] = mixing_fixed_ns
       self.electrons['mixing_beta'] = mixing_beta
       self.electrons['mixing_mode'] = mixing_mode
       self.electrons['mixing_ndim'] = mixing_ndim

       self.dict['ELECTRONS'] = self.electrons
    
    def setIons(self):
        self.ions = {}
        self.dict['IONS'] = self.ions
    
    def setCell(self):
        self.cell = {}
        self.dict['CELL'] = self.cell
    
    def setStructure(self, struct):
        # This function will help build: ATOMIC_SPECIES, ATOMIC_POSITIONS, CELL_PARAMETERS 
        self.dict['ATOMIC_SPECIES'] = 
        self.dict['ATOMIC_POSITIONS'] = 
        self.dict['CELL_PARAMETERS'] = 

    def setKpoints(self, kpts, disp):
        self.dict['KPOINTS'] = [kpts, disp]

    def output(self):
        import pprint
        pp = pprint.PrettyPrinter(indent = 4)
        pp.pprint(self.dict)
    
    def export(self, filename):
        import json
        f = open(filename, 'w')
        json.dump(self.dict, f, indent=4)