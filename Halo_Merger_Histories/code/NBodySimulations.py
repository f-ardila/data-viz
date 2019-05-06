import socket

from colossus.utils import constants
from colossus.cosmology import cosmology

###################################################################################################

# Define a few cosmologies that are used in the N-body simulations but not listed as default 
# cosmologies in the Cosmology class. The 'planck1-nbody' label was used in some codes before it
# was changed to 'planck13-nbody'.

sim_cosmos = {}
sim_cosmos['s01']            = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.0469, 'sigma8': 0.900, 'ns': 0.9500, 'relspecies': False}
sim_cosmos['s02']            = {'flat': True, 'H0': 70.0, 'Om0': 0.40, 'Ob0': 0.0469, 'sigma8': 0.820, 'ns': 0.9500, 'relspecies': False}
sim_cosmos['planck13-nbody'] = {'flat': True, 'H0': 67.0, 'Om0': 0.32, 'Ob0': 0.0491, 'sigma8': 0.834, 'ns': 0.9624, 'relspecies': False}
sim_cosmos['planck1-nbody']  = {'flat': True, 'H0': 67.0, 'Om0': 0.32, 'Ob0': 0.0491, 'sigma8': 0.834, 'ns': 0.9624, 'relspecies': False}

for c in sim_cosmos:
	cosmology.addCosmology(c, sim_cosmos[c])

###################################################################################################

sims = {}
sims['L2000']       	    = {'cosmo': 'bolshoi', 	        'box_size': 2000.0, 'Np': 1024, 'epsilon': 65.0,  'folder': 'Box_L2000_N1024_CBol',  	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L1000']       	    = {'cosmo': 'bolshoi', 	        'box_size': 1000.0, 'Np': 1024, 'epsilon': 33.0,  'folder': 'Box_L1000_N1024_CBol',  	    'Nfiles': 64,   'Nsnapdigits': 3,    'code': 'Gadget'}
sims['L0500']               = {'cosmo': 'bolshoi', 	        'box_size': 500.0,  'Np': 1024, 'epsilon': 14.0,  'folder': 'Box_L0500_N1024_CBol',  	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0250']      	        = {'cosmo': 'bolshoi', 	        'box_size': 250.0,  'Np': 1024, 'epsilon':  5.8,  'folder': 'Box_L0250_N1024_CBol',  	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0125']      	        = {'cosmo': 'bolshoi', 	        'box_size': 125.0,  'Np': 1024, 'epsilon':  2.4,  'folder': 'Box_L0125_N1024_CBol',  	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0063']       	    = {'cosmo': 'bolshoi',          'box_size': 62.5,   'Np': 1024, 'epsilon':  1.0,  'folder': 'Box_L0063_N1024_CBol', 	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0031']        	    = {'cosmo': 'bolshoi',          'box_size': 31.25,  'Np': 1024, 'epsilon':  0.25, 'folder': 'Box_L0031_N1024_CBol', 	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0500-Planck'] 	    = {'cosmo': 'planck13-nbody',   'box_size': 500.0,  'Np': 1024, 'epsilon': 14.0,  'folder': 'Box_L0500_N1024_CPla',  	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0250-Planck'] 	    = {'cosmo': 'planck13-nbody',   'box_size': 250.0,  'Np': 1024, 'epsilon':  5.8,  'folder': 'Box_L0250_N1024_CPla',  	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0125-Planck'] 	    = {'cosmo': 'planck13-nbody',   'box_size': 125.0,  'Np': 1024, 'epsilon':  2.4,  'folder': 'Box_L0125_N1024_CPla', 	    'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0250-s01']    	    = {'cosmo': 's01',              'box_size': 250.0,  'Np': 1024, 'epsilon':  5.8,  'folder': 'Box_L0250_N1024_Cs01',         'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0250-s02']    	    = {'cosmo': 's02',              'box_size': 250.0,  'Np': 1024, 'epsilon':  5.8,  'folder': 'Box_L0250_N1024_Cs02',         'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0125-CutPk']  	    = {'cosmo': 'bolshoi-cutpk',    'box_size': 125.0,  'Np': 1024, 'epsilon':  2.4,  'folder': 'Box_L0125_N1024_CBol_CutPk',   'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0063-CutPk']  	    = {'cosmo': 'bolshoi-cutpk',    'box_size': 62.5,   'Np': 1024, 'epsilon':  1.0,  'folder': 'Box_L0063_N1024_CBol_CutPk',   'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0031-CutPk']  	    = {'cosmo': 'bolshoi-cutpk',    'box_size': 31.25,  'Np': 1024, 'epsilon':  0.25, 'folder': 'Box_L0031_N1024_CBol_CutPk',   'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0100-plm10']  	    = {'cosmo': 'powerlaw_-1.0',    'box_size': 100.0,  'Np': 1024, 'epsilon':  0.5,  'folder': 'Box_L0100_N1024_CPLm10',       'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0100-plm15']  	    = {'cosmo': 'powerlaw_-1.5',    'box_size': 100.0,  'Np': 1024, 'epsilon':  0.5,  'folder': 'Box_L0100_N1024_CPLm15',       'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0100-plm20']  	    = {'cosmo': 'powerlaw_-2.0',    'box_size': 100.0,  'Np': 1024, 'epsilon':  1.0,  'folder': 'Box_L0100_N1024_CPLm20',       'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0100-plm25']  	    = {'cosmo': 'powerlaw_-2.5',    'box_size': 100.0,  'Np': 1024, 'epsilon':  1.0,  'folder': 'Box_L0100_N1024_CPLm25',       'Nfiles': 512,  'Nsnapdigits': 3,    'code': 'uberLGadget'}
sims['L0063-N0512']         = {'cosmo': 'bolshoi', 	        'box_size': 62.5,   'Np':  512, 'epsilon':  2.4,  'folder': 'Box_L0063_N0512_CBol',  	    'Nfiles': 128,  'Nsnapdigits': None, 'code': 'uberLGadget'}
sims['L0500-N0256']  	    = {'cosmo': 'bolshoi', 	        'box_size': 500.0,  'Np':  256, 'epsilon': 65.0,  'folder': 'Box_L0500_N0256_CBol',  	    'Nfiles': 16,   'Nsnapdigits': 4,    'code': 'uberLGadget'}
sims['L0125-N0256']  	    = {'cosmo': 'bolshoi', 	        'box_size': 125.0,  'Np':  256, 'epsilon': 14.0,  'folder': 'Box_L0125_N0256_CBol',  	    'Nfiles': 16,   'Nsnapdigits': None, 'code': 'uberLGadget'}
sims['L0125-N0256-future']  = {'cosmo': 'bolshoi', 	        'box_size': 125.0,  'Np':  256, 'epsilon': 14.0,  'folder': 'Box_L0125_N0256_CBol_future',  'Nfiles': 16,   'Nsnapdigits': 4,    'code': 'uberLGadget'}
sims['L0063-N0256']         = {'cosmo': 'bolshoi', 	        'box_size': 62.5,   'Np':  256, 'epsilon':  5.8,  'folder': 'Box_L0063_N0256_CBol',  	    'Nfiles': 16,   'Nsnapdigits': 4,    'code': 'uberLGadget'}
sims['TestSim100']          = {'cosmo': 'bolshoi', 	        'box_size': 62.5,   'Np':  256, 'epsilon':  5.8,  'folder': 'Box_L0063_N0256_CBol',  	    'Nfiles': 16,   'Nsnapdigits': 4,    'code': 'uberLGadget'}
sims['Bolshoi']     	    = {'cosmo': 'bolshoi', 	        'box_size': 250.0,  'Np': 2048, 'epsilon':  1.0,  'folder': 'Bolshoi',               	    'Nfiles': None, 'Nsnapdigits': None, 'code': 'ART'}
sims['Multidark']   	    = {'cosmo': 'bolshoi', 	        'box_size': 1000.0, 'Np': 2048, 'epsilon':  7.0,  'folder': 'Multidark',             	    'Nfiles': None, 'Nsnapdigits': None, 'code': 'ART'}
sims['MD-SMDPL']   	        = {'cosmo': 'multidark-planck', 'box_size':  400.0, 'Np': 3840, 'epsilon':  1.5,  'folder': 'MD-SMDPL',             	    'Nfiles': None, 'Nsnapdigits': None, 'code': 'Gadget2'}
sims['MD-MDPL']   	        = {'cosmo': 'multidark-planck', 'box_size': 1000.0, 'Np': 3840, 'epsilon':  5.0,  'folder': 'MD-MDPL',             	        'Nfiles': None, 'Nsnapdigits': None, 'code': 'Gadget2'}
sims['MD-BigMDPL']   	    = {'cosmo': 'multidark-planck', 'box_size': 2500.0, 'Np': 3840, 'epsilon': 10.0,  'folder': 'MD-BigMDPL',             	    'Nfiles': None, 'Nsnapdigits': None, 'code': 'Gadget2'}
sims['MD-HMDPL']   	        = {'cosmo': 'multidark-planck', 'box_size': 4000.0, 'Np': 4096, 'epsilon': 25.0,  'folder': 'MD-HMDPL',             	    'Nfiles': None, 'Nsnapdigits': None, 'code': 'Gadget2'}

###################################################################################################

class SimulationProperties:
	
	def __init__(self):
		
		self.name = ""
		self.box_size = 0.0
		self.Np = 0
		self.epsilon = 0.0
		self.mp = 0.0
		self.cosmo = None
		self.folder = None
		self.Nfiles = None
		self.Nsnapdigits = None
		self.code = None

		return
	
###################################################################################################

# This function defines a set of known simulations, with certain short
# string identifiers.

def getSimulationProperties(sim_name):
	
	if sim_name in sims:
		sim = SimulationProperties()
		sim.name = sim_name
		sim.cosmo_name = sims[sim_name]['cosmo']
		sim.box_size = sims[sim_name]['box_size']
		sim.Np = sims[sim_name]['Np']
		sim.epsilon = sims[sim_name]['epsilon']
		sim.folder = sims[sim_name]['folder']
		sim.Nfiles = sims[sim_name]['Nfiles']
		sim.Nsnapdigits = sims[sim_name]['Nsnapdigits']
		sim.code = sims[sim_name]['code']
		
	else:
		print("Invalid sim name, %s. Valid names are:" % sim_name)
		printBoxes()
		exit()

	try:
		previous_cosmo = cosmology.getCurrent()
	except Exception:
		previous_cosmo = None

	cosmo = cosmology.setCosmology(sim.cosmo_name)
	sim.mp = getParticleMass(sim.box_size, sim.Np, cosmo.Om0)
	
	if previous_cosmo is not None:
		cosmology.setCurrent(previous_cosmo)
	
	return sim

###################################################################################################

def getSimDir():
	
	name = socket.gethostname()
	if 'midway' in name:
		dir = '/home/diemer/Simulations/'
	else:
		dir = '/Users/benedito/University/Data/'
	
	return dir

###################################################################################################

def printBoxes():
	
	for c in list(sims.keys()):
		print(c)
	
	return

###################################################################################################

# Returns the mass of a particle in an N-Body simulation with box size
# L Mpc/h and N^3 particles.

def getParticleMass(L, N, Omega_M):
	
	return constants.RHO_CRIT_0_MPC3 * L**3 / N**3 * Omega_M

###################################################################################################
