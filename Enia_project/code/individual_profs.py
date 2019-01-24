import numpy as np
import struct
import time
import pickle
import os

import NBodySimulations
from colossus.utils import constants
from colossus.cosmology import cosmology
from colossus.halo import mass_so

###################################################################################################
# CONSTANTS
###################################################################################################

profile_dir = '/Users/eniaxhakaj/data/profiles/'

all_sims = ['L1000', 'L0500', 'L0250', 'L0125', 'L0063']

###################################################################################################

# This function works on a structured array of profiles which cannot have dimensionality 0. If only
# one profile is in the array, it must have dimensionality (1) which can be achieved by slicing an
# array as profiles[0:1], for example.

def density(prof):
	
	n_profiles = len(prof)
	n_bins = len(prof[0]['Rbin'])
	rho = np.zeros((n_profiles, n_bins), np.float)
	V = prof['Rbin']**3 * 4.0 * np.pi / 3.0
	dV = V[:, 1:] - V[:, :-1]
	dM = prof['mass'][:, 1:] - prof['mass'][:, :-1]
	rho[:, 0] = prof['mass'][:, 0] / V[:, 0]
	rho[:, 1:] = dM / dV
	
	return rho

###################################################################################################

def potential(prof):
	
	mask_nonzero = (prof['Rbin_av'] > 1E-10)
	integrand = prof['mass'] * constants.G
	integrand[mask_nonzero] /= prof['Rbin_av'][mask_nonzero]**2
	dr = np.array(prof['Rbin'])
	dr[:, 1:] -= prof['Rbin'][:, :-1]
	phi = np.cumsum(integrand * dr, axis = 1)
	phi = phi - phi[:, -1][:, None]
	
	return phi

###################################################################################################

def circularVelocity(prof):
	
	Vc = np.sqrt(prof['mass'] * constants.G / prof['Rbin'])
	
	return Vc

###################################################################################################

def root_s(prof, r1, r2, scale):

	n_prof = len(prof)
	r_scale = prof[scale]
	phi = prof['potential']
	idxh = np.arange(n_prof)
	
	R1 = r1 * r_scale
	R2 = r2 * r_scale
	r = prof['Rbin'][:, :]
	r1_idx = np.argmax(r > R1[:, None], axis = 1)
	r2_idx = np.argmax(r > R2[:, None], axis = 1)
	r1_0 = r[idxh, r1_idx - 1]
	r1_1 = r[idxh, r1_idx]
	r2_0 = r[idxh, r2_idx - 1]
	r2_1 = r[idxh, r2_idx]
	phi1_0 = phi[idxh, r1_idx - 1]
	phi1_1 = phi[idxh, r1_idx]
	phi2_0 = phi[idxh, r2_idx - 1]
	phi2_1 = phi[idxh, r2_idx]

	phi1 = phi1_0 + (R1 - r1_0) / (r1_1 - r1_0) * (phi1_1 - phi1_0)
	phi2 = phi2_0 + (R2 - r2_0) / (r2_1 - r2_0) * (phi2_1 - phi2_0)
	ret = np.sqrt(0.5 * (phi2 - phi1) / np.log(R2 / R1))
	
	return ret

###################################################################################################

# Compute vmax from the profiles simply as the maximum of vcirc within the virial radius. This can
# differ from the catalog values because of unbinding and because of the binning of the profiles.

def vmax(prof):
	
	Vc = circularVelocity(prof)
	mask_Rvir = (prof['Rbin'] > prof['Rvir_cat'][:, None])
	Vc[mask_Rvir] = 0.0
	ret = np.max(Vc, axis = 1)
	
	return ret

###################################################################################################

def load(sim_name, snap, Mvir_min = None, Mvir_max = None, Nvir_min = None, Nvir_max = None,
		mdefs = ['200c', '500c']):
	
	# Generate file name
	simProps = NBodySimulations.getSimulationProperties(sim_name)
	fn = 'Profiles_%s_Np200_%03d.bhp' % (simProps.folder, snap)	

	# Check for pickle
	pickle_file = 'Pickles/' + fn[:-4]
	if os.path.exists(pickle_file):
		
		t = time.clock()
		pFile = open(pickle_file, 'rb')
		dic = pickle.load(pFile)
		pFile.close()
		t = time.clock() - t
		print('Loaded %s from pickle, took %.2f seconds' % (fn, t))	
	
	else:
		
		# Open file, get basic info. Then reset the file pointer to the first struct.
		t = time.clock()
		f = open(profile_dir + fn, 'rb')
		n_profiles = struct.unpack('i', f.read(4))[0]
		n_bins = struct.unpack('i', f.read(4))[0]
		f.seek(4)
		struct_str = '3i13f%df' % (n_bins * 5)
		n_mdefs = len(mdefs)
				
		# Create a structure numpy array to hold all the data
		dtypes = []
		# Fields directly from the binaries
		dtypes.append(('id', np.long))
		dtypes.append(('host_id', np.long))
		dtypes.append(('Rvir_cat', np.float))
		dtypes.append(('rs', np.float))
		dtypes.append(('Mvir_cat', np.float))
		dtypes.append(('vmax_cat', np.float))
		dtypes.append(('a_last_mm', np.float))
		dtypes.append(('x', np.float, (3,)))
		dtypes.append(('v', np.float, (3,)))
		dtypes.append(('Rbin', np.float, (n_bins,)))
		dtypes.append(('Rbin_av', np.float, (n_bins,)))
		dtypes.append(('mass', np.float, (n_bins,)))
		dtypes.append(('vr', np.float, (n_bins,)))
		dtypes.append(('sigmav', np.float, (n_bins,)))
		# Pre-computed fields
		dtypes.append(('vc', np.float, (n_bins,)))
		dtypes.append(('vmax_prof', np.float))
		dtypes.append(('potential', np.float, (n_bins,)))
		# Mass definitions
		for j in range(n_mdefs):
			dtypes.append(('R%s' % mdefs[j], np.float))
			dtypes.append(('M%s' % mdefs[j], np.float))
		prof = np.zeros((n_profiles), dtype = dtypes)
	
		# Read the profiles one by one
		for i in range(n_profiles):
			p_data = struct.unpack(struct_str, f.read(struct.calcsize(struct_str)))
			if i == 0:
				z = p_data[3]
				box_size = p_data[4]
				
				# Prepare computation of mdefs; we need to know z for that
				thresholds = np.zeros((n_mdefs), np.float)
				for j in range(n_mdefs):
					thresholds[j] = mass_so.densityThreshold(z, mdefs[j])

			prof['id'][i] = p_data[1]
			prof['host_id'][i] = p_data[2]
			prof['Rvir_cat'][i] = p_data[5]
			prof['rs'][i] = p_data[6]
			prof['Mvir_cat'][i] = p_data[7]
			prof['vmax_cat'][i] = p_data[8]
			prof['a_last_mm'][i] = p_data[9]
			prof['x'][i] = np.array(p_data[10:13])
			prof['v'][i] = np.array(p_data[13:16])
			prof['Rbin'][i] = np.array(p_data[16 + 0 * n_bins:16 + 1 * n_bins])
			prof['Rbin_av'][i] = np.array(p_data[16 + 1 * n_bins:16 + 2 * n_bins])
			prof['mass'][i] = np.array(p_data[16 + 2 * n_bins:16 + 3 * n_bins])
			prof['vr'][i] = np.array(p_data[16 + 3 * n_bins:16 + 4 * n_bins])
			prof['sigmav'][i] = np.array(p_data[16 + 4 * n_bins:16 + 5 * n_bins])
	
		# Pre-compute various quantities
		prof['vc'] = circularVelocity(prof)
		prof['vmax_prof'] = vmax(prof)
		prof['potential'] = potential(prof)
	
		# Compute overdensity radii and masses in a fast, somewhat simplistic fashion. First, 
		# comptue the enclosed overdensity for all profiles and at all radii. Then, simply 
		# interpolate to find the radii. 
		delta = prof['mass'] / 4.0 / np.pi * 3.0 / prof['Rbin' ]**3
		max_delta = np.max(delta, axis = 1)
		min_delta = np.min(delta, axis = 1)
		
		for j in range(n_mdefs):
			RDelta = np.zeros((n_profiles), np.float)
			MDelta = np.zeros((n_profiles), np.float)
			
			# Filter out profiles that never reach teh density threshold
			mask_fail = (thresholds[j] < min_delta) | (thresholds[j] > max_delta)
			RDelta[mask_fail] = -1.0
			MDelta[mask_fail] = -1.0
			mask_valid = np.logical_not(mask_fail)
			delta_valid = delta[mask_valid, :]
			
			# Now find all bins that bracket the density threshold, and then use the lowest-radius
			# one of those
			mask_bins = (delta_valid[:, :-1] >= thresholds[j]) & (delta_valid[:, 1:] < thresholds[j])
			bin_idx = np.argmax(mask_bins, axis = 1)
			idx_halo = np.arange(0, len(bin_idx))
			r = prof['Rbin'][mask_valid]
			r0 = r[idx_halo, bin_idx]
			r1 = r[idx_halo, bin_idx + 1]
			d0 = delta_valid[idx_halo, bin_idx]
			d1 = delta_valid[idx_halo, bin_idx + 1]
			d_diff = d0 - d1
			mask_flat = (np.abs(d_diff) < 1E-10)
			RDelta[mask_valid][mask_flat] = r0[mask_flat]
			RDelta[mask_valid] = r0 + (r1 - r0) * (d0 - thresholds[j]) / d_diff
			MDelta[mask_valid] = mass_so.R_to_M(RDelta[mask_valid], z, mdefs[j])

			prof['R%s' % mdefs[j]] = RDelta
			prof['M%s' % mdefs[j]] = MDelta
			
			#ratio = RDelta[mask_valid] / prof['Rvir_cat'][mask_valid]
			#print(np.mean(ratio), np.median(ratio), np.min(ratio), np.max(ratio))

		f.close()
		t_load = time.clock() - t

		t = time.clock()
		dic = {}	
		dic['profiles'] = prof
		dic['box_size'] = box_size
		dic['z'] = z

		# Add simulation properties
		if simProps.box_size != box_size:
			raise Exception('Found different box sizes, %.4f and %.4f.' % (simProps.box_size, box_size))
		cosmo = cosmology.setCosmology(simProps.cosmo_name)
		mp = NBodySimulations.getParticleMass(simProps.box_size, simProps.Np, cosmo.Om0)
		dic['sim_name'] = sim_name
		dic['cosmo_name'] = simProps.cosmo_name
		dic['particle_mass'] = mp
		
		output_file = open(pickle_file, 'wb')
		pickle.dump(dic, output_file, pickle.HIGHEST_PROTOCOL)
		output_file.close()
		t_save = time.clock() - t

		print('Loaded %s from binary, took %.2f seconds, saving %.2f seconds.' % (fn, t_load, t_save))

	if Nvir_min is not None:
		Mvir_min = Nvir_min * dic['particle_mass']
	if Nvir_max is not None:
		Mvir_max = Nvir_max * dic['particle_mass']
	if Mvir_min is not None:
		Mvir = dic['profiles']['Mvir_cat'][:]
		dic['profiles'] = dic['profiles'][Mvir >= Mvir_min]
	if Mvir_max is not None:
		Mvir = dic['profiles']['Mvir_cat'][:]
		dic['profiles'] = dic['profiles'][Mvir <= Mvir_max]

	return dic

###################################################################################################
