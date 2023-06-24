from spacepy import pycdf
pycdf.lib.set_backward(False)
import datetime
import numpy as np
import interpret_tools
import os
from datetime import datetime
from math import pi

#global variables specifying CDF data variables:
lab_name = 'Name'
lab_dynamic = 'Dynamic'
lab_auth = 'Author' 
lab_date = 'CreateDate'
lab_axmu = 'axis_mu'
lab_axK = 'axis_K'
lab_axL = 'axis_L'
lab_axt = 'axis_t'
lab_axtd = 'axis_t_date'
lab_map = 'map_KL-aeq'
lab_f = 'f'
lab_j = 'j'
lab_axphi = 'axis_phi'
lab_en = 'energy'

class Preloaded:
	def __init__(self, fname_cdf):
		cdf = pycdf.CDF(fname_cdf)
		self.ax_t = cdf[lab_axt][:]
		self.ax_K = cdf[lab_axK][:]
		self.ax_L = cdf[lab_axL][:]
		self.map_alpha = cdf[lab_map][:]
		self.sol_en = cdf[lab_en][:]
		self.sol_f = cdf[lab_f][:]
		self.dynamic = cdf.attrs[lab_dynamic][0]

def get_phi_2015(L_):
	"""Get the average dipole field strength around Earth's equator and dipole moment. Use like so: B0,m = get_B0_m(2000.0)"""
	#from pyIGRF.loadCoeffs import get_coeffs
	#g, h = get_coeffs(2015)
	g10 = -29441.46 #g[1][0]
	g11 = -1501.77 #g[1][1]
	h11 = 4795.99 #h[1][1]

	B0_2 = g10**2 + g11**2 + h11**2
	B0_ = B0_2**0.5
	B0_ = B0_*1e-9 # = 2.986731323946967e-05 T

	RE = 6.3712e6 #m
	#mu0 = 1.25663706e-6 #

	phi = 2 * pi * B0_ * (RE ** 2) / L_ #T m2
	return phi


def convert_to_cdf(dir_solution_full, dynamic, overwrite = True):

	if dynamic:
		ext_txt = "_dyn.txt"
		ext_cdf = "_dyn.cdf"
	else:
		ext_txt = ".txt"
		ext_cdf = ".cdf"

	#check if a previous CDF file os.path.exists:
	fname_cdf = os.path.join(dir_solution_full, 'solution' + ext_cdf)
	print("converting solution to CDF file at {}".format(fname_cdf))
	if os.path.exists(fname_cdf):
		if overwrite:
			print("","warning: CDF file already exists, overwriting:")
			print("", fname_cdf)
			os.remove(fname_cdf)
		else:
			print("","warning: CDF file already exists, not updating")
			return fname_cdf
	else:
		print("","new file will be output to:")
		print("", fname_cdf)


	fname_axis_K = os.path.join(dir_solution_full, "axis_K"+ext_txt)
	Kaxisdata = interpret_tools.get_axis_file(fname_axis_K)
	solution_files_eachK = {}
	for idx_K, Kvalue in enumerate(Kaxisdata):
	    fileset = {}
	    fileset['K'] = Kvalue
	    fileset['axis mu'] = os.path.join(dir_solution_full, "axis_mu"+ext_txt)
	    fileset['axis K'] = os.path.join(dir_solution_full, "axis_K" + ext_txt)
	    fileset['axis L'] = os.path.join(dir_solution_full, "axis_L" + ext_txt)
	    fileset['axis t'] = os.path.join(dir_solution_full, "axis_t"+ext_txt)
	    fileset['config'] = os.path.join(dir_solution_full, "config"+ext_txt)
	    fileset['map'] = os.path.join(dir_solution_full, "map_iK-aeq"+ext_txt)
	    kslice_prefix = "iK-" + str(idx_K+1).zfill(4) + "_2D_"
	    fileset['f'] = os.path.join(dir_solution_full, kslice_prefix + "f"+ext_txt)
	    fileset['en'] = os.path.join(dir_solution_full, kslice_prefix + "en"+ext_txt)
	    fileset['data'] = os.path.join(dir_solution_full, kslice_prefix + "f_data"+ext_txt)
	    fileset['sdev'] = os.path.join(dir_solution_full, kslice_prefix + "f_sdev"+ext_txt)
	    fileset['prog'] = os.path.join(dir_solution_full, "progress"+ext_txt)
	    solution_files_eachK[idx_K] = fileset


	#read in axis files:
	axis_mu = interpret_tools.get_axis_file(solution_files_eachK[0]['axis mu'])
	axis_K = interpret_tools.get_axis_file(solution_files_eachK[0]['axis K'])
	axis_L = interpret_tools.get_axis_file(solution_files_eachK[0]['axis L'])
	axis_t = interpret_tools.get_axis_file(solution_files_eachK[0]['axis t'])
	if not dynamic: axis_t = np.array([axis_t[0]])
	map_aeq = interpret_tools.get_sol_file(solution_files_eachK[0]['map'])
	# convert timestamps to datetime objects (CDF_EPOCH):
	axis_t_date = []
	for ts in axis_t:
		axis_t_date.append(datetime.fromtimestamp(ts))
	axis_t_date = np.array(axis_t_date)
	axis_phi = np.array([get_phi_2015(L_) for L_ in axis_L])


	#open CDF file:
	cdf = pycdf.CDF(fname_cdf, '')
	#cdf.attrs[lab_name] = dir_local_solution
	cdf.attrs[lab_dynamic] = dynamic
	cdf.attrs[lab_auth] = 'A. R. Lozinski [alezin33@bas.ac.uk]'
	cdf.attrs[lab_date] = datetime.now()
	# copy axes:
	cdf[lab_axmu] = axis_mu
	cdf[lab_axmu].attrs['units'] = 'log10(mu/ (1 MeV/G) )'
	cdf[lab_axK] = axis_K
	cdf[lab_axK].attrs['units'] = 'G^0.5 RE'
	cdf[lab_axL] = axis_L
	cdf[lab_axt] = axis_t
	cdf[lab_axt].attrs['units'] = 'seconds since Jan 01 1970 (UTC)'
	cdf[lab_axtd] = axis_t_date
	cdf[lab_axtd].attrs['units'] = 'datetime CDF object'
	#cdf[lab_map] = map_aeq
	#cdf[lab_map].attrs['units'] = 'degrees'
	cdf[lab_axphi] = axis_phi
	cdf[lab_axphi].attrs['units'] = 'T m2'


	#detect whether or not map accounts for secular variation:
	secular_map = np.shape(map_aeq)[0]//np.size(axis_L) == np.size(axis_t)
	if secular_map:
		print("","detected that aeq map has secular variation")
		print()
		map_nt = np.size(axis_t)
	else:
		map_nt = 1

	#detect whether or not energy accounts for secular variation:
	secular_energy = np.shape(
        interpret_tools.get_sol_file(solution_files_eachK[0]['en']))[0] // np.size(axis_mu) == np.size(axis_t)
	if secular_energy:
		print("","detected that energy grid has secular variation")
		print()
		en_nt = np.size(axis_t)
	else:
		en_nt = 1


	data_map = np.zeros((map_nt, np.size(axis_L), np.size(axis_K)))
	for idx_t in range(map_nt):
		data_map[idx_t, :, :] = map_aeq[idx_t * np.size(axis_L) : (idx_t + 1) * (np.size(axis_L))]
	cdf[lab_map] = data_map
	cdf[lab_map].attrs['units'] = 'degrees'


	#create data grids:
	data_f = np.zeros((np.size(axis_t), np.size(axis_mu), np.size(axis_K), np.size(axis_L)))
	data_en = np.zeros((en_nt, np.size(axis_mu), np.size(axis_K), np.size(axis_L)))
	print("","data grid dimensions (t, mu, K, L): ", np.shape(data_f))

	print("","incorporating f, E grids...")
	for idx_K, Kvalue in enumerate(Kaxisdata):
		engrid_fixedK = interpret_tools.get_sol_file(solution_files_eachK[idx_K]['en'])
		soldata_fixedK = interpret_tools.get_sol_file(solution_files_eachK[idx_K]['f'])
		#print(solution_files_eachK[idx_K]['f'])#solution_files_eachK[idx_K]['en'])

		for idx_t in range(en_nt):
			data_en[idx_t, :, idx_K, :] = engrid_fixedK[idx_t * np.size(axis_mu) : (idx_t + 1) * (np.size(axis_mu))]
		for idx_t in range(len(axis_t)):
			data_f[idx_t, :, idx_K, :] = soldata_fixedK[idx_t * np.size(axis_mu) : (idx_t + 1) * (np.size(axis_mu))]

	cdf[lab_f] = data_f
	cdf[lab_f].attrs['desc'] = 'distribution function: relativistic phase space density multiplied by proton rest mass cubed'
	cdf[lab_f].attrs['units'] = 'km-6 s3'

	cdf[lab_en] = data_en
	cdf[lab_en].attrs['desc'] = 'energy at each data coordinate, assumed constant in time (ignoring secular variation)'
	cdf[lab_en].attrs['units'] = 'MeV'


	#translate pdf to flux too:
	data_j = np.zeros((np.size(axis_t), np.size(axis_mu), np.size(axis_K), np.size(axis_L)))
	for idx_t in range(len(axis_t)):
		idx_t_en = min(en_nt-1,idx_t)
		data_j[idx_t, data_en[idx_t_en]>0] = interpret_tools.f2j(data_en[idx_t_en, data_en[idx_t_en] > 0], data_f[idx_t, data_en[idx_t_en] > 0])
	cdf[lab_j] = data_j
	cdf[lab_j].attrs['desc'] = 'unidirectional differential proton flux'
	cdf[lab_j].attrs['units'] = 'cm-2 s-1 str-1 MeV-1'

	cdf.close()
	print("","done")
	return fname_cdf



# interpret_cdf.lab_name
# interpret_cdf.lab_dynamic
# interpret_cdf.lab_auth
# interpret_cdf.lab_date
# interpret_cdf.lab_axmu
# interpret_cdf.lab_axK
# interpret_cdf.lab_axL
# interpret_cdf.lab_phi
# interpret_cdf.lab_axt
# interpret_cdf.lab_axtd
# interpret_cdf.lab_map
# interpret_cdf.lab_f
# interpret_cdf.lab_j
# interpret_cdf.lab_en
