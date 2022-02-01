import os
import shutil
# from subprocess import call
import sys
import warnings

import math

print("Starting module load...")

# from random import choice
# import zipfile
# Import Sahil Moza's implementation of latin-hypercube sampling
try:
	import lhsmdu
except:
	raise ImportError("Could not import the lhsmdu latin-hypercube sampling module.")

from pycoalescence import Simulation, CoalescenceTree
print("Module load complete")

def setup():
	"""
	Perform the set up routine, creating the necessary directories.
	"""
	if not os.path.exists("$TMPDIR/Data"):
		try:
			os.makedirs("$TMPDIR/Data")
		except:
			print("Could not create Data folder in $TMPDIR. Check write access.")


def choose_sim_variables(job_num, dimensions=2):
	"""
	Choose the relevant simulation variables given the command-line argument
	:param job_num: the command-line argument for setting the task
	:return: a list containing the important dispersal parameters
	"""
	# job_type = int(job_num) % 5 # For latin hypercube sampling every job type is different
	job_type = job_num
	if dimensions == 2:
		param_list = latin_two_dimensional_sampling(0.1, 20, 100, 1000, 25)
		return [job_num, job_type, param_list[0, job_num], param_list[1, job_num]]
	elif dimensions == 3:
		param_list = latin_three_dimensional_sampling(0.1, 8.0, 3.0, 10.0, 10, 6400, 100)
		return [job_num, job_type, param_list[0, job_num], param_list[1, job_num], param_list[2, job_num]]
	raise ValueError("Dimensions must be 2 or 3 currently.")



def latin_two_dimensional_sampling(dim_1_min, dim_1_max, dim_2_min, dim_2_max, number):
	"""
	Samples efficiently and completely from the parameter space specified. Returns a list of pairs of parameters which
	simulations should run with.
	:param dim_1_min: the minimum value of the first dimension
	:param dim_1_max: the maximum value of the first dimension
	:param dim_2_min: the minimum value of the second dimension
	:param dim_2_max: the maximum value of the second dimension
	:param number: the number of samples to draw in each dimension
	:return: a list of parameter pairs generated using the lhsmdu module.
	:rtype list
	"""
	samples = lhsmdu.sample(2, number, randomSeed=1001)
	assert (dim_1_min < dim_1_max)
	assert (dim_2_min < dim_2_max)
	samples[0] = (dim_1_max - dim_1_min) * samples[0] + dim_1_min
	samples[1] = (dim_2_max - dim_2_min) * samples[1] + dim_2_min
	return samples

def latin_three_dimensional_sampling(dim_1_min, dim_1_max, dim_2_min, dim_2_max, dim_3_min, dim_3_max, number):
	"""
	Samples efficiently and completely from the parameter space specified. Returns a list of triplets of parameters which
	simulations should run with.
	:param dim_1_min: the minimum value of the first dimension
	:param dim_1_max: the maximum value of the first dimension
	:param dim_2_min: the minimum value of the second dimension
	:param dim_2_max: the maximum value of the second dimension
	:param dim_3_min: the minimum value of the third dimension
	:param dim_3_max: the minimum value of the third dimension
	:param number: the number of samples to draw in each dimension
	:return: a list containing parameter triplets using the lhsmdu module
	:rtype list
	"""
	samples = lhsmdu.sample(3, number, randomSeed=1001)
	assert(dim_1_min < dim_1_max)
	assert(dim_2_min < dim_2_max)
	assert(dim_3_min < dim_3_max)
	samples[0] = (dim_1_max - dim_1_min) * samples[0] + dim_1_min
	samples[1] = (dim_2_max - dim_2_min) * samples[1] + dim_2_min
	samples[2] = (dim_3_max - dim_3_min) * samples[2] + dim_3_min
	return samples


def copyall(src, dst):
	"""
	Copies the source file or directory to the new file or directory. Directories are not overwritten and are copied
	recursively.
	:param src: source file or directory
	:param dst: destination file or directory
	"""
	if not os.path.exists(os.path.dirname(dst)):
		os.makedirs(os.path.dirname(dst))
	if os.path.isfile(src):
		shutil.copy2(src, dst)
	elif os.path.isdir(src):
		if not os.path.exists(dst):
			os.mkdir(dst)
		for item in os.listdir(src):
			s = os.path.join(src, item)
			d = os.path.join(dst, item)
			if os.path.isdir(s):
				copyall(s, d)
			elif os.path.isfile(s):
				shutil.copy2(s, d)
	else:
		raise RuntimeError("Source is not a file or directory.")


def get_sim_type(sim_reference, data_directory, landscape_ref):
	"""
	For use with the paleo simulations, getting the sample file locations of the relevant files
	:param sim_reference: should be one of either "permian", "early_carboniferous" or "late_carboniferous"
	:param data_directory: location of the Data folder, containing the paleo_maps folder
	:param landscape_ref: landscpe reference number
	:return: list containing the locations of the map files
	"""
	if sim_reference in ["permian", "early_carboniferous", "late_carboniferous"]:
		sample_file = os.path.join(data_directory, "paleo_maps", "paleomask_{}.tif".format(sim_reference))
		cover = cover_list[landscape_ref]
		fine_file = os.path.join(data_directory, "paleo_maps", "{}_fragmented_{}_masked.tif".format(sim_reference,
																									cover))
		return [sample_file, fine_file]
	raise ValueError("{} not one of permian, early_carboniferous or late_carboniferous.".format(sim_reference))

cover_list = [0.1, 0.2, 0.5]
offsets_dict = {
	"permian": {'grid_file_name': 'set',
				 'grid_x_size': 9453,
				 'grid_y_size': 9453,
				 'sample_x_offset': 0,
				 'sample_y_offset': 0},
	"early_carboniferous": {'grid_file_name': 'set',
							 'grid_x_size': 4006,
							 'grid_y_size': 4006,
							 'sample_x_offset': 0,
							 'sample_y_offset': 0},
	"late_carboniferous": {'grid_file_name': 'set',
							'grid_x_size': 5809,
							'grid_y_size': 5809,
							'sample_x_offset': 0,
							'sample_y_offset': 0}
}
max_density = {"permian": 37, "early_carboniferous": 5, "late_carboniferous": 34}

if __name__ == "__main__":
	print("Starting simulation")
	# Change this value to change the type of simulation
	sim_list = ["permian", "early_carboniferous", "late_carboniferous"]
	# Get the tmpdir directory set by the HPC for saving files to the fast disk.
	tmpdir = os.environ['TMPDIR']
	output_dir = os.path.join(tmpdir, "Data")
	data_dir = "/work/set114/Panama/Data/"
	save_dir = "/work/set114/Panama/Results/Paleo_07/"
	if not os.path.exists(save_dir):
		os.mkdir(save_dir)
	if not os.path.exists(output_dir):
		if not os.path.exists("$TMPDIR"):
			raise IOError("$TMPDIR does not exist!")
		os.mkdir(output_dir)
	pbs_index = int(sys.argv[1])
	job_num = int(pbs_index % 25)  # we are running 25 simulations for each sim type
	landscape_ref = int(math.floor(pbs_index % 3))
	job_type = int(math.floor(pbs_index / 25))
	sim_type = sim_list[job_type]
	fragment_dir = "/work/set114/Panama/Configs/fragments_{}.csv".format(sim_type)
	_, _, sigma, density_per_km = choose_sim_variables(job_num, dimensions=2)
	# We will be simulating 53 individuals per cell, so need to scale our sample size accordingly
	sample_size = max_density[sim_type]/density_per_km
	simulation = Simulation(logging_level=10)
	simulation.set_simulation_params(seed=job_num, job_type=job_type, output_directory=output_dir,
									 min_speciation_rate=10**-8, sigma=sigma, tau=0.0, deme=density_per_km,
									 sample_size=sample_size, max_time=60*60*71.5, dispersal_method="normal",
									 uses_spatial_sampling=True)
	sample_file, fine_file = get_sim_type(sim_type, data_dir, landscape_ref)
	simulation.set_map_files(sample_file=sample_file, fine_file=fine_file, coarse_file="none")
	simulation.optimise_ram(ram_limit=18)
	simulation.finalise_setup()
	end_location = os.path.join(save_dir, "SQL_data", "data_{}_{}.db".format(str(job_type), str(job_num)))
	if os.path.exists(end_location):
		try:
			tmp_tree = CoalescenceTree()
			tmp_tree.set_database(end_location)
			print("Complete simulation detected... exiting")
			exit(1)
		except IOError:
			pass
	# Now move any paused files, if they exist.
	paused_files_list = ["Dump_active_", "Dump_data_", "Dump_main_", "Dump_map_"]
	is_pause_sim = True
	for pause_file in paused_files_list:
		if not os.path.exists(os.path.join(save_dir, "Pause", (pause_file + str(job_type) +"_"+str(job_num)+".csv"))):
			is_pause_sim = False
	if is_pause_sim:
		print("Paused sim detected... resuming")
		os.mkdir(os.path.join(output_dir, "Pause"))
		# open(os.path.join(output_dir, "Pause", ("Dump_grid_" + str(job_num) +"_"+str(job_type)+".csv")), 'a').close()
		for pause_file in paused_files_list:
			pause_file_src = os.path.join(save_dir, "Pause", (pause_file + str(job_type) +"_"+str(job_num)+".csv"))
			pause_file_dst = os.path.join(output_dir, "Pause", (pause_file + str(job_type) +"_"+str(job_num)+".csv"))
			shutil.copy(pause_file_src, pause_file_dst)
		print("Files copied successfully")
		simulation.resume_coalescence(pause_directory=output_dir, seed=job_num, job_type=job_type, max_time=60*60*71,
									  out_directory=output_dir)
	else:
		simulation.run_coalescence()

	# Now apply the speciation rates for each fragment
	speciation_file = simulation.output_directory + "/SQL_data/data_" + str(job_type) + "_" + str(job_num) + ".db"
	try:
		speciation_rates = [0.00000001, 0.00000005, 0.0000001, 0.0000005, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001,
							0.0005, 0.001]
		# BCI_sim.set_speciation_params(speciation_file,"T","F",speciation_rates)
		simtree = CoalescenceTree()
		simtree.set_database(simulation)
		simtree.set_speciation_params("T", fragment_dir, speciation_rates)
		simtree.apply_speciation()
	except IOError:
		pass
	except Exception as e:
		warnings.warn(str(e))
		warnings.warn("Could not complete speciation.")
	print("------------------")
	destination_file = save_dir + "data_" + str(job_type) + "_" + str(job_num) + ".db"
	try:
		print("Copying files..."),
		copyall(output_dir, save_dir)
		print('done!')
	except Exception as e:
		print(str(e))
		print('Attempting data save...'),
		try:
			shutil.copy2(speciation_file, destination_file)
		except Exception as ec:
			print('failed! Trying other location.')
			print(str(ec))
			try:
				shutil.copy2(speciation_file,
							 "/work/set114/Panama/Results/Dump/data_" + str(job_type) + "_" + str(job_num) + ".db")
			except Exception as ec:
				print("also failed...")
				print("Broken! Quiting")
		else:
			print('success!')
	# Now try and calculate the tree structure

	print("Simulation complete, exiting")