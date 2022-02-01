 # This version uses log sampling of the deme size.
import os
import shutil
import sys
import warnings
# from subprocess import call
from collections import OrderedDict

print("Starting module load...")

# from random import choice
# import zipfile
# Import Sahil Moza's implementation of latin-hypercube sampling
try:
	import lhsmdu
except:
	raise ImportError("Could not import the lhsmdu latin-hypercube sampling module.")

from pycoalescence import CoalescenceTree, Simulation

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
		param_list = latin_two_dimensional_sampling(0.1, 20, 1.4, 3, 25)
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
	assert (dim_1_min < dim_1_max)
	assert (dim_2_min < dim_2_max)
	assert (dim_3_min < dim_3_max)
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


def get_sim_parameters(job_num, save_dir):
	"""
	Gets the sim parameters, returned as a dictionary.

	:param job_num: the value passed in sys.argv
	:param save_dir: the directory to save the sim output to
	:return: dictionary containing the simulation parameters
	"""
	sim_list = []
	job_type = 1
	for proportion_cover in [20, 40, 80]:
		for interval in intervals:
			for tetra_group in tetrapod_groups:
				job_type += 1
				for seed in range(1, 26, 1):
					sim_list.append([interval, tetra_group, seed, job_type, proportion_cover])
	this_sim = sim_list[job_num]
	_, _, sigma, density_per_km = choose_sim_variables(this_sim[2], dimensions=2)
	density_per_km = 10 ** density_per_km
	return {"interval": this_sim[0], "tetrapod_group": this_sim[1], "sigma": sigma,
			"seed": this_sim[2], "job_type": this_sim[3], "density_per_km": density_per_km,
			"proportion_cover": this_sim[4], "output_directory": os.path.join(save_dir, str(this_sim[4]))}


def get_sim_type(interval, tetrapod_group, data_directory, fragment_dir, proportion_cover):
	"""
	For use with the paleo simulations, getting the sample file locations of the relevant files
	:param interval: should be one the intervals
	:param tetrapod_group: should be one of the tetrapod groups
	:param data_directory: location of the Data folder, containing the paleo_maps folder
	:param proportion_cover: the proportion of habitat cover

	:return: list containing the locations of the map files
	"""
	if interval in intervals and tetrapod_group in tetrapod_groups:
		sample_file = os.path.join(data_directory, "paleo_maps", "paleomask_{}_{}.tif".format(interval, tetrapod_group))
		fine_file = os.path.join(data_directory, "paleo_maps", "fragmented",
								 "{}_{}_fragmented.tif".format(interval, round(proportion_cover / 100, 1)))
		fragment_file = os.path.join(fragment_dir, "fragments_{}_{}.csv".format(interval, tetrapod_group))
		return [sample_file, fine_file, fragment_file]
	raise ValueError(
		"{} not one of {} or {} not one of {}.".format(interval, intervals, tetrapod_group, tetrapod_groups))


max_density = OrderedDict({('artinskian', 'amniote'): 15,
						   ('artinskian', 'amphibian'): 15,
						   ('asselian', 'amniote'): 6,
						   ('asselian', 'amphibian'): 6,
						   ('bashkirian', 'amniote'): 1,
						   ('bashkirian', 'amphibian'): 23,
						   ('gzhelian', 'amniote'): 7,
						   ('gzhelian', 'amphibian'): 7,
						   ('kasimovian', 'amniote'): 7,
						   ('kasimovian', 'amphibian'): 4,
						   ('kungurian', 'amniote'): 15,
						   ('kungurian', 'amphibian'): 17,
						   ('moscovian', 'amniote'): 4,
						   ('moscovian', 'amphibian'): 21,
						   ('sakmarian', 'amniote'): 6,
						   ('sakmarian', 'amphibian'): 15})

intervals = ['artinskian',
			 'asselian',
			 'bashkirian',
			 'gzhelian',
			 'kasimovian',
			 'kungurian',
			 'moscovian',
			 'sakmarian']

tetrapod_groups = ['amniote', 'amphibian']

if __name__ == "__main__":
	print("Starting simulation")
	# Change this value to change the type of simulation
	# Get the tmpdir directory set by the HPC for saving files to the fast disk.
	tmpdir = os.environ['TMPDIR']
	homedir= os.environ['HOME']
	output_dir = os.path.join(tmpdir, "Data")
	data_dir = os.path.join(homedir, "Data")
	global_save_dir = os.path.join(homedir, "Results/PaleoMain/Sim8")
	if not os.path.exists(global_save_dir):
		os.mkdir(global_save_dir)
	if not os.path.exists(output_dir):
		if not os.path.exists("$TMPDIR"):
			raise IOError("$TMPDIR does not exist!")
		os.mkdir(output_dir)
	pbs_index = int(sys.argv[1])
	params = get_sim_parameters(pbs_index, global_save_dir)
	seed = params["seed"]
	job_type = params["job_type"]
	interval = params["interval"]
	tetrapod_group = params["tetrapod_group"]
	if interval != "kasimovian" and tetrapod_group != "amniote": # TODO remove for future sims
		exit(1)
	sigma = params["sigma"]
	density_per_km = params["density_per_km"]
	proportion_cover = params["proportion_cover"]
	save_dir = params["output_directory"]
	if not os.path.exists(save_dir):
		os.mkdir(save_dir)
	fragment_dir = os.path.join(homedir, "Data/paleo_data/configs")
	# We will be simulating 53 individuals per cell, so need to scale our sample size accordingly
	sample_size = min(10*max_density[(interval, tetrapod_group)] / density_per_km, 1.0)
	simulation = Simulation(logging_level=10)
	simulation.set_simulation_parameters(seed=seed, job_type=job_type, output_directory=output_dir,
	                                     min_speciation_rate=10 ** -8, sigma=sigma, tau=0.0, deme=density_per_km,
	                                     sample_size=sample_size, max_time=60 * 60 * 22.5, dispersal_method="normal",
	                                     uses_spatial_sampling=True)
	sample_file, fine_file, fragment_file = get_sim_type(interval, tetrapod_group, data_dir, fragment_dir,
														 proportion_cover)
	simulation.set_map_files(sample_file=sample_file, fine_file=fine_file, coarse_file="none")
	# simulation.optimise_ram(ram_limit=18)
	end_location = os.path.join(save_dir, "data_{}_{}.db".format(str(job_type), str(seed)))
	if os.path.exists(end_location):
		try:
			tmp_tree = CoalescenceTree()
			tmp_tree.set_database(end_location)
			print("Complete simulation detected... exiting")
			exit(1)
		except IOError:
			pass
	# Now move any paused files, if they exist.
	paused_files_list = ["Dump_main_{}_{}.csv"]
	is_pause_sim = True
	for pause_file in paused_files_list:
		pause_file = pause_file.format(job_type, seed)
		if not os.path.exists(os.path.join(save_dir, "Pause", pause_file)):
			is_pause_sim = False
	if is_pause_sim:
		print("Paused sim detected... resuming")
		os.mkdir(os.path.join(output_dir, "Pause"))
		for pause_file in paused_files_list:
			pause_file_src = os.path.join(save_dir, "Pause", pause_file)
			pause_file_dst = os.path.join(output_dir, "Pause", pause_file)
			shutil.copy(pause_file_src, pause_file_dst)
		print("Files copied successfully")
		simulation.resume_coalescence(pause_directory=output_dir, seed=seed, job_type=job_type, max_time=60 * 60 * 71,
									  out_directory=output_dir)
	else:
		simulation.run()

	# Now apply the speciation rates for each fragment
	speciation_file = os.path.join(simulation.output_directory, "data_{}_{}.db".format(job_type, seed))
	try:
		speciation_rates = [0.00000001, 0.00000005, 0.0000001, 0.0000005, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001,
							0.0005, 0.001]
		# BCI_sim.set_speciation_params(speciation_file,"T","F",speciation_rates)
		simtree = CoalescenceTree()
		simtree.set_database(simulation)
		simtree.set_speciation_parameters(speciation_rates, record_fragments=fragment_file, record_spatial=True)
		simtree.apply()
	except IOError:
		pass
	except Exception as e:
		warnings.warn(str(e))
		warnings.warn("Could not complete speciation.")
	print("------------------")
	destination_file = os.path.join(save_dir, "data_{}_{}.db".format(job_type, seed))
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
							 os.path.join("/work/set114/Results/Dump/",
										  "data_{}_{}.db".format(job_type, seed)))
			except Exception as ec:
				print("also failed...")
				print("Broken! Quiting")
		else:
			print('success!')
	# Now try and calculate the tree structure

	print("Simulation complete, exiting")
