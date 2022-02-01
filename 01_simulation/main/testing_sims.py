import os
from collections import OrderedDict
from pycoalescence import Simulation

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

tetrapod_groups = {'amniote', 'amphibian'}
output_dir = "/Volumes/Seagate 3TB/Temp/"
src_dir = "/Volumes/Seagate 3TB/Paleo/Data/paleo_maps/main"


params = {'interval': 'sakmarian',
 'tetrapod_group': 'amniote',
 'sigma': 13.086113705647051,
 'seed': 17,
 'job_type': 16,
 'density_per_km': 40.40260281502476}

def get_sim_type(interval, tetrapod_group, data_directory, fragment_dir):
	"""
	For use with the paleo simulations, getting the sample file locations of the relevant files
	:param interval: should be one the intervals
	:param tetrapod_group: should be one of the tetrapod groups
	:param data_directory: location of the Data folder, containing the paleo_maps folder
	:return: list containing the locations of the map files
	"""
	if interval in intervals and tetrapod_group in tetrapod_groups:
		sample_file = os.path.join(data_directory, "paleomask_{}_{}.tif".format(interval, tetrapod_group))
		fine_file = os.path.join(data_directory, "{}.tif".format(interval, tetrapod_group))
		fragment_file = os.path.join(fragment_dir, "fragments_{}_{}.csv".format(interval, tetrapod_group))
		return [sample_file, fine_file, fragment_file]
	raise ValueError("{} not one of {} or {} not one of {}.".format(interval, intervals, tetrapod_group, tetrapod_groups))


seed = params["seed"]
job_type = params["job_type"]
interval = params["interval"]
# tetrapod_group = params["tetrapod_group"]
tetrapod_group = "amphibian"
sigma = params["sigma"]
density_per_km = params["density_per_km"]
sample_size = max_density[(interval, tetrapod_group)]/density_per_km
fragment_dir = "none"
sample_file, fine_file, fragment_file = get_sim_type(interval, tetrapod_group, src_dir, fragment_dir)
fine_map_x = 32072
fine_map_y = 15486
fine_map_x_offset = 10580
fine_map_y_offset = 4380

print("deme: {}, sampling {}".format(density_per_km, sample_size))
simulation = Simulation(logging_level=10)
simulation.set_simulation_params(seed=seed, job_type=job_type, output_directory=output_dir,
								 min_speciation_rate=10**-5, sigma=sigma, tau=0.0,
								 deme=density_per_km,
								 sample_size=sample_size, max_time=60*60*71.5, dispersal_method="normal",
								 uses_spatial_sampling=True)
simulation.set_map_files(sample_file=sample_file, fine_file=fine_file, coarse_file="none")
# THis should produce 130 individuals!
# simulation.set_map_parameters(sample_file=sample_file,sample_x=10431,
# 							  sample_y=8078,
# 							  fine_file="null", fine_x=fine_map_x, fine_y=fine_map_y,
# 							  fine_x_offset=fine_map_x_offset, fine_y_offset=fine_map_y_offset,
# 							  coarse_file="none",
# 							  coarse_x=fine_map_x, coarse_y=fine_map_y, coarse_x_offset=0, coarse_y_offset=0,
# 							  coarse_scale=1.0, historical_fine_map="none", historical_coarse_map="none")
try:
	simulation.optimise_ram(ram_limit=18)
except Exception as e:
	raise e
simulation.run()