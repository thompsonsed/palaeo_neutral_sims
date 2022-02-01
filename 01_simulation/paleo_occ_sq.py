"""
Contains all the functions for setting up spatially explicit simulations of tetrapod diversity.
"""

# from subprocess import call
import os
import shutil
import math
from collections import OrderedDict
import random

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


max_density = OrderedDict([(('artinskian', 'amniote'), 225),
                           (('artinskian', 'amphibian'), 225),
                           (('asselian', 'amniote'), 36),
                           (('asselian', 'amphibian'), 36),
                           (('bashkirian', 'amniote'), 1),
                           (('bashkirian', 'amphibian'), 529),
                           (('gzhelian', 'amniote'), 49),
                           (('gzhelian', 'amphibian'), 49),
                           (('kasimovian', 'amniote'), 49),
                           (('kasimovian', 'amphibian'), 16),
                           (('kungurian', 'amniote'), 225),
                           (('kungurian', 'amphibian'), 289),
                           (('moscovian', 'amniote'), 16),
                           (('moscovian', 'amphibian'), 441),
                           (('sakmarian', 'amniote'), 36),
                           (('sakmarian', 'amphibian'), 225)])

intervals = ["artinskian", "asselian", "bashkirian", "gzhelian", "kasimovian", "kungurian", "moscovian", "sakmarian"]

tetrapod_groups = ["amniote", "amphibian"]

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
    assert dim_1_min < dim_1_max
    assert dim_2_min < dim_2_max
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
    assert dim_1_min < dim_1_max
    assert dim_2_min < dim_2_max
    assert dim_3_min < dim_3_max
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


def get_sim_parameters(job_num, save_dir, sim_type="pristine"):
    """
	Gets the sim parameters, returned as a dictionary.

	:param job_num: the value passed in sys.argv
	:param save_dir: the directory to save the sim output to
	:return: dictionary containing the simulation parameters
	"""
    sim_list = []
    job_type = 1
    proportion_covers = [20, 40, 80] if sim_type == "fragmented" else [100]
    for proportion_cover in proportion_covers:
        for interval in intervals:
            for tetra_group in tetrapod_groups:
                job_type += 1
                for seed in range(1, 26, 1):
                    sim_list.append([interval, tetra_group, seed, job_type, proportion_cover])
    this_sim = sim_list[job_num]
    _, _, sigma, density_per_km = choose_sim_variables(this_sim[2], dimensions=2)
    density_per_km = 10 ** density_per_km
    return {
        "interval": this_sim[0],
        "tetrapod_group": this_sim[1],
        "sigma": sigma,
        "seed": this_sim[2],
        "job_type": this_sim[3],
        "density_per_km": density_per_km,
        "proportion_cover" : this_sim[4],
        "output_directory": os.path.join(save_dir, str(this_sim[3])),
    }


def get_sim_type(interval, tetrapod_group, data_directory, fragment_dir, sim_type, proportion_cover = 1.0):
    """
	For use with the paleo simulations, getting the sample file locations of the relevant files
	:param interval: should be one the intervals
	:param tetrapod_group: should be one of the tetrapod groups
	:param data_directory: location of the Data folder, containing the paleo_maps folder

	:return: list containing the locations of the map files
	"""
    if interval in intervals and tetrapod_group in tetrapod_groups:
        sample_file = os.path.join(data_directory, "paleo_maps",
                                   "paleomask_occ_sq_{}_{}.tif".format(interval, tetrapod_group))
        if sim_type == "clustered":
            fine_file = os.path.join(
                data_directory, "paleo_maps", "clustered", "{}_{}.tif".format(interval, tetrapod_group)
            )
        elif sim_type == "fragmented":
            fine_file = os.path.join(data_directory, "paleo_maps", "fragmented",
                                     "{}_{}_fragmented.tif".format(interval, round(proportion_cover / 100, 1)))
        else:
            fine_file = os.path.join(data_directory, "paleo_maps", "{}.tif".format(interval, tetrapod_group))
        fragment_file = os.path.join(fragment_dir, "fragments_occ_{}_{}.csv".format(interval, tetrapod_group))
        return [sample_file, fine_file, fragment_file]
    raise ValueError(
        "{} not one of {} or {} not one of {}.".format(interval, intervals, tetrapod_group, tetrapod_groups)
    )


generational_time = 10
timespan_years = 1000


def get_sample_size(interval, tetrapod_group, density_per_km):
    # distribute the samples evenly in time.

    # then chose a random number of times between the min and maximum.
    min_number_times = max(math.ceil(max_density[(interval, tetrapod_group)]/density_per_km), 1)
    number_times = random.randint(min_number_times, max(10, min_number_times))
    max_individuals_simulated = density_per_km * number_times
    # Reduce our sample size at each point in time so that we don't simulate too many individuals
    sample_size = min((1 + 1.05*max_density[(interval, tetrapod_group)])/max_individuals_simulated, 1.0)
    random_times = [random.randint(0, int(timespan_years/generational_time)) for _ in range(0, number_times)]
    return sample_size, random_times