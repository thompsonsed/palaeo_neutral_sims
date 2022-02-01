# This version uses log sampling of the deme size.
import math
import os
import random
import shutil
import sys
import warnings

# from subprocess import call
from collections import OrderedDict
from paleo_occ_sq import *
from pycoalescence import Simulation, CoalescenceTree

if __name__ == "__main__":
    print("Starting simulation")
    # Change this value to change the type of simulation
    # Get the tmpdir directory set by the HPC for saving files to the fast disk.
    tmpdir = os.environ["TMPDIR"]
    homedir = os.environ["HOME"]
    output_dir = os.path.join(tmpdir, "Data")
    data_dir = os.path.join(homedir, "Data")
    global_save_dir = os.path.join(homedir, "Results/PaleoMainOcc/Fragmented3")
    if not os.path.exists(global_save_dir):
        os.mkdir(global_save_dir)
    if not os.path.exists(output_dir):
        if not os.path.exists("$TMPDIR"):
            raise IOError("$TMPDIR does not exist!")
        os.mkdir(output_dir)
    pbs_index = int(sys.argv[1])
    params = get_sim_parameters(pbs_index, global_save_dir, sim_type="fragmented")
    seed = params["seed"]
    job_type = params["job_type"]
    interval = params["interval"]
    tetrapod_group = params["tetrapod_group"]
    sigma = params["sigma"]
    density_per_km = params["density_per_km"]
    proportion_cover = params["proportion_cover"]
    save_dir = params["output_directory"]
    print("Parameters are:")
    for k, v in params.items():
        print(f"{k}: {v}")
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    fragment_dir = os.path.join(homedir, "Data/paleo_data/configs")
    # We will be simulating 53 individuals per cell, so need to scale our sample size accordingly
    random.seed(seed)
    sample_size, random_times = get_sample_size(
        interval, tetrapod_group, density_per_km
    )
    simulation = Simulation(logging_level=10)
    simulation.set_simulation_parameters(
        seed=seed,
        job_type=job_type,
        output_directory=output_dir,
        min_speciation_rate=10 ** -9,
        sigma=sigma,
        tau=0.0,
        deme=density_per_km,
        sample_size=sample_size,
        max_time=60 * 60 * 21.5,
        dispersal_method="normal",
        uses_spatial_sampling=True,
    )
    sample_file, fine_file, fragment_file = get_sim_type(
        interval,
        tetrapod_group,
        data_dir,
        fragment_dir,
        sim_type="fragmented",
        proportion_cover=proportion_cover,
    )
    try:
        simulation.set_map_files(
            sample_file=sample_file, fine_file=fine_file, coarse_file="none"
        )
    except IOError as ioe:
        print(
            "Map files don't exist - probably simulation is not required: {}".format(
                ioe
            )
        )
        exit(1)
    for t in random_times:
        simulation.add_sample_time(time=t)
    # simulation.optimise_ram(ram_limit=18)
    end_location = os.path.join(
        save_dir, "data_{}_{}.db".format(str(job_type), str(seed))
    )
    if os.path.exists(end_location):
        try:
            tmp_tree = CoalescenceTree()
            tmp_tree.set_database(end_location)
            print("Complete simulation detected... exiting")
            exit(1)
        except IOError:
            pass
            # Now move any paused files, if they exist.
    paused_files_list = ["Dump_main_{}_{}.csv".format(job_type, seed)]
    is_pause_sim = True
    for pause_file in paused_files_list:
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
        ct = CoalescenceTree()
        ct.speciate_remaining(simulation)
        # simulation.resume_coalescence(
        #     pause_directory=output_dir,
        #     seed=seed,
        #     job_type=job_type,
        #     max_time=60 * 60 * 71,
        #     out_directory=output_dir,
        # )
    else:
        simulation.run()

        # Now apply the speciation rates for each fragment
    speciation_file = os.path.join(
        simulation.output_directory, "data_{}_{}.db".format(job_type, seed)
    )
    try:
        speciation_rates = [
            0.000_000_01,
            0.000_000_05,
            0.000_000_1,
            0.000_000_5,
            0.000_001,
            0.000_005,
            0.00001,
            0.00005,
            0.0001,
            0.0005,
            0.001,
        ]
        # BCI_sim.set_speciation_params(speciation_file,"T","F",speciation_rates)
        simtree = CoalescenceTree()
        simtree.set_database(simulation)
        simtree.set_speciation_parameters(
            speciation_rates, record_fragments=fragment_file, record_spatial=True
        )
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
        print("done!")
    except Exception as e:
        print(str(e))
        print("Attempting data save..."),
        try:
            shutil.copy2(speciation_file, destination_file)
        except Exception as ec:
            print("failed! Trying other location.")
            print(str(ec))
            try:
                shutil.copy2(
                    speciation_file,
                    os.path.join(
                        "/work/set114/Results/Dump/",
                        "data_{}_{}.db".format(job_type, seed),
                    ),
                )
            except Exception as ec:
                print("also failed...")
                print("Broken! Quiting")
        else:
            print("success!")
            # Now try and calculate the tree structure

    print("Simulation complete, exiting")

