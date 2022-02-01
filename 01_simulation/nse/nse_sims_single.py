import os
import sys
from operator import itemgetter

import numpy as np
import pandas as pd
import random
from pycoalescence import CoalescenceTree

intervals = {'artinskian',
 'asselian',
 'bashkirian',
 'gzhelian',
 'kasimovian',
 'kungurian',
 'moscovian',
 'sakmarian'}

tetrapod_groups = {'amniote', 'amphibian'}

immigration_rates = 10.0**-np.linspace(0.1, 4, 20)
speciation_rates = [10**-x for x in range(1, 9)]

# The location to store simulations in
output_dir = "/rds/general/user/set114/home/Results/PaleoMain/nse_01_single" # TODO modify
if not os.path.exists(output_dir):
	os.makedirs(output_dir)
# The location of the csv containing location data for the fossil sites
csv_directory = "/rds/general/user/set114/home/Data/paleo_data" # TODO modify
# The location to output the results csv to
output_csv_directory = "/rds/general/user/set114/home/Results/PaleoMain/nse_01_single" # TODO modify
if not os.path.exists(output_csv_directory):
	os.makedirs(output_csv_directory)

# Read the data for each paleo coord
info_per_pcoord = pd.read_csv(os.path.join(csv_directory, "info_per_pcoord_main.csv"))
info_per_pcoord['lat'] = np.NaN
info_per_pcoord['long'] = np.NaN
# split into lat and long for shapefile
for index, row in info_per_pcoord.iterrows():
	lat, long = row["pcoords"].split(",")
	info_per_pcoord.loc[index, 'lat'] = pd.to_numeric(lat) # plus some small modifier\n",
	info_per_pcoord.loc[index, 'long'] = pd.to_numeric(long)
info_per_pcoord["interval"] = [x.lower() for x in info_per_pcoord["interval"]]
total_per_interval = info_per_pcoord.groupby(["tetrapod_group"],
										   squeeze=True)["individuals_total"].sum().reset_index()
local_community_size = max(total_per_interval.individuals_total)


def generate_group_list():
	output = []
	job_type = 0
	for tet_group in tetrapod_groups:
		job_type += 1
		for seed in range(1, 11, 1):
			output.append({"tetrapod_group" : tet_group, "seed" : seed, "job_type" : job_type})
	output.sort(key=itemgetter("job_type", "seed"))
	return output

if __name__ == "__main__":

	params = generate_group_list()[int(sys.argv[1])]
	seed = params["seed"]
	job_type = params["job_type"]
	tet_group = params["tetrapod_group"]
	total_per_tet_group = int(total_per_interval[total_per_interval.tetrapod_group == tet_group].individuals_total)
	print("Starting {}: {}".format(tet_group, seed))
	sim = os.path.join(output_dir, "data_{}_{}.db".format(job_type, seed))
	# print("data_{}_{}.db".format(job_type, seed))
	tree = CoalescenceTree(sim, logging_level=40)
	# Drop existing tables
	tree._check_database()
	tree.cursor.execute("DROP TABLE IF EXISTS FRAGMENT_RICHNESS")
	tree.cursor.execute("DROP TABLE IF EXISTS BIODIVERSITY_METRICS")
	tree.cursor.execute("DROP TABLE IF EXISTS ALPHA_DIVERSITY")
	tree.cursor.execute("DROP TABLE IF EXISTS BETA_DIVERSITY")
	tree.database.commit()
	tmp_richness = []
	output_alpha = []
	output_richness = []
	# Loop over the speciation parameters (including the potential metacommunities)
	for reference in tree.get_community_references():
		params = tree.get_community_parameters(reference)
		metacommunity_reference = params["metacommunity_reference"]
		if metacommunity_reference == 0:
			meta_params = {"metacommunity_size": 0, "speciation_rate": 0.0}
		else:
			meta_params = tree.get_metacommunity_parameters(metacommunity_reference)
		sad = tree.get_species_abundances(reference=reference)
		individuals = []
		for species_id, abundance in sad:
			individuals.extend([species_id] * abundance)
		# Shuffle the list, repeated 10 times to randomly draw species richness values.
		tmp_tmp_richness = []
		all_gammas = []
		# print("Shuffling species abundances")
		for i in range(10):
			random.shuffle(individuals)
			index = 0
			for interval in intervals:
				interval_group_fragments = info_per_pcoord[
					(info_per_pcoord["interval"] == interval) & (info_per_pcoord["tetrapod_group"] == tet_group)]
				tmp_inds = 0
				for _, row in interval_group_fragments.iterrows():
					fragment_name = row["fragment_name"]
					no_inds = row["individuals_total"]
					tmp_inds += no_inds
					fragment_individuals = []
					fragment_individuals.extend(individuals[index:index + no_inds])
					index += no_inds
					rich = len(set(fragment_individuals))
					tmp_tmp_richness.append({"interval": interval, "tetrapod_group": tet_group, "iteration": i,
											 "fragment": fragment_name, "richness": rich})
				all_gammas.append({"interval": interval, "gamma": len(set(individuals[index - tmp_inds:index]))})
		for interval in intervals:
			interval_group_fragments = info_per_pcoord[
				(info_per_pcoord["interval"] == interval) & (info_per_pcoord["tetrapod_group"] == tet_group)]
			interval_sel = [x for x in tmp_tmp_richness if x["interval"] == interval]
			interval_richness = []
			for _, row in interval_group_fragments.iterrows():
				fragment_name = row["fragment_name"]
				sel = [x["richness"] for x in interval_sel if x["fragment"] == fragment_name]
				mean_richness = sum(sel) / len(sel)
				interval_richness.append(
					{"seed": seed, "reference": reference, "speciation_rate": params["speciation_rate"],
					 "metacommunity_size": meta_params["metacommunity_size"],
					 "metacommunity_speciation_rate": meta_params["speciation_rate"],
					 "interval": interval, "tetrapod_group": tet_group,
					 "fragment": fragment_name, "richness": mean_richness})

			r = [x["richness"] for x in interval_richness]
			interval_gammas = [x["gamma"] for x in all_gammas if x["interval"] == interval]
			gamma = sum(interval_gammas) / len(interval_gammas)
			alpha = sum(r) / len(r)
			beta = gamma / alpha
			output_alpha.append({"interval": interval, "tetrapod_group": tet_group,
								 "seed": seed, "reference": reference,
								 "speciation_rate": params["speciation_rate"],
								 "metacommunity_size": meta_params["metacommunity_size"],
								 "metacommunity_speciation_rate": meta_params["speciation_rate"],
								 "alpha": alpha, "beta": beta, "species_richness": gamma})
			tmp_richness.extend(interval_richness)
	output_richness.extend(tmp_richness)
	output_alpha_df = pd.DataFrame(output_alpha)
	output_alpha_df.to_csv(os.path.join(output_csv_directory, "nse_metrics1_{}_{}.csv".format(job_type, seed)))
	output_fragment_richness_df = pd.DataFrame(output_richness)
	# This is to reduce the file size
	output_fragment_richness_df = output_fragment_richness_df[output_fragment_richness_df.metacommunity_size > 0]
	output_fragment_richness_df.to_csv(os.path.join(output_csv_directory,
													"nse_fragment_richness1_{}_{}.csv".format(job_type, seed)))