#./bin.sh
#PBS -l walltime=71:50:00
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -J 0-19

module load anaconda3/personal
module load intel-suite
module load gdal
module load boost
module load gcc
#module load anaconda
#module load sqlite # this isn't required - have no idea why.
python $HOME//rds/general/user/set114/home/Paleo/nse_sims_single.py $PBS_ARRAY_INDEX
