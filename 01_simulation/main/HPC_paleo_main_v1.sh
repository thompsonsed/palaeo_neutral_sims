#./bin.sh
#PBS -l walltime=71:50:00
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -J 0-399

module load intel-suite
module load python ## unfortunately can't be run with anaconda due to conflicts between sqlite versions in gdal.
module load gdal
module load boost
#module load anaconda
#module load sqlite # this isn't required - have no idea why.
echo $PBS_ARRAY_INDEX
echo $WORK
mkdir $TMPDIR/Data
python $WORK/Panama/Code/HPC_paleo_main_v1.py $PBS_ARRAY_INDEX
