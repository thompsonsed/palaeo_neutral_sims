#./bin.sh
#PBS -l walltime=71:50:00
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -J 0-399

module load intel-suite
module load gcc
module load anaconda3/personal
source activate SimulationEnv
echo $PBS_ARRAY_INDEX
mkdir $TMPDIR/Data
python $HOME/Paleo/HPC_paleo_main_v2.py $PBS_ARRAY_INDEX
