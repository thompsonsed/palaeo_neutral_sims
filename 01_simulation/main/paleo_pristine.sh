#./bin.sh
#PBS -l walltime=23:50:00
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -J 0-399

module load intel-suite
module load gcc
module load anaconda3/personal
source activate SimulationEnv

mkdir $TMPDIR/Data
python $HOME/Paleo/paleo_pristine.py $PBS_ARRAY_INDEX
