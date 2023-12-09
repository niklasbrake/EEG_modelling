ID1=$(sbatch --array=24,50,10,14,8,7,43,28,26,54 --parsable synthetic_spikes.sh)
ID1=$(sbatch --array=7,54 --parsable synthetic_spikes.sh)
sbatch --depend=afterok:$ID1 --array=7,10,14,28,54 shuffle_synapses.sh






sbatch --array=7 synthetic_spikes.sh
sbatch --array=10 synthetic_spikes.sh
sbatch --array=14 synthetic_spikes.sh
sbatch --array=28 synthetic_spikes.sh
sbatch --array=54 synthetic_spikes.sh

scancel 38104361_50
scancel 38104362_10
scancel 38104363_14
scancel 38104364_8
scancel 38104365_43
scancel 38104366_28
scancel 38104367_26
scancel 38104368_54


sbatch --depend=afterok:38124127 --array=7,10,14,28,54 shuffle_synapses.sh