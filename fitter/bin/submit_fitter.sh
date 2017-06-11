#!/bin/bash


PARTITION=coms
TIME=48:00:00
NCPUS=4
NNODES=1

# Get bin dir
BIN=`realpath $0`
BIN=${BIN%/*}

JOB="ppqm-${PWD##*/}"
SUBMIT=qsub.tmp
PWD=`pwd`

cat > $SUBMIT <<!EOF
#!/bin/sh
#SBATCH --job-name=$JOB
#SBATCH --nodes=$NNODES
#SBATCH --cpus-per-task=$NCPUS
#SBATCH --ntasks=$NNODES
#SBATCH --error=$PWD/$JOB\_%j.stderr
#SBATCH --output=$PWD/$JOB\_%j.stdout
#SBATCH --time=$TIME
#SBATCH --partition=$PARTITION
#SBATCH --no-requeue

# Create scratch folder
mkdir /scratch/\$SLURM_JOB_ID
cd /scratch/\$SLURM_JOB_ID

# copy input files
cp $PWD/* .

# copy mndo.exe
cp /home/andersx/mndo_stuff/mndo99/mndo99_20121112_intel64_ifort-11.1.080_mkl-10.3.12.361 mndo
chmod +x mndo

\$PWD = `pwd`

# Run fitter
$BIN/fit.py -r "\$PWD/mndo < "


# Remove scratch folder
rm -rf /scratch/\$SLURM_JOB_ID

!EOF

sbatch $SUBMIT

