#!/bin/bash

pat='[^0-9]+([0-9]+)'

rm job* std* *hzd* out*
rm -rf hzd*/

# job_scf.sh
cat << EOF > job_scf.sh
#!/bin/bash
export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}

module use \$OTHERSTAGES
module load Stages/2020
module load Intel/2020.2.254-GCC-9.3.0
module load ParaStationMPI/5.4.7-1
module load QuantumESPRESSO/6.6

srun pw.x < inp_scf > out_scf 
EOF

# do the scf simulation
job_scf=$(sbatch --job-name=scf \
	   --output=stdout --error=stderr \
	   --partition=batch --account=fzj-mac \
	   --nodes=1 --ntasks-per-node=48 --time=10:00:00 \
	   job_scf.sh)

[[ $job_scf =~ $pat ]]

# do the nscf calculation

cat << EOF > job_nscf.sh
#!/bin/bash
export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}

module use \$OTHERSTAGES
module load Stages/2020
module load Intel/2020.2.254-GCC-9.3.0
module load ParaStationMPI/5.4.7-1
module load QuantumESPRESSO/6.6

srun pw.x < inp_nscf > out_nscf 
EOF

job_nscf=$(sbatch --dependency=afterok:${BASH_REMATCH[1]}\
	   --partition=batch --account=fzj-mac \
	   --job-name=nscf --output=stdout_nscf --error=stderr_nscf \
	   --nodes=1 --ntasks-per-node=48 --time=10:00:00 \
	   job_nscf.sh)

[[ $job_nscf =~ $pat ]]

# inp_ldos 
cat << EOF > inp_ldos 
&projwfc
 outdir='./',
 lsym=.true.
 prefix='hzd',
 filpdos='hzd',
 filproj='hzd.proj.dat',
 ngauss=0, degauss=0.015,
 Emin=-40.0, Emax=40.0, DeltaE=0.01
/
EOF

# job_ldos.sh
cat << EOF > job_ldos.sh
#!/bin/bash

module use \$OTHERSTAGES
module load Stages/2020
module load Intel/2020.2.254-GCC-9.3.0
module load ParaStationMPI/5.4.7-1
module load QuantumESPRESSO/6.6

projwfc.x < inp_ldos > out_ldos 
EOF

# calculate the pdos
job_ldos=$(sbatch --dependency=afterok:${BASH_REMATCH[1]} --job-name=pdos\
	       --output=stdout_pdos --error=stderr_pdos\
		   --partition=batch --account=fzj-mac\
		   --nodes=1 --ntasks-per-node=48 --time=00:30:00 job_ldos.sh)