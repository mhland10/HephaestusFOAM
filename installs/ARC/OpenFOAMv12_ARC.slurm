#!/bin/bash
#SBATCH --job-name=openfoam_build
#SBATCH --output=openfoam_build.log
#SBATCH --error=openfoam_build.err
#SBATCH --time=04:00:00        # adjust based on cluster speed
#SBATCH --ntasks=1             # single process for building
#SBATCH --cpus-per-task=8      # adjust to number of cores for parallel build
#SBATCH --partition=compute2    # adjust to your HPC partition

# --- Load modules ---
module purge
module load gcc/12.4.0
module load cmake/3.24.2
module load openmpi/4.1.6

# --- Ensure OpenMPI paths are first ---
export PATH=/apps/openmpi/4.1.6/bin:$PATH
export LD_LIBRARY_PATH=/apps/openmpi/4.1.6/lib:$LD_LIBRARY_PATH
export MPI_LIB_DIR=/apps/openmpi/4.1.6/lib

# --- Source OpenFOAM environment ---
source $HOME/OpenFOAM-13/etc/bashrc

# --- Build Zoltan so OpenFOAM can use it ---
##cd $HOME/ThirdParty-13/Zoltan-3.90/build
##make clean
##../configure --prefix=$HOME/ThirdParty-13/platforms/linux64GccDPInt32/Zoltan-3.90 \
##             --libdir=$HOME/ThirdParty-13/platforms/linux64GccDPInt32/lib/openmpi-system \
##             --with-mpi=/apps/openmpi/4.1.6 \
##             --with-mpi-compilers=no \
##             --with-mpi-libdir=/apps/openmpi/4.1.6/lib \
##             --disable-zoltan-cppdriver \
##             --disable-shared --enable-static \
##             CFLAGS="-fPIC" \
##             CXXFLAGS="-fPIC -m64 -std=c++0"
##make -j $SLURM_CPUS_PER_TASK
##make install

# --- Move to OpenFOAM directory ---
cd $HOME/OpenFOAM-13

# --- Build OpenFOAM ---
echo "Starting OpenFOAM v13 build at $(date)"
./Allwmake -j $SLURM_CPUS_PER_TASK > build.log 2>&1
echo "OpenFOAM v13 build finished at $(date)"
