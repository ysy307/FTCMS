#!/bin/zsh

set -e

# PETSc must be configured against the very mpi.h it will later be compiled
# with: petscsys.h refuses to compile when the I_MPI version it was configured
# with no longer matches the one on the include path. An oneAPI update is
# therefore enough to break every C file that includes a PETSc header, and the
# fix is to reconfigure here rather than to relax the check.

PETSC_VERSION="3.25.4"

ROOT_DIR=$(realpath .)
THIRD_PARTY_DIR="$ROOT_DIR/third_party"

PETSC_DIR="$THIRD_PARTY_DIR/petsc-src"
PETSC_ARCH="arch-linux-c-opt"

# Only the Intel toolchain is built: the rest of the project finds PETSc
# through pkg-config under third_party/.local/intel.
INSTALL_PREFIX="$THIRD_PARTY_DIR/.local/intel"

MAKE_NP=$(nproc)

if [ ! -d "$PETSC_DIR" ]; then
    echo "Error: PETSc source not found at $PETSC_DIR."
    echo "Clone it with: git clone -b v${PETSC_VERSION} https://gitlab.com/petsc/petsc.git $PETSC_DIR"
    exit 1
fi

for COMPILER in mpiicx mpiicpx mpiifx; do
    if ! command -v "$COMPILER" >/dev/null 2>&1; then
        echo "Error: $COMPILER was not found. Source the oneAPI environment first."
        exit 1
    fi
done

# The NVIDIA HPC SDK also puts an mpiexec on the path, and configure rejects an
# Open MPI launcher next to Intel MPI's MPICH-derived libraries. Take the one
# that ships beside the compiler wrappers rather than whatever the path finds.
MPIEXEC="$(dirname "$(command -v mpiifx)")/mpiexec"
if [ ! -x "$MPIEXEC" ]; then
    echo "Error: no mpiexec beside mpiifx at $MPIEXEC."
    exit 1
fi

echo "--- Installing PETSc ${PETSC_VERSION} ---"

cd "$PETSC_DIR"

export PETSC_DIR
export PETSC_ARCH

# The external packages already sit in $PETSC_ARCH/externalpackages, so a
# reconfigure reuses them and needs no network.
./configure \
    --configModules=PETSc.Configure \
    --optionsModule=config.compilerOptions \
    COPTFLAGS=-O3 \
    CXXOPTFLAGS=-O3 \
    FOPTFLAGS=-O3 \
    --prefix="$INSTALL_PREFIX" \
    --with-cc=mpiicx \
    --with-cxx=mpiicpx \
    --with-fc=mpiifx \
    --with-mpiexec="$MPIEXEC" \
    --with-debugging=0 \
    --with-make-np="$MAKE_NP" \
    --with-blaslapack-dir=/opt/intel/oneapi/mkl/2026.0 \
    --download-hypre \
    --download-mumps \
    --download-scalapack \
    --download-superlu_dist \
    --download-metis \
    --download-parmetis \
    --with-hdf5 \
    --download-hdf5

make -j"$MAKE_NP" all
make install

echo "--- Completed PETSc ${PETSC_VERSION} ---"
echo "PETSC_DIR=$INSTALL_PREFIX"
