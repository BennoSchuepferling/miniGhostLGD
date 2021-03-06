#@header
# ************************************************************************
#
#      miniGhost: stencil computations with boundary exchange.
#              Copyright (2012) sandia corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Richard F. Barrett (rfbarre@sandia.gov) or
#                    Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
#@header

# Simple hand-tuned makefile, configured for the gnu compiler and LibGeoDecomp.
# Modify as necessary for your environment.

#-----------------------------------------------------------------------

# Option: -D_MG_MPI
PROTOCOL = -D_MG_MPI
#PROTOCOL = -D_MG_SERIAL
# Compilers
FC=mpifort
#FC=ftn
CC=mpicc
#CC=cc
CXX=g++

CFLAGS = $(PROTOCOL)
# C main calling Fortran subroutine:
CFLAGS += -Df2c_

FFLAGS = $(PROTOCOL)

# Variable precision: -D_INT8 and/or -D_REAL8.
FFLAGS += -D_MG_INT4 -D_MG_REAL8
# Optimization
# Next line PGI compiler
#OPT_F = -fast -fastsse
#OPT_F = -g
#FFLAGS += $(OPT_F)

# Free-form Fortran source code:
# Next line PGI compiler:
#FFLAGS += -Mfree
# Next line Intel compiler:
#FFLAGS += -free -w
# Next line Cray CCE compiler:
#FFLAGS += -f free -m 4
# Next line gfortran
FFLAGS += -ffree-form

#disables line length restrictions which somehow can give us trouble compiling
FFLAGS += -ffree-line-length-0

# Array bounds checking: (expensive!)
#FFLAGS += -Mbounds
# Compile to include checkpointing capability.
#FFLAGS += -D_MG_CHECKPT

LD = $(CXX)
LDFLAGS = $(CFLAGS)

CPPFLAGS = -std=c++11  -Wall -Wno-sign-promo -Wnon-virtual-dtor -march=native -fopenmp -O3
LIBS = -lgfortran 

# linking with g++ so we need this for fortran libs
LIBS += $(shell mpifort --showme:link)

# get dependencies for lgd from pkg-config 
# just expecting people to build lgd in ../libgedecomp/build/
CPPFLAGS += $(shell pkg-config --cflags ../libgeodecomp/build/libgeodecomp.pc)
LIBS += $(shell pkg-config --libs ../libgeodecomp/build/libgeodecomp.pc)


include make_targets

# End makefile

