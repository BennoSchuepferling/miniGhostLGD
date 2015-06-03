! ************************************************************************
!
!               miniGhost: stencil computations with boundary exchange.
!                 Copyright (2012) Sandia Corporation
!
! Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
! license for use of this work by or on behalf of the U.S. Government.
!
! This library is free software; you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation; either version 2.1 of the
! License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
! USA
! Questions? Contact Richard F. Barrett (rfbarre@sandia.gov) or
!                    Michael A. Heroux (maherou@sandia.gov)
!
! ************************************************************************

MODULE DRIVER_LGD_MOD

CONTAINS

   SUBROUTINE DRIVER_LGD ( IERR )

   ! Purpose
   ! =======
   ! Driver for the bulk synchronous parallel with message aggregation 
   ! communication strategy. Called by the Fortran driver subroutine MINI_GHOST.
   
   ! no more BSPMA or SVAF - its just for LibGeoDecomp 
   ! Author : Benno Schuepferling
   
   USE MG_CONSTANTS_MOD
   USE MG_STENCIL_MOD
   USE MG_PROFILING_MOD
   USE MG_SUM_GRID_MOD

   IMPLICIT NONE

      INTEGER, INTENT(OUT) :: &
         IERR                       ! Return status

      ! ------------------
      ! Local Declarations
      ! ------------------
   
      INTEGER ::                  &
         IVAR,                    &  ! Counter (over variables)
         ISPIKE                      ! Counter (over NUM_SPIKES)
   
      REAL(KIND=MG_REAL) ::       &
         ERROR_ITER,              &  ! Difference between new and old GRIDi sum.
         GSUM                        ! Global SUM across GRIDs.

      REAL(KIND=MG_REAL8) ::      &
         TIME_START,              &
         TIME_START_2,            &
         TIME_START_ALL
   
      ! ---------------------
      ! Executable Statements
      ! ---------------------
   
      IERR = 0

      CALL simulate(GLOBAL_NX, GLOBAL_NY, GLOBAL_NZ, DEBUG_GRID, STENCIL, NUM_VARS, NUM_SPIKES, NUM_TSTEPS, ERROR_TOL, SOURCE_TOTAL, SPIKES, SPIKE_LOC, GRIDS_TO_SUM)

   END SUBROUTINE DRIVER_LGD

END MODULE DRIVER_LGD_MOD
