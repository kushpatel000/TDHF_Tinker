c     1;4205;0c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module tdhfvars  --  contents of the PITDHF            ##
c     ##                       calculation                       ##
c     ##                                                         ##
c     #############################################################
c
c
c     maxkey    maximum number of lines in the keyword file
c
c     nkey      number of nonblank lines in the keyword file
c     keyline   contents of each individual keyword file line
c
c     usetdhf     logical that indicates whether TDHF is used
c     tdhffirst   logical switch that indicates whether the dynamics
c                 are in the first time step
c     tdhfdebug   logical for deciding whether to print debugging
c                 information
c     tdhfdt      time step for dynamics
c     tdhfout     output file unit number
c     tdhf_ nfill fill level for the orbitals (# of electrons / 2)
c     tdhfeda     tdhf electron density matrix (spin up)
c     tdhfedb     tdhf electron density matrix (spin down)
c     tdhffocka   tdhf fock matrix (spin up)
c     tdhffockb   tdhf fock matrix (spin down)
c     tdhfgamma   tdhf gamma matrix
c     tdhfhc      tdhf core matrix
c     bonded      logical that indicates whether a pair of atoms
c                 are bonded
c     tdhfcr      coupling radius for the exponetial scaling term (Angstrom)

      
      module tdhfvars
      implicit none
      logical usetdhf, tdhffirst, tdhfdebug, printq
      logical probcurr, useurhf, res1
      integer tdhfout, debugout, tdhf_nfill, outiter
      integer outfreq
      real*8  tdhfdt, superevmin, superevmax
      complex*16, allocatable :: tdhfeda(:,:)
      complex*16, allocatable :: tdhfedb(:,:)
      complex*16, allocatable :: tdhffocka(:,:)
      complex*16, allocatable :: tdhffockb(:,:)
      complex*16, allocatable :: tdhfgamma(:,:)
      complex*16, allocatable :: tdhfhc(:,:)

c-----Test Vars 11.03.16-----
      logical bonded
      real*8  tdhfcr
c----------------------------

c-----Test Vars 11 Jan 18----      
      complex*16, allocatable :: tdhfed(:,:)
      complex*16, allocatable :: tdhffock(:,:)

      save
      end
