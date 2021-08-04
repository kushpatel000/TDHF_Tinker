c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module civars  --  contents of the CI singlet and      ##
c     ##                     triplet  calculations               ##
c     ##                                                         ##
c     #############################################################
c
c
c     maxkey    maximum number of lines in the keyword file
c
c     nkey      number of nonblank lines in the keyword file
c     keyline   contents of each individual keyword file line
c
c
      module civars
      implicit none
      logical usetriplet,usesinglet
      integer ciout

c     copy variables from other calculations
      integer ci_nfill, iex
      real*8, allocatable :: ci_gamma(:,:)
      real*8, allocatable :: ci_v(:,:)
      real*8, allocatable :: ci_en(:)
      real*8, allocatable :: ci_ed(:,:)
      real*8, allocatable :: ci_hc(:,:)
      save
      end
