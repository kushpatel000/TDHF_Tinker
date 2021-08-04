c
c
c     ###################################################
c     ##       Written by Kush Patel - 11/19/17        ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module laservars  --  contents of the laser pulse      ##
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
      module laservars
      implicit none
      logical uselaser

c     copy variables from other calculations
      real*8 laservec(3)
      real*8, allocatable :: lasershape(:)
      complex*16, allocatable :: cowp(:,:)
      save
      end
