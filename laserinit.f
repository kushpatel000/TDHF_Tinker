c
c     ###################################################
c     ##         Written by Kush Patel - 11/19/17      ##
c     ###################################################
c     
c     ################################################################
c     ##                                                            ##
c     ##  subroutine laserinit -- checks for keywords and           ##
c     ##                          initializes the TD laser pulses   ##  
c     ##                                                            ##
c     ################################################################
c
c
c
      subroutine laserinit(nstep)
      use sizes
      use keys
      use laservars
      logical exist
      integer i, next, freeunit, lunit
      character*20 keyword
      character*120 record, string
      character*120 laserfile
      integer nstep
      complex*16, II
      
      uselaser = .false.

      do i=1,nkey
         next = 1
         record = keyline(i)
         call gettext (record, keyword,next)
         call upcase (keyword)
         if(keyword(1:9) .eq. 'USELASER ') then
            uselaser = .true.
         end if

         if(keyword(1:10) .eq. 'LASERFILE ') then
            allocate(lasershape(nstep))
            string = record(next:120)
            read(string,*,err=10,end=10) laserfile
            lunit = freeunit()
            II = (0.0d0,1.0d0)

            open(unit=lunit,file=filename,action='read')

            do i=1,nstep
               pulse(i) = (0.0d0,0.0d0)
            end do
            do i=1,nstep
               read(9,*,end=11) xx,yy
               pulse(i) = xx+II*yy
            end do
 11         continue
            

            
 10         continue
         end if
         
      end do
      




      end subroutine laserinit
