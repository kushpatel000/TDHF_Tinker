c
c     ###################################################
c     ##         Written by Kush Patel - 2/29/16       ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ciinit  -- checks for keyword "EXCITES" or   ##
c     ##                        "EXCITET" for single or triplet   ##
c     ##                        calculations respectively, then   ##
c     ##                        initializes necessary values      ##
c     ##                                                          ##
c     ##############################################################
c
c     "ciinit" checks for the keyword "EXCITES" or "EXCITET" which
c     decides whether to do excited state calculations
c        'EXCITES' - singlet
c        'EXCITET' - triplet

      subroutine ciinit
      use sizes
      use keys
      use civars
      logical exist
      integer i, next, freeunit
      character*20 keyword
      character*120 record, string
      character*120 cifile

      usesinglet = .false.
      usetriplet = .false.
      do i=1,nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if(keyword(1:8) .eq. 'EXCITES ') then
            usesinglet = .true.
            usetriplet = .false.
            ciout = freeunit()
            cifile = 'ci.sing.debug'
            string = record(next:120)
            iex = 0
            read(string,*,err=10,end=10) iex
 10         continue
         else if(keyword(1:8) .eq. 'EXCITET ') then
            usesinglet = .false.
            usetriplet = .true.
            ciout = freeunit()
            cifile = 'ci.trip.debug'
            string = record(next:120)
            iex = 0
            read(string,*,err=11,end=11) iex
 11         continue
         else if(keyword(1:7) .eq. 'GROUND ') then
            usesinglet = .false.
            usetriplet = .false.
            iex = 0
 12         continue
         end if
      end do

c     open the CI debug file
      if(usesinglet.or.usetriplet) then
         inquire(file=cifile,exist=exist)
         if(exist) then
            open(unit=ciout,file=cifile,status='old')
            rewind(ciout)
         else
            open(unit=ciout,file=cifile,status='new')
         end if
      end if
      
      if(usesinglet.or.usetriplet) then
         write(ciout,*) 'ci sing: ', usesinglet
         write(ciout,*) 'ci trip: ', usetriplet
      else
         iex = 0
      end if

      return
      end
