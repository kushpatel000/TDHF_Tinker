c
c     ###################################################
c     ##         Written by Kush Patel - 2/01/16       ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine tdhfinit  --  checks for keyword tdhf and    ##
c     ##                           and initializes necessary      ##
c     ##                           values                         ##
c     ##                                                          ##
c     ##############################################################
c
c     "tdhfinit" checks for the keyword "TDHF" which decides
c     whether the tchebychev propagator will be used or not
c
      subroutine tdhfinit
      use sizes
      use keys
      use tdhfvars
      logical exist, resumeq
      integer i, next, freeunit
      character*20 keyword
      character*120 record, string


      usetdhf = .false.
      tdhfdebug = .false.
      printq = .false.
      probcurr = .false.
      useurhf = .false.
      resumeq = .false.
      res1 = .false.
      outiter = 0
      outfreq = 1

      do i=1,nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if(keyword(1:5) .eq. 'TDHF ') then
            usetdhf  = .true.
            tdhffirst = .true.
         else if(keyword(1:7) .eq. 'TDHFCR ') then
            string = record(next:120)
            tdhfcr = 0.0d0
            read(string,*,err=10,end=10) tdhfcr
 10         continue
         else if(keyword(1:10).eq.'TDHFDEBUG ') then
            tdhfdebug = .true.
            debugout = freeunit()
c     open the debugging output
            inquire(file='tdhf.debug',exist=exist)
            if(exist) then
               open(unit=debugout,file='tdhf.debug',status='old')
               rewind(unit=debugout)
            else
               open(unit=debugout,file='tdhf.debug',status='new')
            end if
         else if(keyword(1:8).eq.'USEURHF ') then
            useurhf = .true.
         else if(keyword(1:11).eq.'PRINTEVERY ') then
            string = record(next:120)
            read(string,*,err=20,end=20) outfreq
            if(outfreq.lt.1) outfreq = 1
 20         continue
         else if(keyword(1:9).eq.'PROBCURR ') then
            probcurr = .true.
         else if(keyword(1:8).eq.'RESUME ') then
            resumeq = .true.
            res1 = .true.
         end if
      end do

      if(resumeq) call tdhfload()

      return
      end
