c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine cicalc  --  Calculates CI energies and states    ##
c     ##                                                              ##
c     ##     Originally written by Bittner(?) for an older version    ##
c     ##     of Tinker. Rewritten by Kush Patel for Tinker 7.1.2      ##
c     ##                                                              ##
c     ##################################################################

      subroutine cicalc

      use sizes
      use atoms
      use bndstr
      use civars
      use iounit
      use piorbs
      use units
c----------------------------------
c     to see if TDHF is requested
      use tdhfvars
c----------------------------------
      
      implicit none
      
c     cicalc variables
c     index variables
      integer i,j,k,l,m,ma,iip,iorb,jorb,iir
      integer ih,ie
      
      integer nv,nc,nci,nci_max,nl,nu
c        nv  - num valence orbitals 
c        nc  - num condunting orbitals
c        nci - total CI interactions
c        nci_max - maximum calculated CI interactions
c        nl  - lower limit (for systems w/ <10 oribtals)
c        nu  - upper limit (for systems w/ <10 orbitals)
      parameter(nv = 20, nc = 20, nci_max = nv*nc)
      integer iconfig(nci_max,2)
c        keeps track of configurations
c        (n-9,n+1)
c        (n-9,n+2)
c           ...
c        (n-9,n+10)
c        (n-8,n+1)
c           ...
c        (n,n+10)
c      real*8 hfoc(norbit)
c      real*8 exciton(norbit),tranvect(norbit)
c      real*8 hodens(norbit),eldens(norbit)
c      real*8 nco1(norbit,norbit),nco2(norbit,norbit)
c      real*8 nco3(norbit,norbit),nco4(norbit,norbit)
      real*8 nhf(norbit),nh(norbit),ne(norbit)
      real*8 dens(norbit,norbit)
      real*8 cibcm(norbit,norbit)
      real*8 hcst(norbit,norbit),gst(norbit,norbit)

      
c     Summation and Integral Variables
      real*8 qx(nci_max),qy(nci_max),qz(nci_max)
      real*8 qxsum,qysum,qzsum,Jint,Kint
      real*8 s1,s2,qtot,fosci,sume,sumh
      real*8 sum1,sum2,sum3,sum4
      real*8 sum,p

c     CI Hamiltonian matrices
      real*8, allocatable :: sing(:,:)
      real*8, allocatable :: trip(:,:)
      real*8, allocatable :: A1(:,:),A2(:,:)

      
c     LAPACK/BLAS work variables
      real*8, allocatable :: d(:)
      real*8, allocatable :: work(:)
      real*8, allocatable :: evec(:,:)
      integer lwork,info
c        d    - eigenvalues of CI hamiltonian
c        work - work variable

c     Density matrix (State Rep)
      real *8 stds(norbit,norbit)
c     New Coherences variables
      integer r,s,o,kd

      nl = max(ci_nfill - nv + 1, 1)
      nu = min(ci_nfill + nc, norbit)
      nci = (ci_nfill - nl + 1)*(nu - ci_nfill)
      lwork = 32*nci
      
 307  format(X,A9,I3)
      write(ciout,307) 'ci_nfill: ', ci_nfill
      write(ciout,307) 'nl:       ', nl
      write(ciout,307) 'nu:       ', nu
      write(ciout,307) 'norbit:   ', norbit
      write(ciout,307) 'nci:      ', nci
      write(ciout,307) 'iex:      ', iex
      write(ciout,*)
      
      if(iex.lt.1 .or. iex.gt.nci) then
         write(ciout,*) 'Invalid excitation number' 
         call fatal
      end if

c     Allocate orbitals according to the CI interactions
      allocate (sing(nci,nci))
      allocate (trip(nci,nci))
      allocate (d(nci))
      allocate (work(lwork))
      allocate (evec(nci,nci))

c-----Test Print 10/8/16----------
      write(ciout,*) 'Input Energies (ci_en)'
      do i=1,norbit
         write(ciout,*) i, ci_en(i)
      end do
c-----End Test Print--------------      
      
c     sum over molecular orbitals
      m = 0
      do i = nl, ci_nfill
         do j=ci_nfill+1, nu
            m = m+1
            iconfig(m,1) = i
            iconfig(m,2) = j

            qxsum = 0.0d0
            qysum = 0.0d0
            qzsum = 0.0d0
c     transitions before CI
            qx(m) = 0.0d0
            qy(m) = 0.0d0
            qz(m) = 0.0d0

            do iip = 1,norbit
               iorb = iorbit(iip)
               qxsum = qxsum + ci_v(iip,i)*ci_v(iip,j)*x(iorb)
               qysum = qysum + ci_v(iip,i)*ci_v(iip,j)*y(iorb)
               qzsum = qzsum + ci_v(iip,i)*ci_v(iip,j)*z(iorb)
            end do
            qx(m) = qxsum
            qy(m) = qysum
            qz(m) = qzsum

            ma = 0

c     construct CI Hamiltonian for singlet and triplet
            do k = nl,ci_nfill
               do l = ci_nfill+1,nu
                  ma = ma + 1
                  Jint = 0.0d0
                  Kint = 0.0d0

                  do iip=1,norbit
                     s1 = 0.0d0
                     s2 = 0.0d0
                     do iir = 1,norbit
                        s1 = s1 + ci_gamma(iip,iir)*ci_v(iir,j)*
     >                       ci_v(iir,l)
                        s2 = s2 + ci_gamma(iip,iir)*ci_v(iir,k)*
     >                       ci_v(iir,l)
                     end do
                     Jint = Jint + ci_v(iip,i)*ci_v(iip,k)*s1
                     Kint = Kint + ci_v(iip,i)*ci_v(iip,j)*s2
                  end do
                  sing(m,ma) = -Jint + 1.0d0*Kint
c     Original    sing(m,ma) = -Jint + 2.0d0*Kint
                  sing(ma,m) = sing(m,ma)
                  trip(m,ma) = -Jint
                  trip(ma,m) = -Jint
               end do
            end do
            sing(m,m) = sing(m,m) + (ci_en(j)-ci_en(i))
            trip(m,m) = trip(m,m) + (ci_en(j)-ci_en(i))
         end do
      end do
c-----Test Print--------
c     Need to see what all the values in iconfig mean
      write(ciout,*) "  m     i   j"
 999  format(i5,i5,' -> ',i5)
      do i=1,m
         write(ciout,999) i, iconfig(i,1), iconfig(i,2)
      end do
      write(ciout,*) ''
c-----End Test Print----

c-----Test Code 04 Jan 18------------------------------------------------
c     Testing my calculation of the Tamm-Dankoff Approximation

c (1) Need to find convert Hartree Core and Gamma matrices to state
c     representation
c      do i=1,norbit
c         do j=1,norbit
c            s1 = 0.0d0
c            s2 = 0.0d0
c            do k=1,norbit
c               do l=1,norbit
c                  s1 = s1 + ci_v(k,i)*ci_v(l,j)*ci_gamma(k,l)
c                  s2 = s2 + ci_v(k,i)*ci_v(l,j)*ci_hc(k,l)
c               end do
c            end do
c            gst(i,j)  = s1
c            hcst(i,j) = s2
c         end do
c      end do

c (2) Generate the matrices A1 and A2, then copy them to sing
c      allocate(A1(nci,nci))
c      allocate(A2(nci,nci))
c      do m=1,nci
c         ih = iconfig(m,1)
c         ie = iconfig(m,2)
c         do ma=1,nci
c            A1(m,ma) = 0.0d0
c            iir = iconfig(ma,1)
c            iip = iconfig(ma,2)
c            s1 = 0.0d0
c            if(m.eq.ma) then
c               do i = 1,nl,ci_nfill
c                  s1 = s1 + hcst(i,i)
c               end do
c            end if
c            if(iir.eq.ih) s1 = s1 + hcst(ie,iip)
c            if(iip.eq.ie) s1 = s1 - hcst(ih,iir)
c            A1(m,ma) = s1
c         end do
c      end do

c      do m=1,nci
c         do ma=1,nci
c            A2(m,ma) = 0.0d0
c         end do
c         ih = iconfig(m,1)
c         ie = iconfig(m,2)
c         s1 = 0.0d0
c         s2 = 0.0d0
c         do i=nl,ci_nfill
c            do j=1,nl,ci_nfill
c               s1 = s1 + gst(i,j)
c            end do
c            s2 = s2 + gst(i,ie) + gst(ie,i) - gst(i,ih) - gst(ih,i)
c         end do
c         sum = s1 + s2 - gst(ih,ie) - gst(ie,ih) + 2.0d0*gst(ih,ih)
c         A2(m,m) = 0.5d0*sum
c      end do

c      do m=1,nci
c         do ma=1,nci
c            sing(m,ma) = 0.0d0 - A1(m,ma) - A2(m,ma)
c         end do
c      end do
c - well, lets give it a shot...      

c      write(ciout,*) 'Hcore matrix (state)'
c      do i=1,norbit
c         write(ciout,306) (hcst(i,j),j=1,norbit)
c      end do
c      write(ciout,*) ''
      
c      write(ciout,*) 'Gamma matrix (state)'
c      do i=1,norbit
c         write(ciout,306) (gst(i,j),j=1,norbit)
c      end do
c      write(ciout,*) ''

c------------------------------------------------------------------------
      
c
c     Print the CI Hamiltonians
c
 302  format(E13.5)
      if(usesinglet) then
         write(ciout,*) "CI Singlet Hamiltonian"
         do i=1,nci
            do j=1,nci
               write(ciout,302,advance='no') sing(i,j)
            end do
            write(ciout,*) ''
         end do   
      end if
      if(usetriplet) then
         write(ciout,*) "CI Triplet Hamiltonian"
         do i=1,nci
            do j=1,nci
               write(ciout,302,advance='no') trip(i,j)
            end do
            write(ciout,*) ''
         end do
      end if
      write(ciout,*)
c      
c     Get the eigensystem, print
c
      if(usesinglet) then
         call dsyev('V','U',nci,sing,nci,d,work,lwork,info)
c     At this point 'sing' holds the eigenvectors
c     Copy them to 'evec'
         do i=1,nci
            do j=1,nci
               evec(i,j) = sing(i,j)
            end do
         end do
         
         write(ciout,*) 'CI Singlet Energies (ev)'
      end if
      if(usetriplet) then
         call dsyev('V','U',nci,trip,nci,d,work,lwork,info)
c     At this point 'trip' holds the eigenvectors
c     Copy them to 'evec'
         do i=1,nci
            do j=1,nci
               evec(i,j) = trip(i,j)
            end do
         end do

         write(ciout,*) 'CI Triplet Energies (ev)'
      end if

c     Print out the energies
 303  format(i5,2f12.6)
      do i=1,min(20,nci)
         write(ciout,303) i,d(i)*evolt
      end do
      write(ciout,*)
c
c     Analyse CI states
c
      if(usesinglet) then
         write(ciout,*) 'Analysis of Singlet States'
      else if(usetriplet) then
         write(ciout,*) 'Analysis of Triplet States'
      end if

c      do i=1,min(20,nci)
      do i=1,nci
         qxsum = 0.0d0
         qysum = 0.0d0
         qzsum = 0.0d0

         do j=1,nci
            qxsum = qxsum + evec(j,i)*qx(j)
            qysum = qysum + evec(j,i)*qy(j)
            qzsum = qzsum + evec(j,i)*qz(j)
         end do

         qtot = qxsum*qxsum + qysum*qysum + qzsum*qzsum
         fosci = 0.0875161*qtot*d(i)
         if(fosci.lt.10.0d0**(-90)) fosci = 0.0d0

         write(ciout,304) i,d(i)*evolt,fosci
 304     format('E(',i3,') = ',f12.6,2x,'fosci = ',e12.6)
         write(ciout,305) qxsum,qysum,qzsum,sqrt(qtot)
 305     format('qx = ',e12.6, 2x, 'qy = ',e12.6, 2x, 'qz = ',e12.6,
     >        2x,'qtot = ',e12.6)
               
         write(ciout,*) 'v -> c          Ec-Ev          coef'
         do j=1,nci
            if(evec(j,i)**2.gt.0.1) then
               write(ciout,'(i3," ->",i3,3x,f12.6,3x,f12.6)')
     >              iconfig(j,1),iconfig(j,2),
     >              (ci_en(iconfig(j,2))-ci_en(iconfig(j,1)))*evolt,
     >              evec(j,i)
            end if
         end do
         write(ciout,*)
      end do

c      
c     Modify electron densitites for excitations.
c     He were use the assumption that the CI(S)
c     state add and subtract single electron
c     density from the ci_ed(i,j) matrix
c
c     May need iex here.
c
c     First construct the density matrix using
c     the CI coefficients for iex in MO basis
c
      sum1 = 0.0d0
      do i=1,norbit
         do j=1,norbit
            dens(i,j) = 0.0d0
         end do
      end do

      do m=1,nci
         i = iconfig(m,1)
         j = iconfig(m,2)
         dens(i,j) = evec(m,iex)
      end do

      write(ciout,*) 'CI Eigenvectors'
      do i=1,nci
         write(ciout,306) (evec(i,j),j=1,nci)
      end do
      write(ciout,*) 

 306  format(8f12.6)
      write(ciout,*) 'CI expansion coefficients'
      do i=1,norbit
         write(ciout,306) (dens(i,j),j=1,norbit)
      end do
      write(ciout,*)
c     At this point, we have the CI density matrix in the local orbital basis.
c
c
c
c     Now we need to write the populations in the orbital basis. For this, we
c     first determine the orbital population changes
c
      do i=1,norbit
         nhf(i) = 0.0d0
         nh(i)  = 0.0d0
         ne(i)  = 0.0d0
         if(i.le.ci_nfill) then
            if(usetdhf) then
               nhf(i) = 1.0d0
            else
               nhf(i) = 2.0d0
            end if
         end if 
      end do

      do ih = nl,ci_nfill
         sum = 0.0d0
         do ie = ci_nfill+1,nu
            sum = sum + dens(ih,ie)**2
         end do
         nh(ih) = sum
      end do
      do ie=ci_nfill+1,nu
         sum = 0.0d0
         do ih=nl,ci_nfill
            sum = sum + dens(ih,ie)**2
         end do
         ne(ie) = sum
      end do

c     Check: write the occupations over the range nl->nu
 308  format(I3,4F13.5)
      write(ciout,*) "State Populations Post-CI over CI range"
      write(ciout,*) " k      nhf          nh           ne          adj"
      do k=nl,nu
         write(ciout,308) k,nhf(k),nh(k),ne(k),nhf(k)+ne(k)-nh(k)
      end do

c     This is the CI modified bond-charge matrix
c      do i=1, norbit
c         do j=1, norbit
c            cibcm(i,j) = 0.0d0
c            sum = 0.0d0
c            do k=1, norbit
c               sum = sum + (nhf(k)+ne(k)-nh(k))*ci_v(i,k)*ci_v(j,k)
c            end do
c            cibcm(i,j) = sum
c         end do
c      end do
c      write(ciout,*)



c-----Test Code 4-18-17-------------------------------------
c     Setting the state representation density matrix
c     equal to state populations on the diagonal and the CI
c     eigenvector elements as the coherences
c      do i=1, norbit
c         do j=1, norbit
c            stds(i,j) = 0.0d0
c            if(i.eq.j) stds(i,i) = nhf(i)+ne(i)-nh(i)
c         end do
c      end do
c      do m=1, nci
c         i = iconfig(m,1)
c         j = iconfig(m,2)
c         stds(i,j) = evec(m,iex)
c         stds(j,i) = evec(m,iex)
c      end do
c      write(ciout,*) 'State Rep Density Matrix'
c      do i=1, norbit
c         write(ciout,306) (stds(i,j),j=1,norbit)
c      end do
c      write(ciout,*) ''
      
c      do i=1, norbit
c         do j=1, norbit
c            sum = 0.0d0
c            do k=1, norbit
c               do l=1, norbit
c                  sum = sum + stds(k,l)*ci_v(i,k)*ci_v(j,l)
c               end do
c            end do
c            cibcm(i,j) = sum
c         end do
c      end do
c-----Notes 13 Nov 2017-------------------------------------
c     This definition of coherence causes systems to exhibit
c     non physical behavior, namely that there are
c     populations greater than 1 and less than 0.
c     Very easy to see in the 2 level ethylene case.
c     I am temporarily turning it off.
c-----------------------------------------------------------

c-----Test Code 20 Nov 2017---------------------------------
c     The state populations might be inherently included
c     in the CI eigenvectors. So lets take the eigenvector
c     and use that to populate coherences and diagonals.

c     start with the state representation
      do o=1, norbit                    ! state  
         do m=1, norbit                 ! state
            sum1 = 0.0d0
            stds(o,m) = 0.0d0
            do i=nl,ci_nfill            ! hole 1
               do j=ci_nfill+1,nu       ! elec 1
                  do k=nl,ci_nfill      ! hole 2
                     do l=ci_nfill+1,nu ! elec 2
                        kd = 0
                        if(l.eq.o .and. i.eq.k .and. m.eq.j) then
                           kd = kd+1
                        end if
                        if(l.eq.j .and. k.eq.m .and. o.eq.i) then
                           kd = kd-1
                        end if
                        if(l.eq.j .and. i.eq.k .and. o.eq.m) then
                           if(o.le.ci_nfill) kd = kd+1
                        end if
                        sum1= sum1 + float(kd)*dens(i,j)*dens(k,l)
                     end do
                  end do
               end do
            end do
            stds(o,m) = sum1
         end do
      end do

c     unitary transform into site represenation
      do r=1,norbit
         do s=1,norbit
            sum = 0.0d0
            do o=1,norbit
               do m=1,norbit
                  sum = sum + stds(o,m)*ci_v(r,o)*ci_v(s,m)
               end do
            end do
            cibcm(r,s) = sum
         end do
      end do
      
      write(ciout,*) 'Density Matrix (State Rep)'
      do i=1,norbit
         do j=1,norbit
            write(ciout,309,advance='no') stds(i,j)
         end do
         write(ciout,*) ''
      end do
      write(ciout,*) ''      
c-----------------------------------------------------------                  

      
c---------------------------------------------------------------------
c     This is Allen's code here
c---------------------------------------------------------------------      
c
c     create electron and hole in orbital representation
c
c      do i=1,norbit
c         sume = 0.0d0
c         sumh = 0.0d0
c         do j=1,norbit
c            sume = sume + dens(j,i)**2
c            sumh = sumh + dens(i,j)**2
c         end do
c         el(i) = sume
c         ho(i) = sumh
c      end do
c
c      write(iout,*) 'molecular orbitals'
c      do i=1,norbit
c         write(iout,302) (ci_v(i,j),j=1,norbit)
c      end do
c
c      write(iout,*) 'dens(i,j)'
c      do i=1,norbit
c         write(iout,302) (dens(i,j),j=1,norbit)
c      end do
c
c      write(iout,*) 'nci = ', nci
c      write(iout,*) 'nv  = ', nv
c      write(iout,*) 'nc  = ', nc
c
c      do i=1,norbit
c         sum1 = 0.0d0
c         do j=1,norbit
c            do k=1,norbit
c               sum1 = sum1 + ci_v(i,k)*ci_v(j,k)
c            end do
c         end do
c         hfoc(i) = sum1
c      end do
c      sum1 = 0.0d0
c
c      write(iout,*) 'hfoc (orbital density?)'
c      do i=1, norbit
c         write(iout,*) hfoc(i)
c      end do
c
c      write(iout,*) 'hole orbital'
c      do i=1,norbit
c         write(iout,*) ho(i)
c      end do
c
c      write(iout,*) 'electron orbital'
c      do i=1,norbit
c         write(iout,*) el(i)
c      end do
c
c      do i=1,norbit
cc        calc updated density
c         tranvect(i) = hfoc(i)-ho(i)+el(i)
cc        calc excition density
c         exciton(i)  = el(i)-ho(i)
cc        calc electron density
c         eldens(i)   = hfoc(i)+el(i)
cc        calc hole density
c         hodens(i)   = hfoc(i)-ho(i)
c      end do
c
cc     change the bond charge matrix (tranvect) to from
cc     orbital represenation to site represenation
c
c      do i=1,norbit
c         do j=1,norbit
c            sum1 = 0.0d0
c            sum2 = 0.0d0
c            sum3 = 0.0d0
c            sum4 = 0.0d0
c            do k=1,norbit
c               sum1 = sum1 + tranvect(k)*ci_v(i,k)*ci_v(j,k)
c               sum2 = sum2 + exciton(k)*ci_v(i,k)*ci_v(j,k)
c               sum3 = sum3 + eldens(k)*ci_v(i,k)*ci_v(j,k)
c               sum4 = sum3 + hodens(k)*ci_v(i,k)*ci_v(j,k)
c            end do
c            nco1(i,j) = sum1
c            nco2(i,j) = sum2
c            nco3(i,j) = sum3
c            nco4(i,j) = sum4
c         end do
c      end do
c
cc     Center of Mass calculations (in some weird way)
c      sum3 = 0.0d0
c      sum4 = 0.0d0
c      do i=1,norbit
c         sum3 = sum3 + (nco3(i,i)-1)*i
c         sum4 = sum4 + abs((nco4(i,i)-1)*i)
c      end do
c
c      write(iout,*) 'Center of masses'
c      write(iout,*) 'Electron: ', sum3
c      write(iout,*) 'Hole:     ', sum4
c-----------------------------------------------------------------------


c     Print out the CI density/bond-charge matrix
 309  format(F13.5)
      write(ciout,*) 'Density/Bond-Charge Matrix'
      do i=1,norbit
         do j=1,norbit
            write(ciout,309,advance='no') cibcm(i,j)
         end do
         write(ciout,*) ''
      end do
      write(ciout,*) ''
      
c     As a check, take a look at the changes in the bond-orders.
c     These will go in to the forcefield calculation next.
      write(ciout,*) 'bond   i  j  old new'
      do k=1, nbpi
         p = pnpl(k)
         i = ibpi(2,k)
         j = ibpi(3,k)
         if(usetdhf) then
            if(useurhf) then
               pnpl(k) = realpart(cibcm(i,j) + tdhfedb(i,j))
            else
               pnpl(k) = realpart(cibcm(i,j) + 0.5d0*tdhfed(i,j))
            end if
c-----Test 03.20.17----
            pbpl(k) = pnpl(k)*ci_hc(i,j)/(-0.0757d0)
c----------------------            
         else
            pnpl(k) = cibcm(i,j)
c-----Test 03.20.17----
            pbpl(k) = cibcm(i,j)*ci_hc(i,j)/(-0.0757d0)
c----------------------
         end if
         i = ibnd(1,ibpi(1,k))
         j = ibnd(2,ibpi(1,k))
         write(ciout,'(5i5,2f12.2)') k,ibpi(2,k),ibpi(3,k),i,j,p,pnpl(k)
      end do
      write(ciout,*)
      

c------------------------------------------------------------------------
c     If needed, update the TDHF electron density matrix
      if(usetdhf) then
         write(ciout,*)
         write(ciout,*) 'Updating TDHF electron density matrix'
         write(ciout,*)
         do i=1, norbit
            do j=1, norbit
               if(useurhf) then
                  tdhfeda(i,j) = cibcm(i,j)
               else
                  tdhfed(i,j) = cibcm(i,j) + 0.5d0*tdhfed(i,j)
               end if
            end do
         end do
      end if
c------------------------------------------------------------------------      

      call flush(ciout)
      
      end subroutine

      
