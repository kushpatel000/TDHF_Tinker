c     #########################################################
c     ##  This code is written by Kush Patel.                ##
c     ##                                                     ##
c     ##  Prints various matrices and calls the TDHF         ##
c     ##  iterator to update rho according to the            ##
c     ##  equation:                                          ##
c     ##                                                     ##
c     ##      (ih)(d rho/dt) = <[F(rho),rho]> = i L rho      ##
c     ##                                                     ##
c     #########################################################
      
      subroutine pitduhf
      use sizes
      use atoms
      use atomid
      use files
      use iounit
      use orbits
      use piorbs
      use bndstr
      use units
      use tdhfvars
      
      implicit none
      
c     pitdhf variables
      integer i,j,k,m
      integer info
      complex*16 tau(norbit),work(8*norbit),fockcopy(norbit,norbit)
      real*8 d(norbit),e(norbit-1), p
c     variables for reconstructing fock
      integer iorb,jorb
      integer iatn,jatn
      complex*16 s1a,s2a,s1b,s2b
      real*8 hcii, hcij
      real*8 brij,erij
      real*8 bebond,eebond
      real*8 ebb,ebe,abnz,aeth,ble,blb
      real*8 xij,yij,zij,rij,rijsq
      real*8 gii,gij,g11,g11sq,g12,g14
      real*8 ovlap,covlap
      real*8 cionize,iionize,jionize
      real*8 ip(norbit)
c     variables for population analysis
      real*8 Smatrix(norbit,norbit), pop(norbit)
c     variables for iterative output
      integer lext
      character*7 ext
      integer freeunit
      logical exist
c     variables for total energy calculation
      complex*16 energy, paij, pbij
c     Dipole calculation
      real*8 hcom(3)
      real*8 ecom(3)
      real*8 qii
      real*8 dipl(3)
c-------Testing-----------
c     Seeing if original method of calculating energy
c     results in regular fluctuations
      complex*8 xi,xj,xk,xg,xcor
      complex*8 pii,pjj,pij
c-------------------------
c     Electron current 
      real*8 currax, curray, curraz
      real*8 currbx, currby, currbz
      real*8 cc1,cc2,cc3
      
 300  format(2E17.6, 3X)
 301  format(8E16.5)


c     open an output file
      outiter = outiter + 1
      printq = mod(outiter,outfreq).eq.0
      printq = printq.or.(outiter.eq.1)
      lext = 7
      if(printq) then
         call numeral(outiter,ext,lext)
         tdhfout = freeunit()
         inquire(file=filename//'.tdhf.'//ext(1:lext),exist=exist)
         if(exist) then
            open(unit=tdhfout,
     &       file=filename(1:leng)//'.tdhf.'//ext(1:lext),status='old')
            rewind(unit=tdhfout)
         else
            open(unit=tdhfout,
     &       file=filename(1:leng)//'.tdhf.'//ext(1:lext),status='new')
         end if
         
 309     format(A18,E16.5)
         write(tdhfout,*) '## Begin piTDHF iteration ##'
         write(tdhfout,*) ''
         write(tdhfout,309) ' Coupling radius: ', tdhfcr
         write(tdhfout,309) ' Time Step:       ', tdhfdt
         write(tdhfout,*) ''
      end if
      
c----------------------------------------------------------
c     Following code develops the gamma matrix as well
c     as the hc matrix. These are the various integrals
c     and the Hatree Core matrix
c     Taken straight from the piscf subroutine
c     of picalc.f
c----------------------------------------------------------
c     initialize some constants and parameters
      cionize = -11.16d0 / evolt

c     set the bond energies, alpha values and ideal bond length
c     parameter for carbon-carbon pibond type parameters
      ebe  = 129.37d0
      ebb  = 117.58d0
      aeth =   2.309d0
      abnz =   2.142d0
      ble  =   1.338d0
      blb  =   1.397d0
     
c     assign empirical one-center Coulomb integrals, and
c     first or second ionization potential depending on
c     whether the orbital contribures one or two electrons
      do i = 1, norbit
         iorb = iorbit(i)
         tdhfgamma(i,i) = emorb(iorb)
         ip(i) = worb(iorb) + (1.0d0-qorb(iorb))*emorb(iorb)
      end do
      
c     calculate two-center repulsion integrals
c     according to Ohno's semi-empirical formula
      do i = 1, norbit
         iorb = iorbit(i)
         gii = tdhfgamma(i,i)
         do j = i+1, norbit
            jorb = iorbit(j)
            g11 = 0.5d0 * (gii+tdhfgamma(j,j))
            g11sq = 1.0d0 / g11**2
            xij = x(iorb) - x(jorb)
            yij = y(iorb) - y(jorb)
            zij = z(iorb) - z(jorb)
            rijsq = (xij**2 + yij**2 + zij**2) / bohr**2
            g12 = 1.0d0 / sqrt(rijsq + g11sq)
            tdhfgamma(i,j) = g12
            tdhfgamma(j,i) = g12
         end do
      end do

c     the first term in the sum to find alpha is the first
c     or second ionization potential, then the two-center
c     repulsion integrals are added
      do i = 1, norbit
         hcii = ip(i)
         do j = 1, norbit
            if(i.ne.j) then
               jorb = iorbit(j)
               hcii = hcii - qorb(jorb)*tdhfgamma(i,j)
            end if
         end do
         tdhfhc(i,i) = hcii
      end do

c     get two-center repulsion integrals via Ohno's formula
c      do k = 1, nbpi
c         i = ibpi(2,k)
c         j = ibpi(3,k)
       do i = 1, norbit-1
       do j = i+1, norbit
         iorb = iorbit(i)
         jorb = iorbit(j)
         iatn = atomic(iorb)
         jatn = atomic(jorb)
         xij = x(iorb) - x(jorb)
         yij = y(iorb) - y(jorb)
         zij = z(iorb) - z(jorb)
         rij = sqrt(xij**2 + yij**2 + zij**2)
         rijsq = rij**2 / bohr**2
         g11 = 0.5d0 * (tdhfgamma(i,i) + tdhfgamma(j,j))
         g11sq = 1.0d0 / g11**2
         g12 = tdhfgamma(i,j)

c     compute the bond energy using a Morse potential
         erij = aeth * (ble-rij)
         brij = abnz * (blb-rij)
         eebond = (2.0d0*exp(erij)-exp(2.0d0*erij)) * ebe / hartree
         bebond = (2.0d0*exp(brij)-exp(2.0d0*brij)) * ebb / hartree

c     compute the carbon-carbon resonance integral using
c     the Whitehead and Lo formula
         g14 = 1.0d0 / sqrt(4.0d0*rijsq+g11sq)
         hcij = 1.5d0*(bebond-eebond) - 0.375d0*g11
     &        + (5.0d0/12.0d0)*g12 - g14/24.0d0

c     if either atom is non-carbon, then factor the resonance
c     integral by overlap ratio and ionization potential ratio
         if(iatn.ne.6 .or. jatn.ne.6) then
            call overlap(iatn,jatn,rij,ovlap)
            call overlap(6,6,rij,covlap)
            hcij = hcij * (ovlap/covlap)
            iionize = ip(i)
            if(qorb(iorb) .ne. 1.0d0) then
               if(iatn .eq.  7) iionize = 0.595d0 * iionize
               if(iatn .eq.  8) iionize = 0.525d0 * iionize
               if(iatn .eq. 16) iionize = 0.890d0 * iionize
            end if
            jionize = ip(j)
            if(qorb(jorb) .ne. 1.0d0) then
               if(jatn .eq.  7) jionize = 0.595d0 * jionize
               if(jatn .eq.  8) jionize = 0.525d0 * jionize
               if(jatn .eq. 16) jionize = 0.890d0 * jionize
            end if
            hcij = hcij * (iionize+jionize)/(2.0d0*cionize)
         end if
c-----Test Code 11.04.16--------------------------------------
c     This code is the pitdhf copy of the similar thing in fullpi.f.
c     Applies Yukowa scaling to hcij terms of nonbonded atoms
c-------------------------------------------------------------
         bonded = .false.
         do k=1, nbpi
            iorb = ibpi(2,k)
            jorb = ibpi(3,k)
            bonded = (iorb.eq.i .and. jorb.eq.j) .or. bonded
         end do
         if(.not.bonded) then
            hcij = hcij*(exp(-1.0d0*rij/tdhfcr))/rij            
         end if
c-------------------------------------------------------------

c     set symmetric elements to the same value
         tdhfhc(i,j) = hcij
         tdhfhc(j,i) = hcij
      end do
      end do
            
cc     reconstruct the Fock matrix
c      do i=1, norbit
c         do j=1, norbit
c            gij = tdhfgamma(i,j)
c            s1 = (0.0d0,0.0d0)
c            s2 = -1.0d0*tdhfed(i,j)*gij
c
c            if(i.eq.j) then
c               do k=1, norbit
cc                 s1 = s1 + 2.0d0*tdhfgamma(i,k)*tdhfed(k,k) original
c                  s1 = s1 + 1.0d0*tdhfgamma(i,k)*tdhfed(k,k)
c               end do
c            end if
c            tdhffock(i,j) = tdhfhc(i,j) + s1 + s2
c         end do
c      end do     
      
c     reconstruct the Fock Matrix
      do i=1, norbit
         do j=1, norbit
            s1a = (0.0d0, 0.0d0) ! Total Density matrix element
            s1b = (0.0d0, 0.0d0) ! rho(k,k)*J(i,k)
            s2a = -1.0d0*tdhfgamma(i,j)*tdhfeda(i,j)
            s2b = -1.0d0*tdhfgamma(i,j)*tdhfedb(i,j)
     
            if(i.eq.j) then
               do k=1, norbit
                  s1a = tdhfeda(k,k) + tdhfedb(k,k)
                  s1b = s1b + s1a*tdhfgamma(i,k)
               end do
            end if

            tdhffocka(i,j) = tdhfhc(i,j) + s2a + s1b
            tdhffockb(i,j) = tdhfhc(i,j) + s2b + s1b   
         end do
      end do

c     
c     diagonalize the fock matrices to obtain eigenvalues
c     for the superoperator
c
      do m=1,2
c
c     print density matrices
c
         if(printq) then
            if(m.eq.1) then
               write(tdhfout,*) 'TDHF Electron Density Matrix (alpha)'
            else
               write(tdhfout,*) 'TDHF Electron Density Matrix (beta)'
            end if
            do i=1, norbit
               do j=1, norbit
                  if(m.eq.1) then
                     write(tdhfout,300,advance='no') tdhfeda(i,j)
                  else
                     write(tdhfout,300,advance='no') tdhfedb(i,j)
                  end if
               end do
               write(tdhfout,*) ''
            end do
            write(tdhfout,*) ''

c     Print the diagonal elements
 307        format(2E20.9,3X)
            if(m.eq.1) then
               write(tdhfout,*) 'TDHF Density Matrix Diagonal (alpha)'
            else
               write(tdhfout,*) 'TDHF Density Matrix Diagonal (beta)'
            end if
         end if
         s1a = (0.0d0,0.0d0)
         do j=1, norbit
            if(m.eq.1) then
               s1a = s1a + tdhfeda(j,j)
               if(printq) write(tdhfout,307) tdhfeda(j,j)
            else
               s1a = s1a + tdhfedb(j,j)
               if(printq) write(tdhfout,307) tdhfedb(j,j)
            end if  
         end do

         if(printq) then
 302        format(A12,2F17.6)
            write(tdhfout,*) ''
            write(tdhfout,302) 'Trace: ', s1a
            write(tdhfout,*) ''
            
            s1a = (0.0d0,0.0d0)
            do i=1, norbit
               s1a = s1a + tdhfeda(i,i) + tdhfedb(i,i)
            end do
            
            if(m.eq.2) then
 303           format(A19,2F17.6)
               write(tdhfout,303) 'Trace (full): ', s1a
               write(tdhfout,*) ''
            end if
         end if
c
c     copy fock matrix to tdhfcopy
c
         if(m.eq.1) then
            if(printq) write(tdhfout,*) 'Fock matrix (alpha)'
               do i=1, norbit
                  do j=1, norbit
                     fockcopy(i,j) = tdhffocka(i,j)                  
                     if(printq) then
                        write(tdhfout,300,advance='no') tdhffocka(i,j)
                     end if
                  end do
                  if(printq) write(tdhfout,*) ''
               end do
            else
               if(printq) write(tdhfout,*) 'Fock matrix (beta)'
               do i=1, norbit
                  do j=1, norbit
                     fockcopy(i,j) = tdhffockb(i,j)                  
                     if(printq) then
                        write(tdhfout,300,advance='no') tdhffockb(i,j)
                     end if
                  end do
                  if(printq) write(tdhfout,*) ''
               end do
            end if
            if(printq) write(tdhfout,*) ''
c
c     occasionally recalculate the eigenvalues
c
            if(mod(outiter,100).eq.1 .or. res1) then
               res1 = .false.
               if(tdhfdebug) then
                  write(debugout,*) 'Recalculating Super Eigenvalues'
               end if
c
c     tridiagonalize
c     
            call zhetrd('U',norbit,fockcopy,norbit,d,e,tau,
     &           work,8*norbit,info)
            if(tdhfdebug) then
               write(debugout,*) 'Info: ', info
               write(debugout,*) ''
               write(debugout,*) 'Post Tridiagonal Fock'
               do i=1,norbit
                  do j=1, norbit
                     if(m.eq.1) then
                        write(debugout,300,advance='no') tdhffocka(i,j)
                     else
                        write(debugout,300,advance='no') tdhffockb(i,j)
                     end if
                  end do
                  write(debugout,*)''
               end do
               write(debugout,*)''
               
               write(debugout,*) 'tri-diagonal matrix'
               write(debugout,*) 'diagonal     superdiagonal'
               do i=1, norbit
                  write(debugout,301,advance='no') d(i)
                  if(i.ne.norbit) then
                     write(debugout,301,advance='no') e(i)
                  end if
                  write(debugout,*) ''
               end do
               write(debugout,*) ''
            end if
            
c     get the eigenvalues
            if(m.eq.1) then
               superevmax = 0.0d0
               superevmin = 0.0d0
            end if
            
            call zsteqr('N',norbit,d,e,'null',norbit,'null',info)
            if(printq) then 
               if (m.eq.1) then
                  write(tdhfout,*) 'Fock (alpha) Eigenvalues: '
               else
                  write(tdhfout,*) 'Fock (beta) Eigenvalues: '
               end if
               write(tdhfout,301) (d(i),i=1,norbit)
               write(tdhfout,*)''
            end if
            superevmax = max(superevmax,d(norbit)-d(1))
            superevmin = min(superevmin,d(1)-d(norbit))
            if(printq) then
               write(tdhfout,*) 'Supereigenvalues: {', superevmin,', ',
     &              superevmax, '}'
               write(tdhfout,*) ''
            end if
            superevmax = 1.5d0*superevmax
            superevmin = 1.5d0*superevmin
            
         end if

c
c     calculate total electronic energy
c
         if(printq) then
            if(m.eq.2) then
               energy = (0.0d0,0.0d0)
               s1a = (0.0d0,0.0d0)
               s1b = (0.0d0,0.0d0)
               do i=1, norbit
                  do j=1, norbit
                     paij = tdhfeda(i,j)
                     pbij = tdhfedb(i,j)
                     s1a = s1a + 0.5d0*paij*(tdhfhc(j,i)+tdhffocka(j,i))
                     s1b = s1b + 0.5d0*pbij*(tdhfhc(j,i)+tdhffockb(j,i))
                  end do
               end do
               energy = s1a + s1b
 305           format(A25,2E16.5)
               write(tdhfout,305) 'Total alpha Energy: ', s1a
               write(tdhfout,305) 'Total beta Energy: ' , s1b
               write(tdhfout,305) 'Total Electronic Energy: ', energy
               write(tdhfout,*) ''
            end if
c
c     calculate the effective center of mass of electron and hole
c        
            if(m.eq.2) then
               do i=1, 3
                  hcom(i) = 0.0d0
                  ecom(i) = 0.0d0
                  dipl(i) = 0.0d0
               end do
               do i=1, norbit
c                 qii = density - charge = net population
                  qii = realpart(tdhfeda(i,i)+tdhfedb(i,i)) - qorb(i)
                  if(qii.gt.0.0d0) then
                     ecom(1) = ecom(1) + qii*x(i) !/norbit
                     ecom(2) = ecom(2) + qii*y(i) !/norbit
                     ecom(3) = ecom(3) + qii*z(i) !/norbit
                  else
                     hcom(1) = hcom(1) + abs(qii)*x(i) !/norbit
                     hcom(2) = hcom(2) + abs(qii)*y(i) !/norbit
                     hcom(3) = hcom(3) + abs(qii)*z(i) !/norbit           
                  end if
                  dipl(1) = dipl(1) + qii*x(i)
                  dipl(2) = dipl(2) + qii*y(i)
                  dipl(3) = dipl(3) + qii*z(i)
               end do
 304           format(A12,3E15.5)
               write(tdhfout,304) 'Hole CoM: ', (hcom(i),i=1,3)
               write(tdhfout,304) 'Elec CoM: ', (ecom(i),i=1,3)
               write(tdhfout,304) 'Dipole:   ', (dipl(i),i=1,3)
               write(tdhfout,*)
            end if
         end if

         if(m.eq.2 .and. probcurr .and. printq) then
            write(tdhfout,*) 'Probability Current'
            write(tdhfout,*) ' qx   qy    qz    ax     ay   az bx by bz'
            do i=1, norbit
               iorb = iorbit(i)
               currax = 0.0d0
               curray = 0.0d0
               curraz = 0.0d0
               currbx = 0.0d0
               currby = 0.0d0
               currbz = 0.0d0
               do j=1, norbit
                  jorb = iorbit(j)
                  cc1 = 2.0d0*tdhfhc(j,i)*aimag(tdhfeda(j,i))
                  cc2 = 2.0d0*tdhfhc(j,i)*aimag(tdhfedb(j,i))
                  currax = currax + cc1* (x(jorb)-x(iorb))
                  curray = curray + cc1* (y(jorb)-y(iorb))
                  curraz = curraz + cc1* (z(jorb)-z(iorb))
                  currbx = currbx + cc2* (x(jorb)-x(iorb))
                  currby = currby + cc2* (y(jorb)-y(iorb))
                  currbz = currbz + cc2* (z(jorb)-z(iorb))
               end do
 306           format(9E16.5)
               write(tdhfout,306) x(iorb),y(iorb),z(iorb),currax,curray,
     &              curraz,currbx,currby,currbz
            end do
         end if
c     
c     call tdhfiterate and update rho
c
c         if(m.eq.1) then
c            call tdhfiterate(tdhfeda,tdhffocka,superevmax,
c     &           superevmin,tdhfdt,norbit,15)
c         else
c            call tdhfiterate(tdhfedb,tdhffockb,superevmax,
c     &           superevmin,tdhfdt,norbit,15)
c         end if         
      end do

c     call tdhfiterate and update rho
      call tdhfiterate(tdhfeda,tdhffocka,superevmax,
     &     superevmin,tdhfdt,norbit,9)
      call tdhfiterate(tdhfedb,tdhffockb,superevmax,
     &     superevmin,tdhfdt,norbit,9)


c-----Test Code 10 Dec 17-----------------------------------
c     Impose zero for imaginary parts on diagonal
      do i=1,norbit
         tdhfeda(i,i) = realpart(tdhfeda(i,i))
         tdhfedb(i,i) = realpart(tdhfedb(i,i))
      end do      
c-----------------------------------------------------------      
      
cc     Calculate the Overlap matrix to get the population
c      do i=1,norbit
c         do j=1,norbit 
c            Smatrix(i,j) = 0.0d0
c         end do
c      end do
c      do k=1,nbpi
c         i = ibpi(2,k)
c         j = ibpi(3,k)
c         iorb = iorbit(i)
c         jorb = iorbit(j)
c         iatn = atomic(iorb)
c         jatn = atomic(jorb)
c         xij = x(iorb) - x(jorb)
c         yij = y(iorb) - y(jorb)
c         zij = z(iorb) - z(jorb)
c         rij = sqrt(xij**2 + yij**2 + zij**2)
c         call overlap(iatn,jatn,rij,ovlap)
c         Smatrix(i,j) = ovlap
c         call overlap(jatn,iatn,rij,ovlap)
c         Smatrix(j,i) = ovlap
c      end do
c      do i=1,norbit
c         Smatrix(i,i) = 1.0d0
c      end do
c            
c      write(tdhfout,*) 'Overlap Matrix'
c      do i=1,norbit
c         do j=1,norbit
c            write(tdhfout,301,advance='no') Smatrix(i,j)
c         end do
c         write(tdhfout,*) ''
c      end do
c      write(tdhfout,*) ''
c
c      call popanal('M',norbit,tdhfed,Smatrix,pop)
c      write(tdhfout,*) 'Orbital Populations'
c      do i=1, norbit
c         write(tdhfout,301) pop(i)
c      end do
c      write(tdhfout,*)
      
c     
c     update the nonplanar pi bond orders (pnpl)
c     this is the same code that's at the end of picalc
c     
 401  format(5i5,2f12.2)
 402  format(3A5,2A15)
      if(printq) write(tdhfout,402) 'Bond','i','j','old','new'
      do k=1, nbpi
         i = ibpi(2,k)
         j = ibpi(3,k)
         p = pnpl(k)
         
c        take just the real part of electron density
c         pnpl(k) = realpart( tdhfed(i,j) )
c        take the real part of the electron density
c         pnpl(k) = zabs( tdhfeda(i,j) + tdhfedb(i,j))
         pnpl(k) = realpart( tdhfeda(i,j) + tdhfedb(i,j) )
         pbpl(k) = pnpl(k) * tdhfhc(i,j)/(-0.0757d0)
         
         i = ibnd(1,ibpi(1,k))
         j = ibnd(2,ibpi(1,k))
         if(printq)
     &        write(tdhfout,401) k,ibpi(2,k),ibpi(3,k),i,j,p,pnpl(k)
      end do
      if(printq) write(tdhfout,*)''
         
c-------------------------
      if(printq) then
         xi = 0.0d0
         xj = 0.0d0
         xk = 0.0d0
         xg = 0.0d0
         xcor = 0.0d0
         
         write(tdhfout,*) 'Original Type Energy Calculation'
         do i=1, norbit
            do j=1, norbit
               pii = tdhfeda(i,i) + tdhfedb(i,i)
               pjj = tdhfeda(j,j) + tdhfedb(j,j)
               pij = tdhfeda(i,j) + tdhfedb(i,j)
               xi = xi + pij*tdhfhc(i,j)
               xj = xj + 0.50d0*pij*tdhfgamma(i,j)
               xk = xk - 0.25d0*pij*tdhfgamma(i,j)
               xcor = xcor - 0.5d0*(pij*pij - pii*pjj)*tdhfgamma(i,j)
            end do
         end do
         do i=1, norbit-1
            do j=i+1, norbit
               xg = xg + tdhfgamma(i,j)
            end do
         end do
         write(tdhfout,*) 'OEnergy:  ', xi+xj+xk+xg+xcor
         write(tdhfout,*) 'Core:     ', xi
         write(tdhfout,*) 'Coulomb:  ', xj
         write(tdhfout,*) 'Exchange: ', xk
         write(tdhfout,*) 'Nuclear:  ', xg
         write(tdhfout,*) 'Xcorrect: ', xcor
c-------------------------

         write(tdhfout,*)
         write(tdhfout,*)'## End of piTDHF iteration ##'
         write(tdhfout,*)''
         write(tdhfout,*)''
      
c     Close the output file
         flush(tdhfout)
         close(unit=tdhfout)

c     Write restart information
         call writerestart()

      end if
         
      end subroutine
