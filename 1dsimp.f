      program lahyc
c*******************************************************            
c                                                      *
c  This is a Lagrangian 1D hydrodynamics code.         *
c  adds particles                                      *
c  this is the main driver of the integration.         *
c                                                      *
c*******************************************************
c
      implicit double precision (a-h,o-z)
c
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
c
      common /bstuf/ rb, dumrb, f1rb, f2rb
      common /carac/ deltam(idim), abar(idim)
      common /cases/ t9nse, rhoswe, rhonue, rhonux
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim)
      common /celle/ x(0:idim),v(0:idim),f(0:idim)
      common /cgas / gamma
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne
      common /ener2/ tkin, tterm
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /epcap/ betafac, c2cu, c3cu
      common /numb / ncell, ncell1
      common /propt/ dtime,tmax
      common /shock/ cq,cl
      common /tempe/ temp(idim)
      common /therm/ xmu(idim)
      common /timej / time, dt
      common /typef/ iextf, ieos
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu
c
c--read initial condition
c
      call readini
c
c--define code's units
c
      call unit
c
c--set all quantities needed for run
c
      call preset
c
c--main loop
c
      nstep = 0
 42   continue
c
      nstep = nstep + 1
c
c--step runs the runge-kutta operator
c
      call step
c
c--produce output if desired
c
      lu=61
      call printout(lu)
      print *, 'savetime=',time
      if (time.gt.1.0) dtime = .1
      if (time.gt.10.0) dtime = 1.
      if (time.gt.100.0) dtime = 10.
      if (time.gt.1000.0) dtime = 100.
      if (time.gt.10000.0) dtime = 1000.
      if (time.gt.100000.0) dtime = 10000.
      if (time.gt.1000000.0) dtime = 100000.
      if (time.gt.10000000.0) dtime = 1000000.
c
c--check if simulation is finished
c
      print *,time,tmax
      if(time.gt.tmax) go to 90
c
      go to 42
c
   90 call printout(lu)
      stop
      end
      subroutine hydro(time,ncell,x,v,
     1           u,rho,ye,f,du,dye,q)
c
c****************************************************************
c                                                               *
c  this subroutine advances the system of hydro equations by    *
c  one time step.                                               *
c                                                               *
c  Description of the variables :                               *
c  ------------------------------                               *
c                                                               *
c  ncell        :  number of cells                              *
c  ncell1       :  number of cell edges                         *
c  rho          :  density (cell centered)                      *
c  pr           :  pressure (cell centered)                     *
c  u            :  specific internal energy (cell centered)     *
c  q            :  artificial viscous stress (cell centered)    *
c  deltam       :  mass of the cell (cell centered)             *
c  prold,qold,rhold  : pressure artificial viscous stress and   *
c                     density at previous time step             *
c  x            :  cell boundaries (edge centered)              *
c  v            :  velocity (edge centered)                     *
c  vold         :  velocity at previous time step               *
c  gamma        :  ratio of specific heat                       *
c  cq, cl       :  quadratic and linear artificial viscous      *
c                  stress coefficients                          *
c  tkin, uint   :  kinetic and specific internal energy         *
c  time         :  current time                                 *
c  dt, dtold    :  current and previous time step               *
c  tmax         :  maximum time for the simulation              *
c                                                               *
c****************************************************************
c
c  cell
c center          1         2         3
c------------|---------|---------|---------|---
c  cell      0         1         2         3
c  edge
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
      parameter (tiny=-1e-5)
c
      dimension x(0:idim), v(0:idim), f(0:idim)
      dimension u(idim), rho(idim), ye(idim)
      dimension dye(idim), du(idim), q(idim)
c
      logical ebetaeq, pbetaeq
      common /beta/ ebetaeq(idim), pbetaeq(idim)
      common /bstuf/ rb, dumrb, f1rb, f2rb
      common /carac/ deltam(idim), abar(idim)
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne
      common /ener1/ dq(idim), dunu(idim)
      common /ener2/ tkin, tterm
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /epcap/ betafac, c2cu, c3cu
      common /shock/ cq,cl
      common /tempe/ temp(idim)
      common /therm/ xmu(idim)
      common /typef/ iextf, ieos
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu
c      common /nuout/ rlumnue, rlumnueb, rlumnux,
c     1               enue, enueb, enux, e2nue, e2nueb, e2nux
c
c--compute density
c ---------------------------------------------------------
c
      call density(ncell,x,rho)
c
c
c--compute thermodynamical properties
c  ----------------------------------
c  (this allows the use of simple eos)
c
      print *, 'time=',time,ieos
      if(ieos.eq.1)call eospg(ncell,rho,u)
      if(ieos.eq.2)call eospgr(ncell,rho,u)
c
c--do gravity
c  --------------------------------------------------------
c
      call gravity(ncell,deltam,x,f)
c
c--compute q values
c
      call artvis(ncell,x,rho,v,q)
c
c--compute forces on the particles
c
      call forces(ncell,x,f,q,v,rho)
c
      call energ(ncell,x,v,dye,du,rho,time)
c
   99 return
      end
c
      subroutine artvis(ncell,x,rho,v,q)
c***********************************************************
c                                                          *
c  This subroutine updates the q-value for artificial      *
c  viscosity.                                              *
c                                                          *
c***********************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
c
      dimension x(0:idim),v(0:idim)
      dimension rho(idim),q(idim)
c
      common /ener1/ dq(idim), dunu(idim)
      common /shock/ cq,cl
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /carac/ deltam(idim), abar(idim)
c
      data pi4/12.56637d0/
c
c--update q value
c
      do kp05=1,ncell
         k=kp05 - 1
         k1=kp05
         q(kp05)=0.d0
         akp1=pi4*x(k1)*x(k1)
         ak=pi4*x(k)*x(k)
         akp05=0.5d0*(ak+akp1)
         gradv=v(k1)*(3.d0*akp05-akp1) - v(k)*(3.d0*akp05-ak)
         if(gradv.lt.0.d0)then
c
c--quadratic term
c
            dv=v(k1) - v(k)
            alpha=0.5d0*cq*cq*rho(kp05)*dabs(dv)/akp05
            q(kp05)=-alpha*gradv
c
c--linear term
c
            cs=vsound(kp05)
            alphal=0.5d0*cl*rho(kp05)*cs/akp05
            q(kp05)=q(kp05) - alphal*gradv
         end if
         dq(kp05)=-q(kp05)*gradv/deltam(kp05)
      enddo
c
      return
      end
c
      subroutine density(ncell,x,rho)
c****************************************************************
c                                                               *
c  This subroutine calculates the density using the continuity  *
c  equation.                                                    *
c                                                               *
c****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
c
      dimension x(0:idim),rho(idim)
      common /carac/ deltam(idim), abar(idim)
c
      data pi4/12.566371/
c
c--update density
c
      do kp05=1,ncell
         k1=kp05
         k=kp05 - 1
         rho(kp05)=3.d0*deltam(kp05)/
     1              (pi4*(x(k1)*x(k1)*x(k1)-x(k)*x(k)*x(k)))
      enddo
      return
      end
c
      subroutine energ(ncell,x,v,dye,du,rho,time)
c************************************************************
c                                                           *
c  this routine computes the change in internal energy      *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
      parameter (avokb=6.02e23*1.381e-16)
c
      dimension x(0:idim),v(0:idim)
      dimension du(idim), dye(idim), rho(idim)
c
      double precision lc
      logical ebetaeq, pbetaeq
      common /beta/ ebetaeq(idim), pbetaeq(idim)
      common /tempe/ temp(idim)
      common /typef/ iextf, ieos
      common /cpots/ xmue(idim), xmuhat(idim)
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /ener1/ dq(idim), dunu(idim)
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /carac/ deltam(idim), abar(idim)
      double precision dupp
      common /lum/ lc
c
      data pi4/12.566371/
c
c--compute change in specific internal energy
c
c--entropy conversion factor
      sfac=avokb*utemp/uergg      
c
      do kp05=1,ncell
         k=kp05-1
         k1=kp05
         akp1=pi4*x(k1)*x(k1)
         akp=pi4*x(k)*x(k)
         pdv=pr(kp05)*(akp1*v(k1)-akp*v(k))/deltam(kp05)
c
         du(kp05)=-pdv+0.5*dq(kp05)
         if (kp05.eq.1) then
c            du(kp05)=du(kp05)+3.4d44*(1.d0+time/5.d2)**(-1.0)
c     $           /2.d49/deltam(kp05)
c            du(kp05)=du(kp05)+9.6d42*(1.d0+time/9.1d2)**(-1.0)
c     $           /2.d49/deltam(kp05)
         end if
c         if (kp05.lt.75) then
c            du(kp05)=du(kp05)+2.d10/1.d16*((time+1.0)/8640)**(-1.3)
c         end if
         if (du(kp05).gt.1.d15) then
            print *, 'large du',kp05,pdv,dq(kp05)
            stop
         end if
      enddo
c
      return
      end
c
      subroutine eospg(ncell,rho,u)     
cc************************************************************
c                                                           *
c  This subroutine computes the pressure and sound speed    *
c  for all particles on a list assuming a perfect gas       *
c  equation of state.                                       *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
c
      dimension rho(idim), u(idim)
c
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /prev/ xold(0:idim),vold(0:idim),rhold(idim),
     1              prold(idim), tempold(idim), yeold(idim)
      common /carac/ deltam(idim), abar(idim)
      common /cgas / gamma
c
c--initialize quantities
c
      vsmax=0.
      gama1=gamma - 1.
c
c--isothermal equation of state
c
      if(gamma.eq.1.)then
         do k=1,ncell
            prold(k) = pr(k)
            pr(k)=u(k)*rho(k)
            vsound(k)=dsqrt(pr(k)/rho(k))
            vsmax=dmax1(vsmax,vsound(k))
         enddo
         return
      end if
c
c--perfect gas equation of state
c
      do k=1,ncell
         prold(k) = pr(k)
         pr(k)=u(k)*gama1*rho(k)
c         if (k.eq.1) print *, 'energy total',u(k),gama1
         vsound(k)=dsqrt(gamma*pr(k)/rho(k))
         vsmax=dmax1(vsound(k),vsmax)
      enddo
      return
      end
c
      subroutine eospgr(ncell,rho,u)
c************************************************************
c                                                           *
c  This subroutine computes the pressure, and sound speed   *
c  according to an equation of state that includes gas and  *
c  radiation pressure.                                      *
c  This part of the code has not been debugged and most     *
c  likely won't run.                                        *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
c
      dimension rho(idim), u(idim)
c
      double precision umass
      common /prev/ xold(0:idim),vold(0:idim),rhold(idim),
     1              prold(idim), tempold(idim), yeold(idim)
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /carac/ deltam(idim), abar(idim)
      common /tempe/ temp(idim)
      common /therm/ xmu(idim)
      common /cgas / gamma
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne
      common /units/ umass, udist, udens, utime, uergg, uergcc
c
      double precision rhok, uk, adrho, t9, rdmu, pgas, prad
c
      vsmax=0.
      gama1=gamma - 1.
c
      tempmx=-1e20
      tempmn=1e20
      do k=1,ncell
         rhok=dble(rho(k))
         uk=dble(u(k))
         adrho=dble(arad)/rhok
         t9=dble(temp(k))
         rdmu=dble(bigr/xmu(k))
         call rootemp1(t9,uk,rdmu,adrho,iflag)
         if(iflag.eq.1) then
            write(*,*)'temperature did not converge!',k
            write(*,*)t9,ui, rhok
         end if
         temp(k)=t9
         tempmx=dmax1(temp(k),tempmx)
         tempmn=dmin1(temp(k),tempmn)
c
c--compute the various pressures
c
         prad=dble(arad)*t9*t9*t9*t9*0.3333333333333d0
         pgas=rhok*t9*rdmu
         prold(k) = pr(k)
         pr(k)=pgas+prad
c
c--compute thermal energy and sound speed
c
         vsound(k)=dsqrt(1.66666666*real(pgas)/rho(k))
         vsmax=dmax1(vsound(k),vsmax)
      enddo
      print*,'eospgr:max,min temp(K)',tempmx*1e9,tempmn*1e9
c
      return
      end
      subroutine forces(ncell,x,f,q,v,rho)
c****************************************************************
c                                                               *
c  This subroutine computes the force on the cells that need to *
c  have their forces evaluated.  We also do  neutrino diffusion *
c  here.                                                        *
c                                                               *
c****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
c
      dimension x(0:idim),f(0:idim),v(0:idim)
      dimension q(idim),rho(idim)
c
      common /eosnu/ prnu(idim1)
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /carac/ deltam(idim), abar(idim)
      common /damping/ damp, dcell
c
      data pi4/12.56637d0/
c
c--acceleration by pressure
c
      prnu(1)=0.d0
      prnu(ncell+1)=0.d0
      do k=1,ncell
         km05=k
         kp05=k+1 
         ak=pi4*x(k)*x(k)
         xk3=(.5d0*(x(k-1)+x(k)))**2.5
         xkp3=(.5d0*(x(k)+x(k+1)))**2.5
	 akp1=pi4*x(k+1)*x(k+1)
	 akm1=pi4*x(k-1)*x(k-1)
	 akp05=0.5d0*(akp1+ak)
	 akm05=0.5d0*(ak+akm1)
c         deltamk=0.5d0*(deltam(km05)+deltam(kp05))
c
c--pressure gradients
c
         pressp = pr(kp05)
         pressm = pr(km05)
         gradp=ak*(pressp - pressm)
c
c--artificial viscosity pressure
c
         gradq=0.5d0*q(kp05)*(3.d0*akp05-ak) - 
     1        0.5d0*q(km05)*(3.d0*akm05-ak)
c
c         write (99,199) k,x(k),f(k),
c     $        gradp/deltam(km05),gradq/deltam(km05)
c         if (km05.eq.1) print *, 'pressure grad',km05,gradp,gradq,
c     $        pr(km05),pr(kp05)
         if (gradp.ne.gradp) stop
         f(k)=f(k)-(gradp+gradq)/deltam(km05)
c
c--damping
c
         if (k.lt.int(dcell)) then
            f(k)=f(k)-damp*v(k)
         end if
      enddo
c 199  format(I4,4(1x,1pe12.4))
c
      return
      end
c 
      subroutine gravity(ncell,deltam,x,f)
c****************************************************************
c                                                               *
c  This subroutine adds the gravitational force on the          *
c  particles.  This will have the option of a neutron star      *
c  core or it can allow the particles to make up the core.      *
c                                                               *
c****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
c
      dimension x(0:idim),f(0:idim)
      dimension gpot(idim),deltam(idim)
      dimension xmi(0:idim)
c
      common /core / dcore, xmcore
      common /rshift/ gshift(idim)
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne
c
c--compute internal mass for all edges
c--and calculate forces (actually accelerations)
c
      if(dcore.gt.1.d0) then
         xmi(0)=xmcore
         r2=x(0)**2
         f(0)=-gg*xmi(0)/r2
      else
         xmi(0)=0
         f(0)=0
      end if
      do k=1,ncell
         xmi(k)=xmi(k-1)+deltam(k)
         r2=x(k)**2
         f(k)=-gg*xmi(k)/r2
      enddo
c
c--calculate gravitational potential     
c
      const=gg/(clight*clight)
      r=x(ncell)
      gpot(ncell)=-const*xmi(ncell)/r
      rold=r
      do k=ncell-1,1,-1
         r=x(k)
         gpot(k)=gpot(k+1)-(rold-r)*const*xmi(k)/(r*r)
         rold=r
      enddo
c
c--calculate gravitational redshift (w.r.t. r=infinity)
c
      do k=1,ncell
         gshift(k)=1.d0/dsqrt(1.d0-2.0d0*gpot(k))
c         gshift(k)=1.d0
      enddo
c
      return
      end
c
      subroutine mmw
c**************************************************************
c                                                             *
c  subroutine sets the mean molecular weight of the gas       *
c  assuming complete ionization.                              *
c                                                             *
c**************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
c
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim)
      common /numb/ ncell, ncell1
      common /therm/ xmu(idim)
      common /carac/ deltam(idim), abar(idim)
c
      do i=1,ncell
         xmu(i)=abar(i)/(abar(i)*ye(i)+1.)
      enddo
      return
      end
c
      subroutine preset
c****************************************************************
c                                                               *
c  This subroutine sets up all quantities needed before         *
c  starting a simulation.                                       *
c                                                               *
c****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
c
c
      common /carac/ deltam(idim), abar(idim)
      common /cases/ t9nse, rhoswe, rhonue, rhonux
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim)
      common /celle/ x(0:idim),v(0:idim),f(0:idim)
      common /etnus/ etanue(idim),etanueb(idim),etanux(idim)
      common /numb/ ncell,ncell1
      common /therm/ xmu(idim)
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /uswest/ usltemp, uslrho, uslu, uslp, u2sluncell1
      common /nsat/ satc,xtime
      common /heat/ iheat
c
c
c
c
c
c--set up mean molecular weight according to matter types
c
      call mmw
c
c--set the critical densities for Swesty's eos, nu trapping
c
      satc=0.d0
      xtime=1000000.
c
c--remove heated cells
      iheat=7
c
      return
      end
c
      subroutine rootemp1(temp,utr,rdmu,adrho,iflag)
c***************************************************************
c                                                              *
c     subroutine computes temperature using a Newton-Raphson   *
c     procedure found in the numerical recipes, p.254          *
c                                                              *
c***************************************************************
      implicit double precision (a-h,o-z)
c
      data itmax/80/ , tol/1.d-2/ 
c
c--find root allowing a maximum of itmax iteration
c
      iflag=0
      do i=1,itmax
         at93=adrho*temp*temp*temp
         df=4.d0*at93 + 1.5d0*rdmu
         f=at93*temp + 1.5d0*rdmu*temp - utr
         dtemp=f/df
         temp=temp - dtemp
         if(abs(dtemp/temp).lt.tol) return
      enddo
      iflag=1
c
c--iteration went out of bound or did not converge
c
c      print *,'rootemp: iteration did not converge'
      write(*,*)temp,dtemp
c
      return
      end
c
      subroutine readini
c*************************************************************
c
c     reads initial conditions
c
c************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
      parameter (iqn=17)
      real ycc,yccave
c
      common /cc   / ycc(idim,iqn), yccave(iqn)
      common /ener1/ dq(idim), dunu(idim)
      common /numb/ ncell, ncell1
      common /celle/ x(0:idim),v(0:idim),f(0:idim)
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim)
      common /eosq/ pr(idim1), vsound(idim), u2(idim), vsmax
      common /carac/ deltam(idim), abar(idim)
      common /core / dcore, xmcore
      common /typef/ iextf, ieos
      common /bstuf/ rb, dumrb, f1rb, f2rb
      common /shock/ cq,cl
      common /tempe/ temp(idim)
      common /cgas/ gamma
      common /timej/ time, dt
      common /propt/ dtime,tmax
      common /ener2/ tkin, tterm
      common /outp/ rout, p1out, p2out
      common /nuprint/ nups,nupk,tacr
      common /damping/ damp,dcell
      common /neutm/ iflxlm, icvb
      common /ufactor/ ufact,yefact
      logical ifign
      common /ign  / ifign(idim)
      double precision lc
c
      character*9 filin,filout
      data pi4/12.56637d0/
      gg=13.34
      tacr=1.d2
c
c--read options
c
      open(11,file='inlahyc')
      read(11,10) filin
      read(11,10) filout
   10 format(A)
      read(11,*) idump
      read(11,*) dtime,tmax
      read(11,*) cq,cl
      read(11,*) iextf,ieos,dcore
      read(11,*) ncell,delp,nups,damp,dcell
      read(11,*) iflxlm, icvb, ufact, yefact
      print *,'cq,cl',cq,cl
      print *,'iextf,ieos',iextf,ieos
c
c--open binary file containing initial conditions
c
     
      open(60,file=filin,form='unformatted')
      open(61,file=filout,form='unformatted')
c
c--position pointer in binary file
c
      do i=1,idump-1
         read(60) idummy
      enddo
c     
c--read data
c
      nqn=17
c
      read(60) nc,t,xmcore,rb,lc,
     1   (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),
     2      (u(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),
     3      (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),
     7      (pr(i),i=1,nc),(u2(i),i=1,nc),
     $     ((ycc(i,j),j=1,nqn),i=1,nc)
c
      time = t
c
      print *, nc, t, xmcore, r
      print *, lc
      do k=1,nc
         write(44,*)k,temp(k),u(k),x(k),v(k),lc,xmcore
      end do
      p1out=0.
cpr(ncell)
      p2out=1.39*gg*deltam(ncell)/x(ncell)**4/pi4
      rout=x(ncell)
      print *, 'ncell = ', ncell
      print *, 'xmcore = ', xmcore
      print *, 'rho(1) = ', rho(1)
      print *, x(ncell),pr(ncell)
      print *, 'deltam(1) = ', deltam(1)
      print *, 'time = ',time
      gamma=4.d0/3.d0
      ncell1=ncell+1
c
      return
      end
c
      subroutine printout(lu)
c**************************************************************
c                                                             *
c  This subroutine prints out all the results.                *
c                                                             *
c**************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
      parameter (iqn=17)
c
      real ycc,yccave
c
      common /cc   / ycc(idim,iqn), yccave(iqn)
      common /ener1/ dq(idim), dunu(idim)
      common /numb/ ncell, ncell1
      common /celle/ x(0:idim),v(0:idim),f(0:idim)
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim)
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /core / dcore, xmcore
      common /typef/ iextf, ieos
      common /bstuf/ rb, dumrb, f1rb, f2rb
      common /tempe/ temp(idim)
      common /cgas / gamma
      common /timej / time, dt
      common /carac/ deltam(idim), abar(idim)
      common /ener2/ tkin, tterm
      dimension uint(idim), s(idim)
      double precision lc
      common /lum/ lc
c
      do i=1,ncell
         uint(i)=u(i)
      enddo
      t=time
      nc=ncell
      nqn=17
c
c
c      do i=1,1
c         write (12,*) i, nc,t,gamma,tkin,tterm
c      end do
c
      write(lu)nc,t,xmcore,rb,lc,
     1   (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),
     2      (uint(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),
     3      (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),
     7      (pr(i),i=1,nc),(s(i),i=1,nc),
     $     ((ycc(i,j),j=1,nqn),i=1,nc)
c
      return
c
c--an error as occured while writting
c
 10    print *,'an error has occured while writing'
c

      return
      end
c
      subroutine step
c************************************************************
c                                                           *
c  This subroutine integrates the system of equations using *
c  a Runge-Kutta-Fehlberg integrator of second order.       *
c  Particles are allowed to have individual time-steps.     *
c  All particles are synchronized every dtime at which time *
c  the subroutine is exited.                                *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (idim=4000)
      parameter (idim1=idim+1)
c
      logical ifirst, reset(idim)
c
      logical ebetaeq, pbetaeq
      common /beta/ ebetaeq(idim), pbetaeq(idim)
      common /bstuf/ rb, dumrb, f1rb, f2rb
      common /carac/ deltam(idim), abar(idim)
      common /cellc/ u(idim),rho(idim),ye(idim),q(idim)
      common /celle/ x(0:idim),v(0:idim),f(0:idim)
      common /cgas / gamma
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne
      common /core / dcore, xmcore
      common /ener2/ tkin, tterm
      common /eosq / pr(idim1), vsound(idim), u2(idim), vsmax
      common /epcap/ betafac, c2cu, c3cu
      common /numb / ncell, ncell1
      common /prev / xold(0:idim),vold(0:idim),rhold(idim),
     1              prold(idim), tempold(idim), yeold(idim)
      common /shock/ cq,cl
      common /tempe/ temp(idim)
      common /therm/ xmu(idim)
      common /timej / time, dt
      common /typef/ iextf, ieos
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /dum  / dumx(0:idim),dumv(0:idim),dumu(idim),
     1               dumye(idim)
      common /f1   / f1v(0:idim),f1u(idim),f1ye(idim)
      common /f2   / f2v(0:idim),f2u(idim),f2ye(idim)
c
      common /timei/ istep(idim),t0(idim),steps(idim), 
     1               dum2v(idim)
      common /propt/ dtime,tmax
      common /nsat/ satc,xtime
      common /nuprint/ nups,nupk,tacr
      common /outp/ rout, p1out, p2out
c
      save ifirst
      data ifirst/.true./
      data tiny/3.e-5/
c
c--Compute next dump time and initialise variables.
c
      tol=1.d-3
      dt0=1.0d5
      tnext=time+dtime
c      xlog2=0.30103
c
c--define coefficients for Runge-Kutta integrator
c
      f11=0.5*256./255.
      f21=1./256.
      f22=255./256.
      e1=1./512.
c
c--first time set time step to default
c
      if(ifirst) then
         ifirst=.false.
         do i=1,ncell
            steps(i)=dt0/1.d8
            istep(i)=nint(dtime/steps(i))
            t0(i)=time
            tempold(i)=temp(i)
            rhold(i)=rho(i)
            yeold(i)=ye(i)
         enddo
         call hydro(time,ncell,x,v,
     1           u,rho,ye,f1v,f1u,f1ye,q)
      end if
c
c--set step counter
c
      do i=1,ncell
         istep(i)=nint(2.*dtime/steps(i))
      enddo
c
c--integration
c  -----------
c
c  a) Unique time step case
c  ------------------------
c
      ivar=0
      if(ivar.eq.0) then
c
c--get predictions at half time step
c
 99      continue
c 
c--gas particles
c
         dtf11=f11*steps(1)
         dumx(0)=x(0) + dtf11*v(0)
         do i=1,ncell
            dtf11=f11*steps(i)
            dumx(i)=x(i) + dtf11*v(i)
            dumv(i)=v(i) + dtf11*f1v(i)
            dumu(i)=u(i) + dtf11*f1u(i)
            dumye(i)=ye(i) + dtf11*f1ye(i)
            reset(i)=.false.
            if (dumye(i).lt.0.02) then
               dumye(i)=.02
            end if
         enddo
         dumv(ncell)=0.

c--get forces at half the time step, using predictions
c
         thalf=time + f11*steps(1)
         call hydro(thalf,ncell,dumx,dumv,
     $         dumu,rho,dumye,f2v,f2u,f2ye,q)
c
c--advance all the gas particles
c
         dtf21=f21*steps(1)
         dtf22=f22*steps(1)
         dumx(0)=x(0) + dtf21*v(0) + dtf22*v(0)
         do i=1,ncell
            dtf21=f21*steps(i)
            dtf22=f22*steps(i)
            dumx(i)=x(i) + dtf21*v(i) + dtf22*dumv(i)
            dumv(i)=v(i) + dtf21*f1v(i) + dtf22*f2v(i)
            dumu(i)=u(i) + dtf21*f1u(i) + dtf22*f2u(i)
            dumye(i)=ye(i) + dtf21*f1ye(i) + dtf22*f2ye(i)
            if (dumye(i).lt.0.02) then
               dumye(i)=.02
            end if
         enddo
         dumv(ncell)=0.
c
c     set saturation const=0
         satc=0
c
c--get forces at end of time step
c
         tfull=time + steps(1)
         call hydro(tfull,ncell,dumx,dumv,
     $         dumu,rho,dumye,f2v,f2u,f2ye,q)
c
c--estimate integration error and set time step. If reduction
c  the maximum time step is set to dtime.
c  Error checking is supressed for particles having been
c  reset.
c 
         rapmax=0.
         ermax=-1.e10
         stepmin=1.e30
         stepmax=0.
         ierr=0
         nreset=0
         denmax=0.
         tempmax=0.
         do i=1,ncell
            if(reset(i)) then
               steps(i)=1.e20
               nreset=nreset + 1
            else
               dxi=(x(i)-x(i-1))
               vsi=0.05*vsound(i)
               ysi=ye(i)
               usi=u(i)
               erx=abs(dumv(i))/dxi
               erv=abs(f1v(i)-f2v(i))/(abs(v(i))+vsi)
               eru=abs(f1u(i)-f2u(i))/u(i)
               erye=abs(f1ye(i)-f2ye(i))/ysi
c
c--get maximum error and determine time step
c
               erm=dmax1(erx,erv,eru,erye)
               if(erm.gt.ermax)then
                  if(erm.eq.erx)then
                     ierr=1
                  elseif(erm.eq.erv)then
                     ierr=4
                  elseif(erm.eq.eru)then
                     ierr=7
                  elseif(erm.eq.erye) then
                     ierr=9
                  endif
                  imax=i
                  ermax=erm
               endif 
               rap=steps(i)*erm*e1/tol + tiny
               rapmax=dmax1(rapmax,rap)
               rmod=min(2.0d0,1./dsqrt(rap))
               steps(i)=min(steps(i)*rmod,dtime,dt0)
               tempmax=max(temp(i),tempmax)
               denmax=max(rho(i),denmax)
            end if
         enddo
c         write(*,*)'max error flag: ',ierr,imax,ebetaeq(imax)
c         write(*,*)'temp,x',temp(imax),dumx(imax),dumx(imax-1)
c         write(*,*)'max temperature, max density', tempmax,denmax
c         write(*,*)'v,rho',dumv(imax),rho(imax)
c         write(*,*)'u,ye',u(imax),ye(imax)
c         write(*,*)'new: u,ye',dumu(imax),dumye(imax)
c         write(*,*)'fx',f2v(imax)
c         write(*,*)'du,dye',f2u(imax),f2ye(imax)
c
c--find minimum time step
c
         dtmin=1.e30
         do i=1,ncell
            if (dtmin.gt.steps(i)) then
c            dtmin=dmin1(dtmin,steps(i))
c               print *, steps(i), dtmin
               dtmin=steps(i)
            end if
            stepmax=dmax1(stepmax,steps(i))
            stepmin=dmin1(stepmin,steps(i))
         enddo
c         xtime=xtime*0.8d0**satc
         do i=1,ncell
            steps(i)=dtmin
c            steps(i)=min(dtmin,xtime)
            if (time.lt.1.d-7) then
               steps(i)=min(steps(i),1.d-9)
            end if
         enddo
         write(*,110)'step:t,dt,rapmax,nreset', tfull,dtmin,rapmax,
     1             nreset
  110    format(A30,3(1g12.5,1x),I3)
c
c--check for rejection
c
         if(rapmax.gt.1.5)then
            go to 99
         end if
c
c--update system
c
c         x(0)=dumx(0)
         do i=1,ncell
            x(i)=dumx(i)
            v(i)=dumv(i)
            u(i)=dumu(i)
            ye(i)=dumye(i)
            f1v(i)=f2v(i)
            f1u(i)=f2u(i)
            f1ye(i)=f2ye(i)
         enddo
         if (x(ncell).gt.rout) then 
            pr(ncell+1)=2.*p1out
         else
            pr(ncell+1)=0.1*p1out
         end if
         time=tfull
c
c--do various neutrino flagging and checking
c
         do i=1,ncell
            tempold(i)=temp(i)
            rhold(i)=rho(i)
            yeold(i)=ye(i)
         enddo
         convf=ufoe/utime
         nupk=nupk+1
c
c--start new time step
c
         if(time.lt.tnext)then
            if (mod(nupk,nups).eq.0) then
               lu=72
               open(lu,file='1dtmp',form='unformatted')
               call printout(lu)
               close (lu)
            end if
            go to 99
         end if
         return
      end if
c
      end
      subroutine unit
c************************************************************
c                                                           *
c  this routine computes the transformation between the     *
c  physical units (cgs) and the units used in the code.     *
c  And the value of physical constants in code units        *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
c
      double precision umass
      double precision usltemp, uslrho, uslu, uslp,u2slu
      double precision ud1,ut1,ue1,ueg1,uec1
      double precision gg1, arad1, bigr1
      double precision uopr, uotemp, uorho1, uotemp1, uou1
      double precision ufoe1, sigma1, sigma2, xsecnn1, xsecne1,fermi
c
      common /typef/ iextf, ieos
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common/uocean/ uopr, uotemp, uorho1, uotemp1, uou1
      common/uswest/ usltemp, uslrho, uslu, uslp, u2slu
      common /epcap/ betafac, c2cu, c3cu
      common /unit2/ utemp, utmev, ufoe, umevnuc, umeverg
      common /const/ gg, clight, arad, bigr, xsecnn, xsecne
c
      data ggcgs/6.67e-8/, avo/6.02e23/
      data aradcgs /7.565e-15/, boltzk/1.381e-16/
      data hbar/1.055e-27/, ccgs/3e10/
      data emssmev/0.511e0/, boltzmev/8.617e-11/
      data ergmev /6.2422e5/, sigma1/9d-44/, sigma2/5.6d-45/
      data c2cgs /6.15e-4/, c3cgs /5.04e-10/, fermi/1d-13/
c
c--1) work out code units:
c
c--specifie mass unit (g)
c
      umass=2d33
c
c--specifie distance unit (cm)
c
      udist=1e9
c
c--transformation factor for :
c
c 1a) density
c
      ud1=dble(umass)/dble(udist)**3
      udens=ud1
c
c 1b) time
c
c     ut1=dsqrt(dble(udist)**3/(dble(gg)*umass))
c     utime=ut1
      ut1=1.d1
      utime=ut1
c
c 1c) ergs
c
      ue1=dble(umass)*dble(udist)**2/dble(utime)**2
c--uerg will overflow
c      uerg=ue1
c
c 1d) ergs per gram
c
      ueg1=dble(udist)**2/dble(utime)**2
      uergg=ueg1
c
c 1e) ergs per cc
c
      uec1=dble(umass)/(dble(udist)*dble(utime)**2)
      uergcc=uec1
c
c--2) constants
c
c 2a) gravitational
c
      gg1=dble(ggcgs)*umass*dble(utime)**2/(dble(udist)**3)
      gg=gg1
c
c 2aa) velocity of light
c
      clight=dble(ccgs*utime/udist)

c
c 2b) Stefan Boltzmann (note that T unit is 1e9 K here)
c
      utemp=1e9
      arad1=dble(aradcgs)/ue1*dble(udist)**3*dble(utemp)**4
      arad=arad1
c
c 2c) Perfect Gas "R" constant
c
      bigr1=dble(avo)*dble(boltzk)/ue1*umass*dble(utemp)
      bigr=bigr1
c
c 2d) nucleon+nu x-section/4pi, in udist**2/umass/Mev**2
c
      xsecnn1=sigma1*umass*dble(avo)/dble(udist*udist)
      xsecnn=xsecnn1/(4d0*3.14159d0)
c
c 2e) e+nu x-section/4pi, in Enu(Mev)*udist**2/umass/Mev**2/utemp**4
c
      xsecne1=sigma2*umass*dble(ergmev)*(dble(utemp)*dble(boltzk))**4/
     1     (dble(udist)**2*ud1*3.14159d0*(dble(hbar)*dble(ccgs))**3)
      xsecne=xsecne1/(4d0*3.14159d0)
c
c--3a) Conversion to the Ocean eos units from code units:
c     in the ocean eos; density unit=1e7g/cc
c                       temperature unit= 1e9K
c                       energy/mass=1.0e17cgs
c                       energy/vol =1.0e24cgs
c
      uotemp=1d9/dble(utemp)
      uotemp1=1.d0/uotemp
      uopr=1.d24/dble(uergcc)
      uorho1=ud1/1.d7
      uou1=dble(uergg)/1.d17




c
c--3b) Conversion to the Swesty-Lattimer units from code units:
c     
      usltemp=utemp*boltzmev
      uslrho=ud1*avo*fermi**3      
      uslp=uec1*fermi**3/1.602d-6
      if (ieos.eq.3) then
c--if u is internal energy
         uslu=dble(ergmev)*dble(uergg)/dble(avo)
         u2slu=dble(uergg)/dble(utemp)/(dble(avo)*dble(boltzk))
      elseif (ieos.eq.4) then
c--if u is specific entropy
         uslu=dble(uergg)/dble(utemp)/(dble(avo)*dble(boltzk))
         u2slu=dble(ergmev)*dble(uergg)/dble(avo)
      endif
      print *,'ieos,uslu',ieos,uslu
c
c--4) common unit2 stuff
c
c 4a) temp code unit in Mev
c
      utmev=utemp*boltzmev
c
c 4b) energy code unit in foes
c
      ufoe1=ue1/1d51
      ufoe=ufoe1
c
c 4c) Mev/nucleon in code units
c
      umevnuc=avo/(ergmev*uergg)
c
c 4e) Mev to errgs times avogadro's number
c
      umeverg=avo*1.602e-6
c
c--5) epcapture betafac, c2, and c3. 
c
      betafac=emssmev/utmev
      c2cu=c2cgs*utime
      c3cu=c3cgs*utime*ergmev
c
      write(*,100)umass,udist,udens,utime,uergg,uergcc
  100 format(1x,5(1pe12.5,1x),/,1x,2(1pe12.5,1x))
c
   99 return
      end

c © 2020. Triad National Security, LLC. All rights reserved.

c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos

c National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.

c Department of Energy/National Nuclear Security Administration. All rights in the program are

c reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear

c Security Administration. The Government is granted for itself and others acting on its behalf a

c nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare

c derivative works, distribute copies to the public, perform publicly and display publicly, and to permit

c others to do so.
c This program is open source under the BSD-3 License.

c Redistribution and use in source and binary forms, with or without modification, are permitted

c provided that the following conditions are met:
c 1. Redistributions of source code must retain the above copyright notice, this list of conditions and

c the following disclaimer.

 

c 2.Redistributions in binary form must reproduce the above copyright notice, this list of conditions

c and the following disclaimer in the documentation and/or other materials provided with the

c distribution.

c 3.Neither the name of the copyright holder nor the names of its contributors may be used to endorse

c or promote products derived from this software without specific prior written permission.

c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS

c IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE

c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR

c PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR

c CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,

c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,

c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;

c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,

c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR

c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF

c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

























