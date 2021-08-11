      program setup
**************************************************************
c                                                             *
c  This subroutine prints out all the results.                *
c                                                             *
c**************************************************************
c
      implicit none
c
      integer idim,idim1,iqn
      parameter (idim=4000)
      parameter (idim1=idim+1)
      parameter (iqn=17)
c
      real ycc(idim,iqn),yccave(iqn)
c
      double precision dq(idim),t
      double precision x(0:idim),v(0:idim),f(0:idim)
      double precision u(idim),rho(idim),ye(idim),q(idim)
      double precision pr(idim1), vsound(idim), u2(idim), vsmax
      double precision dcore, xmcore
      double precision rb,rho0
      double precision temp(idim)
      double precision deltam(idim), abar(idim)
      double precision uint(idim),s(idim)
      double precision lc,pi43,mtot
c
      integer nc,nqn,i,j
c
      pi43=4.0*3.14159258/3.0
c
      open(69,file='new',form='unformatted')
      open(79,file="setup_data.dat", status='new')
103   format(12(1pe12.4))
104   format(I5, 4(1pe12.4))
      nc=200
      nqn=17
      rb=1.d13/1.d9
      xmcore=1.4
      x(0)=rb
      v(0)=0.
      rho0=5.d-10/2.d6
      mtot=0.
      
      do i=1,nc
c3.6d9/1.d8
c+(1.008)**dble(i)*rb
         if (mtot.lt.0.001) then
            x(i)=1.04*x(i-1)
            rho(i)=rho0*(x(0)/x(i))**2
            v(i)=6.d9/1.d8*(x(i)/x(0))
            uint(i)=1.d10*(rho0/rho(i))/1.d16
            t=x(i)/v(i)
         else
            x(i)=1.15*x(i-1)
            rho(i)=max(10./6.d23/2.d6,rho0*(x(0)/x(i))**8)
            uint(i)=1.0*(rho0/rho(i))/1.d16
            v(i)=0.
         end if
         deltam(i)=pi43*rho(i)*(x(i)**3-x(i-1)**3)
         q(i)=0.
         dq(i)=0.
         mtot=mtot+deltam(i)

         pr(i)=(5.0/3.0-1.0)*rho(i)*uint(i)
         temp(i)=(uint(i)*1.d16/7.5657d-15/
     $        2.d6/rho(i))**0.25/1.d9
         s(i)=temp(i)**3/rho(i)

         abar(i)=1.
         ye(i)=0.5
         if (mod(i,100).eq.0) then
            print *, i,x(i),temp(i),rho(i),uint(i),
     $           mtot
         end if
      end do
      pr(nc+1)=0.1*pr(nc)
      
      lc=1.
c
      write(79, 104) nc, t, xmcore, rb, lc
      do i=1,nc
        write(79, 103) x(i), v(i), q(i), dq(i), uint(i), deltam(i),
     $    abar(i), rho(i), temp(i), ye(i), pr(i), s(i)

      end do
      write(69)nc,t,xmcore,rb,lc,
     1   (x(i),i=0,nc),(v(i),i=0,nc),(q(i),i=1,nc),(dq(i),i=1,nc),
     2      (uint(i),i=1,nc),(deltam(i),i=1,nc),(abar(i),i=1,nc),
     3      (rho(i),i=1,nc),(temp(i),i=1,nc),(ye(i),i=1,nc),
     7      (pr(i),i=1,nc),(s(i),i=1,nc),
     $     ((ycc(i,j),j=1,nqn),i=1,nc)
c
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
