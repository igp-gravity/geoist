!===============================================================
! SUBROUTINES cpwt_sub
!
!   Provies subroutines to CPWT program 
!
!   (See MakeFILE for compilation)
!
!===============================================================
!   Pascal Audet
!   Department of Earth and Planetary Science
!   U! Berkeley
!   paudet@berkeley.edu      
!
!   July 2009
!================================================================
   
!_____________________________________________________________
!
!     SUBROUTINE LIBRARY
!_____________________________________________________________
!

!--------------------------------------------
!     Function determines smallest power of 2
!--------------------------------------------
      SUBROUTINE npow2(n, n2)
      n2=1
      DO 
        n2=2*n2
        IF(n2.ge.n)RETURN
      END DO
      RETURN 
      END SUBROUTINE npow2


!------------------------------------------------
!     Subroutine defk
!
!     Computes wavenumbers for given nx,ny,dx,dy
!------------------------------------------------
      SUBROUTINE defk(nx,ny,dx,dy,kx,ky)
      IMPLICIT NONE
      REAL :: pi = 3.141592653589793
      INTEGER :: nx,ny,ix,iy
      REAL :: kx(nx),ky(ny),dx,dy,dkx,dky

      dkx=2*pi/(FLOAT(nx)*dx)
      dky=2*pi/(FLOAT(ny)*dy)
      DO iy=1,ny/2+1
        DO ix=1,nx/2+1
          kx(ix)=FLOAT(ix-1)*dkx
          ky(iy)=FLOAT(iy-1)*dky
        END DO
      END DO
      DO iy=ny/2+2,ny
        DO ix=nx/2+2,nx
          kx(ix)=-kx(nx-ix+2)
          ky(iy)=-ky(ny-iy+2)
        END DO
      END DO
      RETURN
      END SUBROUTINE defk


!---------------------------------------------
!     Subroutine wave_function
!
!     !omputes the daughter Morlet wavelet 
!     in Fourier space (kx,ky) with a 
!     preferred orientation in wavevector space
!     given by 'a', at scale 's'
!---------------------------------------------
      SUBROUTINE wave_function(nx, ny, kx, ky, k0, s, a, daughter)
      IMPLICIT NONE
      REAL :: pi = 3.141592653589793
      INTEGER :: nx, ny, ix, iy
      REAL :: k0, s, a 
      REAL :: norm, expntx, expnty, expnt
      REAL :: kx(nx), ky(ny), kx0, ky0
      COMPLEX :: daughter(nx,ny)

!
! Morlet wavelet in Fourier domain
!
      kx0=k0*cos(a)
      ky0=k0*sin(a)
      norm=s/SQRT(pi)
      DO iy=1,ny
      DO ix=1,nx
        expntx = -0.5*((s*kx(ix) - kx0)**2)
        expnty = -0.5*((s*ky(iy) - ky0)**2)
        daughter(ix,iy) = CMPLX(norm*EXP(expntx)*EXP(expnty))
        expnt = -0.5*((kx(ix)**2 + ky(iy)**2) + (kx0**2 + ky0**2))
        daughter(ix,iy) = daughter(ix,iy) - CMPLX(norm*EXP(expnt))
      END DO
      END DO

      RETURN
      END SUBROUTINE wave_function


!----------------------------------------------------------
!     Subroutine jackknife_1d_admit
!
!     Calculates Jackknife error on admittance.
!----------------------------------------------------------
      SUBROUTINE jackknife_1d_sgram(ss,nx,ny,na,ns,ess)
      IMPLICIT NONE
      INTEGER :: nx,ny,na,ns,ix,iy,ia,ia2,i_k,is
      REAL :: ss(nx,ny,na,ns),p1

      REAL :: mss,dss,ssj(na-1),ssj2(na),ess(nx,ny,ns)

      DO is=1,ns
      DO iy=1,ny
      DO ix=1,nx
      DO ia=1,na

        i_k=0
        p1=0.

        DO ia2=1,na

          IF (ia2==ia) THEN
            CYCLE
          END IF

          i_k=i_k+1

          p1=p1+ss(ix,iy,ia2,is)

          IF (p1.ne.0) THEN
            ssj(i_k)=p1
          ELSE
            ssj(i_k)=0.
          END IF

        END DO

        ssj2(ia)=SUM(ssj)/FLOAT(na-1)

      END DO

      mss=SUM(ssj2)/FLOAT(na)
      dss=0.

      DO ia=1,na
        dss=dss+(ssj2(ia)-mss)**2
      END DO

      dss=(SQRT(dss*FLOAT(na-1)/FLOAT(na)))
      IF (ISNAN(dss)) dss=0.

      ess(ix,iy,is)=dss
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE jackknife_1d_sgram


!----------------------------------------------------------
!     Subroutine jackknife_1d_admit
!
!     Calculates Jackknife error on admittance.
!----------------------------------------------------------
      SUBROUTINE jackknife_1d_admit(sxp,ss1,nx,ny,na,ns,sad)
      IMPLICIT NONE
      INTEGER :: nx,ny,na,ns,ix,iy,ia,ia2,i_k,is
      REAL :: ss1(nx,ny,na,ns),p1
      REAL :: sxp(nx,ny,na,ns),xp

      REAL :: mad,dad,ad(na-1),ad2(na),sad(nx,ny,ns)

      DO is=1,ns
      DO iy=1,ny
      DO ix=1,nx
      DO ia=1,na

        i_k=0
        p1=0.
        xp=(0.,0.)

        DO ia2=1,na

          IF (ia2==ia) THEN
            CYCLE
          END IF

          i_k=i_k+1

          p1=p1+ss1(ix,iy,ia2,is)
          xp=xp+sxp(ix,iy,ia2,is)

          IF (p1.ne.0) THEN
            ad(i_k)=xp/p1
          ELSE
            ad(i_k)=0.
          END IF

        END DO

        ad2(ia)=SUM(ad)/FLOAT(na-1)

      END DO

      mad=SUM(ad2)/FLOAT(na)
      dad=0.

      DO ia=1,na
        dad=dad+(ad2(ia)-mad)**2
      END DO

      dad=(SQRT(dad*FLOAT(na-1)/FLOAT(na)))
      IF (ISNAN(dad)) dad=0.

      sad(ix,iy,is)=dad
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE jackknife_1d_admit


!----------------------------------------------------------
!     Subroutine jackknife_1d_corr
!
!     !alculates Jackknife error on coherency.
!----------------------------------------------------------
      SUBROUTINE jackknife_1d_corr(sxp,ss1,ss2,nx,ny,na,ns,sco)
      IMPLICIT NONE
      INTEGER :: nx,ny,na,ns,ix,iy,ia,ia2,i_k,is
      REAL :: ss1(nx,ny,na,ns),ss2(nx,ny,na,ns),p1,p2
      REAL :: sxp(nx,ny,na,ns),xp

      REAL :: mco,dco,co(na-1),co2(na),sco(nx,ny,ns)

      DO is=1,ns
      DO iy=1,ny
      DO ix=1,nx
      DO ia=1,na

        i_k=0
        p1=0.
        p2=0.
        xp=(0.,0.)

        DO ia2=1,na

          IF (ia2==ia) THEN
            CYCLE
          END IF

          i_k=i_k+1

          p1=p1+ss1(ix,iy,ia2,is)
          p2=p2+ss2(ix,iy,ia2,is)
          xp=xp+sxp(ix,iy,ia2,is)

          IF ((p1.ne.0).and.(p2.ne.0)) THEN
            co(i_k)=xp/SQRT(p1*p2)
          ELSE
            co(i_k)=0.
          END IF

        END DO

        co2(ia)=SUM(co)/FLOAT(na-1)

      END DO

      mco=SUM(co2)/FLOAT(na)
      dco=0.

      DO ia=1,na
        dco=dco+(co2(ia)-mco)**2
      END DO

      dco=(SQRT(dco*FLOAT(na-1)/FLOAT(na)))
      IF (ISNAN(dco)) dco=0.

      sco(ix,iy,is)=dco
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE jackknife_1d_corr


!----------------------------------------------------------
!     Subroutine jackknife_1d_corr
!
!     Calculates Jackknife error on coherence.
!----------------------------------------------------------
      SUBROUTINE jackknife_admit_coh(sxp,ss1,ss2,nx,ny,na,ns,sad,sco)
      IMPLICIT NONE
      INTEGER :: nx,ny,na,ns,ix,iy,ia,ia2,i_k,is
      REAL :: ss1(nx,ny,na,ns),ss2(nx,ny,na,ns),p1,p2
      COMPLEX :: sxp(nx,ny,na,ns),xp

      REAL :: mco,dco,co(na-1),co2(na),sco(nx,ny,ns)
      REAL :: mad,dad,ad(na-1),ad2(na),sad(nx,ny,ns)

      DO is=1,ns
      DO iy=1,ny
      DO ix=1,nx
      DO ia=1,na

        i_k=0
        p1=0.
        p2=0.
        xp=CMPLX(0.,0.)

        DO ia2=1,na

          IF (ia2==ia) THEN
            CYCLE
          END IF

          i_k=i_k+1

          p1=p1+ss1(ix,iy,ia2,is)
          p2=p2+ss2(ix,iy,ia2,is)
          xp=xp+sxp(ix,iy,ia2,is)

          IF ((p1.ne.0).and.(p2.ne.0)) THEN
            co(i_k)=REAL(xp)**2/(p1*p2)
            ad(i_k)=REAL(xp)/p1
          ELSE
            co(i_k)=0.
            ad(i_k)=0.
          END IF

        END DO

        co2(ia)=SUM(co)/FLOAT(na-1)
        ad2(ia)=SUM(ad)/FLOAT(na-1)

      END DO

      mco=SUM(co2)/FLOAT(na)
      mad=SUM(ad2)/FLOAT(na)
      dco=0.
      dad=0.

      DO ia=1,na
        dco=dco+(co2(ia)-mco)**2
        dad=dad+(ad2(ia)-mad)**2
      END DO

      dco=(SQRT(dco*FLOAT(na-1)/FLOAT(na)))
      dad=(SQRT(dad*FLOAT(na-1)/FLOAT(na)))
      IF (ISNAN(dco)) dco=0.
      IF (ISNAN(dad)) dad=0.

      sco(ix,iy,is)=dco
      sad(ix,iy,is)=dad
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE jackknife_admit_coh


!----------------------------------------------------------
!     Subroutine mean_std
!
!     !omputes mean and standard deviation from signal
!----------------------------------------------------------
      SUBROUTINE mean_std(hwl,atot,mhwl,shwl)
      IMPLICIT NONE
      INTEGER :: atot,a
      REAL :: hwl(atot),mhwl,shwl

      mhwl=0.
      shwl=0.
      DO a=1,atot
        mhwl=mhwl+hwl(a)
      END DO
      mhwl=mhwl/REAL(atot)
      DO a=1,atot
        shwl=shwl+(hwl(a)-mhwl)**2
      END DO
      shwl=SQRT(shwl/REAL(atot))
      RETURN
      END SUBROUTINE mean_std


!---------------------------------------------
!     Subroutine taper_data
!
!     Puts half taper function on four sides
!     of matrix
!---------------------------------------------
      SUBROUTINE taper_data(dat,maxx,maxy,maxt)
      IMPLICIT NONE
      REAL :: pi = 3.141592653589793
      INTEGER :: maxt,it,ix,iy,maxx,maxy
      REAL :: taperx(maxx/2),tapery(maxy/2)
      REAL :: dat(maxx,maxy)

      taperx=0.
      tapery=0.
      DO it=1,maxt
        tapery(it)=0.5*(1.-COS(pi*(it-1)/maxt))
      END DO
      DO it=1,maxt
        taperx(it)=0.5*(1.-COS(pi*(it-1)/maxt))
      END DO
      DO iy=1,maxy
         DO it=1,maxt
            dat(it,iy)=taperx(it)*dat(it,iy)
            dat(maxx+1-it,iy)=taperx(it)*dat(maxx+1-it,iy)
         END DO
      END DO
      DO ix=1,maxx
         DO it=1,maxt
            dat(ix,it)=tapery(it)*dat(ix,it)
            dat(ix,maxy+1-it)=tapery(it)*dat(ix,maxy+1-it)
         END DO
      END DO
      RETURN
      END SUBROUTINE taper_data


!---------------------------------------------
!     Subroutine mirror_data
!
!     Mirrors 2D data along edges
!---------------------------------------------
      SUBROUTINE mirror_data(dat,nx,ny,nnx,nny)
      IMPLICIT NONE
      INTEGER :: nx,ny,nnx,nny,ix,iy
      REAL :: dat(nnx,nny)

      DO iy=1,ny
      DO ix=1,nx
        dat(ix+nx,iy+ny)=dat(nx-ix+1,ny-iy+1)
        dat(ix+nx,iy)=dat(nx-ix+1,iy)
        dat(ix,iy+ny)=dat(ix,ny-iy+1)
      END DO
      END DO

      RETURN
      END SUBROUTINE mirror_data


!-------------------------------------------------
!     Subroutine cfftswap
!
!     Swaps both directions of a matrix
!-------------------------------------------------
      SUBROUTINE cfftswap(spec,nx,ny)
      IMPLICIT NONE
      COMPLEX :: spec(nx,ny),tmp(nx,ny)
      INTEGER :: nx,ny,ix,iy

      tmp=spec
      DO iy=1,ny
      DO ix=1,nx
        spec(nx-ix+1,ny-iy+1)=tmp(ix,iy)
      END DO
      END DO
      RETURN
      END SUBROUTINE cfftswap


!-------------------------------------------------
!     Subroutine fftswap
!
!     Swaps both directions of a matrix
!-------------------------------------------------
      SUBROUTINE fftswap(spec,nx,ny)
      IMPLICIT NONE
      REAL :: spec(nx,ny),tmp(nx,ny)
      INTEGER :: nx,ny,ix,iy

      tmp=spec
      DO iy=1,ny
      DO ix=1,nx
        spec(nx-ix+1,ny-iy+1)=tmp(ix,iy)
      END DO
      END DO
      RETURN
      END SUBROUTINE fftswap


!-------------------------------------------------
!     Subroutine cfftshift
!
!     Swaps both directions of a matrix
!-------------------------------------------------
      SUBROUTINE cfftshift(spec,nx,ny)
      IMPLICIT NONE
      COMPLEX :: spec(nx,ny),tmp(nx,ny)
      INTEGER :: nx,ny,ix,iy

      tmp=spec
      DO iy=1,ny/2
      DO ix=1,nx/2
        spec(ix,iy)=tmp(nx/2+ix,ny/2+iy)
        spec(nx/2+ix,iy)=tmp(ix,ny/2+iy)
        spec(nx/2+ix,ny/2+iy)=tmp(ix,iy)
        spec(ix,ny/2+iy)=tmp(nx/2+ix,iy)
      END DO
      END DO
      RETURN
      END SUBROUTINE cfftshift


!-------------------------------------------------
!     Subroutine fourn
!
!     !alculates Fast Fourier Transform
!     From Numerical Recipes textbook
!-------------------------------------------------

      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,&
          k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do idim=1,ndim
        ntot=ntot*nn(idim)
      end do
      nprev=1
      do idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do i1=i2,i2+ip1-2,2
              do i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
              end do
            end do
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
        end do
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,ip1
            do i1=i3,i3+ip1-2,2
              do i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
              end do
            end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
          end do
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
      end do
      return
      END SUBROUTINE fourn


