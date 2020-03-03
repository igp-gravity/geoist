! Copyright 2019 Pascal Audet
!
! This file is part of PlateFlex.
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!===========================================================================
!
! MODULE conf_cpwt
!
! Configuration module defining global variables used in modules to 
! interface with the Python codes.
!
!===========================================================================

      MODULE conf_cpwt

      IMPLICIT NONE

      REAL, PARAMETER :: pi = 3.141592653589793
      INTEGER, PARAMETER :: na = 11
!
! Wavelet parameter
!
      REAL :: k0

      END MODULE conf_cpwt

!===========================================================================
!
! MODULE cpwt
!
! Contains subroutines to calculate the wavelet transform, scalogram and 
! cross-spectral quantities with errors
!
!===========================================================================

      MODULE cpwt

      CONTAINS

!---------------------------------------------------------------------------
! Subroutine wlet_transform
!
! Subroutine to compute the complex-valued wavelet transform coefficients
! obtained from directional Morlet wavelets
!---------------------------------------------------------------------------

      SUBROUTINE wlet_transform(grid, nx, ny, dx, dy, kf, ns, wt_grid)     

      USE conf_cpwt

      IMPLICIT NONE

      INTEGER :: nx, ny, nnx, nny, ns, nn(2)
      REAL    :: dx, dy
      REAL    :: grid(nx,ny), grid_mirror(2*nx, 2*ny), kf(ns)
      COMPLEX :: wt_grid(nx,ny,na,ns)

      INTEGER :: is, ia
      REAL    :: da, scales, angle, lam, kk
      CHARACTER(LEN=1) :: trail
      CHARACTER(LEN=30) :: progress

! Allocatable work arrays
      REAL, ALLOCATABLE :: grid_pad(:,:), kx(:), ky(:)
      COMPLEX, ALLOCATABLE :: ft_grid(:,:), ft2_grid(:,:), daughter(:,:) 
!
! Python bindings
!
!f2py REAL, intent(in) :: grid
!f2py INTEGER, intent(hide),depend(grid) :: nx=shape(grid,0), ny=shape(grid,1)
!f2py REAL, intent(in) :: dx, dy
!f2py REAL, intent(in) :: kf
!f2py INTEGER, intent(hide),depend(kf) :: ns=shape(kf,0)
!f2py COMPLEX, intent(out) :: wt_grid

        CALL npow2(2*nx, nnx)
        CALL npow2(2*ny, nny)
        ALLOCATE(grid_pad(nnx,nny))
! 
! Mirror data
!
        grid_mirror = 0.
        grid_mirror(1:nx,1:ny) = grid
        CALL mirror_data(grid_mirror,nx,ny,2*nx,2*ny)
! 
! Taper data
!
        ! CALL taper_data(grid_mirror,2*nx,2*ny,6)
!
! Pad with zeros and remove mean value
!
        grid_pad = 0.
        grid_pad(1:2*nx,1:2*ny) = grid_mirror - sum(grid_mirror)/2./2./nx/ny
!
! Complex form
!
        ALLOCATE(ft_grid(nnx,nny), ft2_grid(nnx,nny), daughter(nnx,nny))
        ft_grid = CMPLX(grid_pad)
!
! Fourier transform
!
        nn(1) = nnx
        nn(2) = nny
        CALL fourn(ft_grid,nn,2,1)
!
! Define wavenumbers
!
        ALLOCATE(kx(nnx), ky(nny))
        CALL defk(nnx,nny,dx,dy,kx,ky)
!
! Define angle parameters as per Kirby (2005)
!
        da = 2.*SQRT(-2.*LOG(0.75))/k0
!_______________________________________________
!     
!     MAIN WAVELET TANSFORM LOOP
!_______________________________________________
!
!
! Loop through scales
!
        trail = '-'
        progress = ''
        WRITE(*, fmt="(1x,a,i2,a)", advance="no")'#loops = ',ns,':'
        DO is = 1,ns
!
! Define wavenumbers and scales
!
          trail = '_'
          kk = kf(is)*1.e3
          lam = 2.*pi/kk
          scales = k0/kk
          ! WRITE(*, fmt="(1x,a)", advance="no")TRIM(progress)//trail
          WRITE(*, fmt="(1x,i2)", advance="no")is
          progress = TRIM(progress)//trail
!
! Loop through angles
!
          DO ia = 1,na
!
! Define angle increments
! 
            angle = REAL(ia-1)*da-pi/2.e0
!
! Calculate daughter wavelet
!
            CALL wave_function(nnx,nny,kx,ky,k0,scales,angle,daughter)
!
! Compute wavelet transform in Fourier space
!
            daughter = daughter*CONJG(ft_grid)
!
! Back to physical space
!
            CALL fourn(daughter,nn,2,-1)
!
! Normalization
!
            ft2_grid = (daughter)/FLOAT(nnx*nny)
!
! Quadrant swap
!
            CALL cfftswap(ft2_grid,nnx,nny)
!
! Store array
!
            wt_grid(:,:,ia,is) = ft2_grid(1:nx,1:ny)

          END DO
        END DO
        WRITE(*,*)

      END SUBROUTINE wlet_transform

!---------------------------------------------------------------------------
! Subroutine wlet_scalogram
!
! Subroutine to compute the real-valued wavelet scalogram obtained from 
! azimuthal-averaging of wavelet transform coefficients
!---------------------------------------------------------------------------

      SUBROUTINE wlet_scalogram(wl_trans, nx, ny, ns, wl_sg, ewl_sg)

      USE conf_cpwt

      IMPLICIT NONE

      INTEGER :: nx, ny, ns
      COMPLEX :: wl_trans(nx,ny,na,ns)
      REAL    :: wl_sg(nx,ny,ns), ewl_sg(nx,ny,ns)
      REAL    :: sg(nx,ny,na,ns)

      INTEGER :: is, ia
!
! Python bindings
!
!f2py COMPLEX, intent(in) :: wl_trans
!f2py INTEGER, intent(hide),depend(wl_trans) :: nx=shape(wl_trans,0)
!f2py INTEGER, intent(hide),depend(wl_trans) :: ny=shape(wl_trans,1)
!f2py INTEGER, intent(hide),depend(wl_trans) :: ns=shape(wl_trans,3)
!f2py REAL, intent(out) :: wl_sg, ewl_sg

        wl_sg = 0.e0
        DO is = 1,ns
          DO ia = 1,na
            sg(:,:,ia,is) = CABS(wl_trans(:,:,ia,is)*CONJG(wl_trans(:,:,ia,is)))
!
! Azimuthal average
!
            wl_sg(:,:,is) = wl_sg(:,:,is) + sg(:,:,ia,is)
          END DO
        END DO
!
! Normalization
!
        wl_sg = wl_sg/FLOAT(na)
!
! Jackknife estimation
!
        CALL jackknife_1d_sgram(sg,nx,ny,na,ns,ewl_sg)

        RETURN
      END SUBROUTINE wlet_scalogram

!---------------------------------------------------------------------------
! Subroutine wlet_xscalogram
!
! Subroutine to compute the complex-valued wavelet cross-scalogram obtained 
! from azimuthal-averaging of wavelet transform coefficients
!---------------------------------------------------------------------------

      SUBROUTINE wlet_xscalogram(wl_trans1, wl_trans2, nx, ny, ns, wl_xsg, ewl_xsg)

      USE conf_cpwt

      IMPLICIT NONE

      INTEGER :: nx, ny, ns
      COMPLEX :: wl_trans1(nx,ny,na,ns), wl_trans2(nx,ny,na,ns)
      COMPLEX :: wl_xsg(nx,ny,ns), xsg(nx,ny,na,ns)
      REAL :: ewl_xsg(nx,ny,ns), ewl_xsg_real(nx,ny,ns), ewl_xsg_imag(nx,ny,ns)

      INTEGER :: is, ia
!
! Python bindings
!
!f2py COMPLEX, intent(in) :: wl_trans1, wl_trans2
!f2py INTEGER, intent(hide),depend(wl_trans1) :: nx=shape(wl_trans1,0)
!f2py INTEGER, intent(hide),depend(wl_trans1) :: ny=shape(wl_trans1,1)
!f2py INTEGER, intent(hide),depend(wl_trans1) :: ns=shape(wl_trans1,3)
!f2py REAL, intent(out) :: wl_xsg, ewl_xsg

        wl_xsg = CMPLX(0.e0,0.e0)
        DO is = 1,ns
          DO ia = 1, na
            xsg(:,:,ia,is) = wl_trans1(:,:,ia,is)*CONJG(wl_trans2(:,:,ia,is))
!
! Azimuthal average
!
            wl_xsg(:,:,is) = wl_xsg(:,:,is) + xsg(:,:,ia,is)
          END DO
        END DO
!
! Normalization
!
        wl_xsg = wl_xsg/FLOAT(na)
!
! Jackknife estimation
!
        CALL jackknife_1d_sgram(REAL(xsg),nx,ny,na,ns,ewl_xsg_real)
        CALL jackknife_1d_sgram(IMAG(xsg),nx,ny,na,ns,ewl_xsg_imag)
        ewl_xsg = SQRT(ewl_xsg_real**2 + ewl_xsg_imag**2)

        RETURN
      END SUBROUTINE wlet_xscalogram

!---------------------------------------------------------------------------
! Subroutine wlet_admit_coh
!
! Subroutine to compute the complex-valued wavelet admittance and real-valued
! wavelet coherence obtained from azimuthal-averaging of wavelet transform 
! coefficients
!---------------------------------------------------------------------------

      SUBROUTINE wlet_admit_coh(wl_trans1, wl_trans2, nx, ny, ns, &
        wl_admit, ewl_admit, wl_coh, ewl_coh)

      USE conf_cpwt

      IMPLICIT NONE

      INTEGER :: nx, ny, ns
      COMPLEX :: wl_trans1(nx,ny,na,ns), wl_trans2(nx,ny,na,ns)
      REAL    :: wl_admit(nx,ny,ns), wl_coh(nx,ny,ns)
      REAL    :: ewl_admit(nx,ny,ns), ewl_coh(nx,ny,ns)

      REAL    :: sg1(nx,ny,na,ns), sg2(nx,ny,na,ns)
      COMPLEX :: xsg(nx,ny,na,ns)

      REAL    :: wl_sg1(nx,ny,ns), wl_sg2(nx,ny,ns)
      COMPLEX :: wl_xsg(nx,ny,ns)
      INTEGER :: is, ia
!
! Python bindings
!
!f2py COMPLEX, intent(in) :: wl_trans1, wl_trans2
!f2py INTEGER, intent(hide),depend(wl_trans1) :: nx=shape(wl_trans1,0)
!f2py INTEGER, intent(hide),depend(wl_trans1) :: ny=shape(wl_trans1,1)
!f2py INTEGER, intent(hide),depend(wl_trans1) :: ns=shape(wl_trans1,3)
!f2py REAL, intent(out) :: wl_admit, ewl_admit, wl_coh, ewl_coh

        DO is = 1,ns
          DO ia = 1, na

            sg1(:,:,ia,is) = CABS(wl_trans1(:,:,ia,is)*CONJG(wl_trans1(:,:,ia,is)))
            sg2(:,:,ia,is) = CABS(wl_trans2(:,:,ia,is)*CONJG(wl_trans2(:,:,ia,is)))
            xsg(:,:,ia,is) = wl_trans1(:,:,ia,is)*CONJG(wl_trans2(:,:,ia,is))
!
! Azimuthal average
!
            wl_sg1(:,:,is) = wl_sg1(:,:,is) + sg1(:,:,ia,is)
            wl_sg2(:,:,is) = wl_sg2(:,:,is) + sg2(:,:,ia,is)
            wl_xsg(:,:,is) = wl_xsg(:,:,is) + xsg(:,:,ia,is)

          END DO
        END DO
!
! Normalization
!
        wl_sg1 = wl_sg1/FLOAT(na)
        wl_sg2 = wl_sg2/FLOAT(na)
        wl_xsg = wl_xsg/FLOAT(na)
!
! Admittance and coherence (squared-real coherency) 
!
        wl_admit = REAL(wl_xsg)/wl_sg1
        wl_coh = REAL(wl_xsg)**2/(wl_sg1*wl_sg2)
        ! wl_coh = wl_xsg/SQRT(wl_sg1*wl_sg2)
!
! Jackknife estimation
!
        PRINT*,'Calculation jackknife error on admittance and coherence'
        CALL jackknife_admit_coh(xsg,sg1,sg2,nx,ny,na,ns,ewl_admit,ewl_coh)

        RETURN
      END SUBROUTINE wlet_admit_coh

      END MODULE cpwt