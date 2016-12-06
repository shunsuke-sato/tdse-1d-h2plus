!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation_RT
  use global_variables
  implicit none
  integer :: ix,iy


  allocate(zwfn(0:Nx,0:NR))
  allocate(dipole_t(0:Nt_iter),norm_t(0:Nt_iter))
!Initial condition
!  zwfn = wfn ! GS wavefunction

!!momentum kick
!  do iy=0,NR
!    do ix=0,Nx
!      zwfn(ix,iy) = exp(zI*kick_mom*xn(ix))*wfn(ix,iy) 
!    end do
!  end do

!!quadrupole kick
  do iy=0,NR
    do ix=0,Nx
      zwfn(ix,iy) = exp(zI*kick_mom*xn(ix)**2)*wfn(ix,iy) 
    end do
  end do

  return
end subroutine preparation_RT
