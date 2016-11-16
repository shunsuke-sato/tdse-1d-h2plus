!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dipole(dipole_t ,norm_t)
  use global_variables
  implicit none
  real(dp),intent(out) :: dipole_t ,norm_t
  integer :: ix,iy
  real(dp) :: tmp


  dipole_t = 0d0; norm_t = 0d0
  do iy = 0,NR
    do ix = 0,Nx
      dipole_t = dipole_t + xn(ix)*abs(zwfn(ix,iy))**2
      norm_t = norm_t + abs(zwfn(ix,iy))**2
    end do
  end do

  dipole_t = dipole_t*dx*dr
  norm_t = norm_t*dx*dr


  return
end subroutine wfn_rho
