!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine quadrupole(qd_pole_m,norm_m)
  use global_variables
  implicit none
!  real(dp),intent(out) :: dipole_t ,norm_t
  real(dp) :: qd_pole_m, norm_m
  integer :: ix,iy


  qd_pole_m = 0d0; norm_m = 0d0
  do iy = 0,NR
    do ix = 0,Nx
      qd_pole_m = qd_pole_m + Rn(iy)**2*abs(zwfn(ix,iy))**2
      norm_m = norm_m + abs(zwfn(ix,iy))**2
    end do
  end do

  qd_pole_m = qd_pole_m*dx*dr
  norm_m = norm_m*dx*dr


  return
end subroutine quadrupole
