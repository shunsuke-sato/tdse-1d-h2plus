!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dipole(dipole_m,norm_m)
  use global_variables
  implicit none
!  real(dp),intent(out) :: dipole_t ,norm_t
  real(dp) :: dipole_m, norm_m
  integer :: ix,iy


  dipole_m = 0d0; norm_m = 0d0
  do iy = 0,NR
    do ix = 0,Nx
      dipole_m = dipole_m + xn(ix)*abs(zwfn(ix,iy))**2
      norm_m = norm_m + abs(zwfn(ix,iy))**2
    end do
  end do

  dipole_m = dipole_m*dx*dr
  norm_m = norm_m*dx*dr


  return
end subroutine dipole
