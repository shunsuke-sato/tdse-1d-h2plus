!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dipole(dipole ,norm)
  use global_variables
  implicit none
  real(dp),intent(out) :: dipole,norm
  integer :: ix,iy
  real(dp) :: tmp


  dipole = 0d0; norm = 0d0
  do iy = 0,NR
    do ix = 0,Nx
      dipole = dipole + xn(ix)*abs(zwfn(ix,iy))**2
      norm = norm + xn(ix)*abs(zwfn(ix,iy))**2
    end do
  end do


  return
end subroutine wfn_rho
