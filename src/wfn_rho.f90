!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine wfn_rho
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp

  rho(:)=0d0

!  do ix = 0,Nx
!    tmp = 0d0
!    do iy = 0,Nx
!      tmp = tmp + wfn(ix,iy)**2
!    end do
!    rho(ix) = tmp
!  end do

  do ix = 0,Nx
    rho(ix) = sum(wfn(ix,:)**2)
  end do
  rho = rho*dr

  do ix = 0,NR
    rho_n(ix) = sum(wfn(:,ix)**2)
  end do
  rho_n = rho_n*dx

  return
end subroutine wfn_rho
