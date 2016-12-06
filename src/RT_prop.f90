!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine RT_prop
  use global_variables
  implicit none
  real(dp) :: dipole_m,norm_m,qd_pole_m
  integer :: it

  call dipole(dipole_m, norm_m)
  dipole_t(0) = dipole_m; norm_t(0) = norm_m
  write(*,"(A,2x,I4,2e16.6e3)")"it=",0,dipole_m,norm_m
  do it = 1,Nt_iter
    call dt_evolve(it)
    call multi_pole(dipole_m, qd_pole_m, norm_m)
    dipole_t(it) = dipole_m; quadrupole_t(it) = qd_pole_m; norm_t(it) = norm_m
    write(*,"(A,2x,I7,2e16.6e3)")"it=",it,it*dt,norm_m
  end do

  call write_td_results

  return
end subroutine RT_prop
