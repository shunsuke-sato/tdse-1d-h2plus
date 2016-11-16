!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine RT_prop
  use global_variables
  implicit none
  integer :: it

  do it = 1,Nt_iter
    call dt_eveolve(it)
  end do

  return
end subroutine RT_prop
