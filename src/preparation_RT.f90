!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation_RT
  use global_variables
  implicit none

  allocate(zwfn(0:Nx,0:NR))

!Initial condition
  zwfn = wfn ! GS wavefunction

  return
end subroutine preparation_RT
