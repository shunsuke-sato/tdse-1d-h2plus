!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer,parameter :: Nexp=4
  integer :: iexp 
  complex(zp) :: zfact

! Propagator = Taylor expantions
  zfact = 1d0
  ztmp_wfn=zwfn
  do iexp = 1,Nexp
    zfact = zfact*(-zI*dt)/dble(iexp)
    call zhpsi
    zwfn = zwfn + zfact*ztmp_hwfn
    ztmp_wfn = ztmp_hwfn
  end do

end subroutine dt_evolve
