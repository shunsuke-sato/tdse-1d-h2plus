!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine hpsi
  use global_variables
  implicit none
! finite difference
  real(dp),parameter :: cN0=-205d0/72d0,cN1=8d0/5d0
  real(dp),parameter :: cN2=-1d0/5d0,cN3=8d0/315d0
  real(dp),parameter :: cN4=-1d0/560d0    
  integer :: ix,iy
  real(dp) :: c0,c1,c2,c3,c4
  real(dp) :: c0e,c1e,c2e,c3e,c4e,c0n,c1n,c2n,c3n,c4n
! nine-points formula  
  c0e=-0.5d0*cN0/(dx**2)/mu_e
  c1e=-0.5d0*cN1/(dx**2)/mu_e
  c2e=-0.5d0*cN2/(dx**2)/mu_e
  c3e=-0.5d0*cN3/(dx**2)/mu_e
  c4e=-0.5d0*cN4/(dx**2)/mu_e

  c0n=-0.5d0*cN0/(dr**2)/mu_n
  c1n=-0.5d0*cN1/(dr**2)/mu_n
  c2n=-0.5d0*cN2/(dr**2)/mu_n
  c3n=-0.5d0*cN3/(dr**2)/mu_n
  c4n=-0.5d0*cN4/(dr**2)/mu_n

  tmp_wfn_b(:,:)=0d0
  tmp_wfn_b(0:Nx,0:NR) = tmp_wfn(0:Nx,0:NR)
  tmp_wfn_b(0:Nx,-1) = tmp_wfn(0:Nx,1)
  tmp_wfn_b(0:Nx,-2) = tmp_wfn(0:Nx,2)
  tmp_wfn_b(0:Nx,-3) = tmp_wfn(0:Nx,3)
  tmp_wfn_b(0:Nx,-4) = tmp_wfn(0:Nx,4)

  do iy=0,NR
  do ix=0,Nx
     tmp_hwfn(ix,iy)=(c0e+c0n)*tmp_wfn_b(ix,iy) &
          + c1e*(tmp_wfn_b(ix+1,iy) + tmp_wfn_b(ix-1,iy)) &
          + c2e*(tmp_wfn_b(ix+2,iy) + tmp_wfn_b(ix-2,iy)) &
          + c3e*(tmp_wfn_b(ix+3,iy) + tmp_wfn_b(ix-3,iy)) &
          + c4e*(tmp_wfn_b(ix+4,iy) + tmp_wfn_b(ix-4,iy)) &
          + c1n*(tmp_wfn_b(ix,iy+1) + tmp_wfn_b(ix,iy-1)) &
          + c2n*(tmp_wfn_b(ix,iy+2) + tmp_wfn_b(ix,iy-2)) &
          + c3n*(tmp_wfn_b(ix,iy+3) + tmp_wfn_b(ix,iy-3)) &
          + c4n*(tmp_wfn_b(ix,iy+4) + tmp_wfn_b(ix,iy-4)) 

  end do
  end do

  tmp_hwfn(:,:) = tmp_hwfn(:,:) + v_all(:,:)*tmp_wfn(:,:)
  return
end subroutine hpsi
