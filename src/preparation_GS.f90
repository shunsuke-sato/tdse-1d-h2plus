!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation_GS
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp,r,rp,rm

  write(*,'(A)')'===== preparatin_GS =============================================================='
  write(*,'(A)')
  write(*,'(A)')

  allocate(wfn(0:Nx,0:NR), v_ext(0:Nx), v_int(0:Nx,0:NR), v_all(0:Nx,0:NR))
  allocate(rho(0:Nx),rho_n(0:NR))

  write(*,'(A)')'=== preparing initial wave-function ===='
! wfn(0,:) == wfn(Nx,:) == 0 
  wfn(:,:) = 0d0

  do ix = 1,Nx-1  
  do iy = 1,NR-1
     call random_number(tmp)
     tmp = tmp-0.5d0
     wfn(ix,iy) = tmp
  end do
  end do

  tmp_wfn(:,:) = wfn(:,:)

! Normalize
  tmp = sum(wfn(:,:)**2)*dx*dr
  wfn(:,:)=wfn(:,:)/sqrt(tmp)


  write(*,'(A)')'=== preparing external potential ===='
  v_ext = 0d0
!  do ix=0,Nx
!!     v_ext(ix) = 0.5d0*xn(ix)**2
!     v_ext(ix) = -2d0/sqrt(1d0+xn(ix)**2)
!  end do

  write(*,'(A)')'=== preparing interaction potential ===='
  do iy=0,NR
  do ix=0,Nx
     rp = xn(ix) + 0.5d0*Rn(iy)
     rm = xn(ix) - 0.5d0*Rn(iy)
     r = Rn(iy)
!     v_int(ix,iy) = 0d0 !1d0/sqrt(1d0+r**2)
     v_int(ix,iy) = -1d0/sqrt(1d0+rp**2)-1d0/sqrt(1d0+rm**2)+1d0/sqrt(0.03d0+r**2)
  end do
  end do

  write(*,'(A)')'=== preparing total potential ===='
  do iy=0,NR
  do ix=0,Nx
!     v_all(ix,iy) = v_ext(ix) + v_ext(iy) + v_int(ix,iy)
     v_all(ix,iy) = v_int(ix,iy)
  end do
  end do
  
  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete preparatin_GS ====================================================='

  return
end subroutine preparation_GS
