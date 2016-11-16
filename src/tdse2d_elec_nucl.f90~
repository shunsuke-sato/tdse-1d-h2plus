module global_variables


! parameter
  integer,parameter :: dp = kind(0d0), zp = kind((1d0,1d0))
  complex(zp),parameter :: zI = (0d0,1d0)
  real(dp),parameter :: pi=3.14159265359d0
  
  real(8),parameter :: mass_H = 1836.15267389d0
  real(8),parameter :: mu_e = 2d0*mass_H/(2d0*mass_H+1d0)
  real(8),parameter :: mu_n = 0.5d0*mass_H


! mesh
  integer :: Nx,NR
  real(dp) :: length_x,dx,length_r,dr
  real(dp),allocatable :: xn(:),Rn(:),xRn(:,:)

! GS: wave-function, density, potential and so on
  integer :: Ncg
  real(dp),allocatable :: wfn(:,:), rho(:), v_ext(:), v_int(:,:), v_all(:,:), rho_n(:)

! TD
  character(4) :: RT_mode
  integer :: Nt_iter
  real(dp) :: T_calc,dt,kick_mom
  complex(zp),allocatable :: zwfn(:,:)
  real(dp),allocatable :: dipole_t(:),norm_t(:)
  real(dp) :: field_max,field_duration,field_omega
  real(dp) :: field_max_eV_per_AA,field_duration_fs,field_omega_eV
  real(dp),allocatable :: field_t(:)

! temporary
  real(dp),allocatable :: tmp_wfn(:,:),tmp_hwfn(:,:),tmp_wfn_b(:,:)
  complex(zp),allocatable :: ztmp_wfn(:,:),ztmp_hwfn(:,:),ztmp_wfn_b(:,:)

end module global_variables
!=========================================================================================
program main
use global_variables
  implicit none

  write(*,'(A)')'Start qm1d'
  call input
  call mesh

  call preparation_GS
  call GS_CG

!  stop
  call preparation_RT
  call RT_prop
!
!  call write_results

end program main
!=========================================================================================
subroutine input
  use global_variables
  implicit none

! == input parameter == !                                                                                                                       
  write(*,'(A)')'===== Input parameter ============================================================='
  write(*,'(A)')

  Nx = 400
  length_x = 40d0
  Nr = 400
  length_R = 10d0
! GS
  Ncg = 2400
! TD
  T_calc = 20d0
  dt = 0.01
  Nt_iter = aint(T_calc/dt)+1

  RT_mode='kick' ! kick or gs
  kick_mom = 1e-3

  field_max_eV_per_AA = 1d0
  field_duration_fs = 10d0
  field_omega_eV = 1.55d0

  field_max = field_max_eV_per_AA * (0.529d0/27.2d0)
  field_duration = field_duration_fs/0.02418d0
  field_omega = field_omega_eV/27.2d0

  write(*,'(A,2x,I4)')'Nx =',Nx
  write(*,'(A,2x,e26.16e3)')'length_x =',length_x
  write(*,'(A,2x,e26.16e3)')'T_calc =',T_calc
  write(*,'(A,2x,e26.16e3)')'dt =',dt
  write(*,'(A,2x,I10)')'Nt_iter =',Nt_iter


! temporary array
  allocate(tmp_wfn(0:Nx,0:NR),tmp_hwfn(0:Nx,0:NR),tmp_wfn_b(-4:Nx+4,-4:NR+4))
  allocate(ztmp_wfn(0:Nx,0:NR),ztmp_hwfn(0:Nx,0:NR),ztmp_wfn_b(-4:NR+4,-4:NR+4))

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete Input parameter ===================================================' 
  return
end subroutine input
!=========================================================================================
subroutine mesh
  use global_variables
  implicit none
  integer :: ix, iy 

  write(*,'(A)')'===== Making mesh ================================================================'
  write(*,'(A)')
  write(*,'(A)')

  allocate(xn(0:Nx),Rn(0:NR),xRn(0:Nx,0:NR))
  dx = length_x/dble(Nx)
  dr = length_R/dble(NR)
  
  do ix = 0,Nx
     xn(ix) = dx*dble(ix) - 0.5d0*length_x
  end do

  do ix = 0,NR
    Rn(ix) = dr*dble(ix)
  end do
  
  write(*,'(A)')'===== Complete Making mesh ========================================================'
  return
end subroutine mesh
!=========================================================================================
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
!=========================================================================================
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
!=========================================================================================
subroutine GS_CG
  use global_variables
  implicit none
  real(dp),allocatable :: xvec(:,:),pvec(:,:),rvec(:,:)
  real(dp),allocatable :: hxvec(:,:),gvec(:,:),hpvec(:,:)

  real(dp) :: xx,pp,xp,xhx,php,xhp,esp,esp_res,gg,gg0
  real(dp) :: ss,lambda,alpha,beta,aa,bb,cc
  integer :: iorb,iorb_t,ix,iy,iter_cg
  real(dp) :: esp_iter(Ncg),esp_res_iter(Ncg)

  allocate( xvec(0:Nx,0:NR),pvec(0:Nx,0:NR),rvec(0:Nx,0:NR) )
  allocate( hxvec(0:Nx,0:NR),gvec(0:Nx,0:NR),hpvec(0:Nx,0:NR) )

  write(*,'(A)')'===== Ground state calculation ===================================================='
  write(*,'(A)')
  write(*,'(A)')

  xvec(:,:)=wfn(:,:)

  tmp_wfn = xvec; call hpsi; hxvec = tmp_hwfn
  
  xx=sum(xvec(:,:)**2)*dx*dr
  xhx=sum(xvec(:,:)*hxvec(:,:))*dx*dr
  lambda=xhx/xx
  do iter_cg=1,Ncg
     gvec(:,:)=2d0*(hxvec(:,:)-lambda*xvec(:,:))/xx
     
     gg0=sum(gvec(:,:)**2)*dx*dr
     select case(iter_cg)
     case(1)
        pvec(:,:)=-gvec(:,:)
     case default
        beta=gg0/gg
        pvec(:,:)=-gvec(:,:)+beta*pvec(:,:)
     end select
     gg=gg0

     tmp_wfn = pvec; call hpsi; hpvec = tmp_hwfn
        
     pp=sum(pvec(:,:)**2)*dx*dr
     php=sum(pvec(:,:)*hpvec(:,:))*dx*dr
     xp=sum(xvec(:,:)*pvec(:,:))*dx*dr
     xhp=sum(hxvec(:,:)*pvec(:,:))*dx*dr

     aa=php*xp-xhp*pp
     bb=php*xx-xhx*pp
     cc=xhp*xx-xhx*xp
     ss=bb**2-4d0*aa*cc
     if(ss > 0d0)then
        alpha=(-bb+sqrt(ss))/(2d0*aa)
     else
        exit
     end if
     
     xvec(:,:)=xvec(:,:)+alpha*pvec(:,:)

     tmp_wfn = xvec; call hpsi; hxvec = tmp_hwfn
     xx=sum(xvec(:,:)**2)*dx*dr
     xhx=sum(xvec(:,:)*hxvec(:,:))*dx*dr
     lambda=xhx/xx
     esp_iter(iter_cg)=lambda
     esp_res_iter(iter_cg)=sum((hxvec(:,:)-lambda*xvec(:,:))**2)*dx*dr
     
  end do

  xvec(:,:)=xvec(:,:)/sqrt(xx)
  tmp_wfn = xvec; call hpsi; hxvec = tmp_hwfn
  esp=sum(xvec(:,:)*hxvec(:,:))*dx*dr
  esp_res=sum((hxvec(:,:)-esp*xvec(:,:))**2)*dx*dr
  wfn(:,:)=xvec(:,:)
!  if(wfn(Nx/2,Nx/2) < 0d0) wfn(:,:) = -wfn(:,:)


  write(*,'(A)')'esp,     esp_res'
  write(*,'(e16.6e3,3x,e16.6e3)')esp,esp_res


  open(90,file = 'GS_data.out')
  write(90,'(A)')'# iter,  esp(iter),  esp_res(iter)'
  do iter_cg =1,Ncg
     write(90,'(I6,2x,e26.16e3,2x,e26.16e3)')iter_cg,esp_iter(iter_cg),esp_res_iter(iter_cg)
  end do
  close(90)

  open(90,file = 'GS_wfn.out')
  write(90,'(A)')'# x, y, wfn(x,y)'
  do ix =1,Nx
  do iy =1,NR
     write(90,'(100e26.16e3)')xn(ix),Rn(iy),wfn(ix,iy)
  end do
     write(90,*)
  end do
  close(90)

  call wfn_rho
  open(90,file = 'GS_rho_e.out')
  write(90,'(A)')'# x,  rho(x)'
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),rho(ix)
  end do
  close(90)
  write(*,*)'rho',sum(rho(:))*dx


  open(90,file = 'GS_rho_e_R2.0.out')
  write(90,'(A)')'# x,  rho(x)'
  ss = sum(wfn(:,80)**2)*dx
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),wfn(ix,80)**2/ss
  end do
  close(90)

  open(90,file = 'GS_rho_e_R2.5.out')
  write(90,'(A)')'# x,  rho(x)'
  ss = sum(wfn(:,100)**2)*dx
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),wfn(ix,100)**2/ss
  end do
  close(90)

  open(90,file = 'GS_rho_e_R3.0.out')
  write(90,'(A)')'# x,  rho(x)'
  ss = sum(wfn(:,120)**2)*dx
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),wfn(ix,120)**2/ss
  end do
  close(90)


  open(90,file = 'GS_rho_n.out')
  write(90,'(A)')'# x,  rho(x)'
  do ix =0,NR
     write(90,'(100e26.16e3)')Rn(ix),rho_n(ix)
  end do
  close(90)
  write(*,*)'rho',sum(rho_n(:))*dr


  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== End Ground state calculation ================================================'  

  return
end subroutine GS_CG
!=========================================================================================
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
!=========================================================================================
subroutine preparation_RT
  use global_variables
  implicit none

  allocate(zwfn(0:Nx,0:NR))

!Initial condition
  zwfn = wfn ! GS wavefunction

  return
end subroutine preparation_RT
!=========================================================================================
subroutine RT_prop
  use global_variables
  implicit none
  integer :: it

  do it = 1,Nt_iter
    call dt_eveolve(it)
  end do

  return
end subroutine RT_prop
!=========================================================================================
subroutine dt_eveolve(it)
  use global_variables
  implicit none
  integer :: it

end subroutine dt_eveolve
