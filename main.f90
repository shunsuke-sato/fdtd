module global_variables
  implicit none

! Mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)

! Physical constants
  real(8),parameter :: c_light = 137.035999139d0
  real(8),parameter :: a_B = 0.52917721067d0
  real(8),parameter :: eV = 27.2116d0

! Grit parameters
  integer,parameter :: NL = -2000, NR = 2000, NM = 300
  integer,parameter :: Nt = 30000
  real(8),parameter :: Hx = 10d0/a_B, dt = 0.02d0
  real(8) :: Lx(NL:NR)

! Fields and materials
  real(8) :: Azc(NL:NR),Azc_new(NL:NR),Azc_old(NL:NR)
  real(8) :: chi(NM),sigma(NM)
  real(8) :: tpulse = 10d0/0.02418d0 ,omega = 1.55d0/eV, I0 = 1d12
  real(8) :: f0 = 5.338d-9*sqrt(I0)
end module global_variables
!===============================================================================
program main
  use global_variables
  implicit none

  call initialization

end program main
!===============================================================================
subroutine initialize
  use global_variables
  implicit none
  integer :: ix

! Grid initialization
  do ix = NL,NR
     Lx(ix) = Hx*(dble(ix)-0.5d0)
  end do

! Material initialization
  chi = 4d0
  sigma = 0d0

  
! Field initialization
  Azc = 0d0; Azc_old = 0d0
  do ix = NL,NR

     tt = -Lx(ix)/c_light - 5d0*Hx/c_light - 0.5d0*tpulse 
     if(abs(tt) < 0.5d0*tpulse) then
        Azc(ix) = - f0/omega*cos(omega*tt)*cos(pi*tt/tpulse)**2
     end if
     tt = -Lx(ix)/c_light - 5d0*Hx/c_light - 0.5d0*tpulse + dt
     if(abs(tt) < 0.5d0*tpulse) then
        Azc_old(ix) = - f0/omega*cos(omega*tt)*cos(pi*tt/tpulse)**2
     end if
  end do
  


  return
end subroutine initialize
!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================
