! The simulation uses a natural system of units, with the atomic diameter, 
! the atomic mass, the depth of the Lennard-Jones potential, and Boltzmann's 
! constant all set equal to 1. For argon (for example), the unit of distance 
! is 3.4 angstroms, the unit of mass is 40 atomic mass units, and the unit 
! of energy is 0.01 electron-volts; the corresponding unit of time is then 
! 2.2 picoseconds, the unit of velocity is 160 meters per second, and the unit 
! of temperature is 120 kelvin. The Molecular size scrollbar determines the
! scale of the image, in screen pixels per unit of distance.

program argon
  use plot
  use randomseed
  implicit none
  
  ! variables that can be adjusted troughout simulations:
  real(8), parameter :: T_IC = 1.000_8                  ! the temperature of the system in natural units
  real(8), parameter :: rho = .88_8                     ! density of Argon in natural units
  real(8), parameter :: r_verlet = 4.5_8                ! cut off radius for outer sphere
  real(8), parameter :: r_cut =  2.5_8                  ! cut off radius for inner sphere
  real(8), parameter :: dt = 0.004                      ! time step size
  integer, parameter :: N_cell = 4                      ! number of unitcells in one direction    
  integer, parameter :: N_bin = 5000                    ! number of bins for pair correlation function
  integer, parameter :: time_end = 3000                 ! final time instant
  integer, parameter :: hist_time = 100                 ! timesteps between calculating the histogram
  integer, parameter :: eq_time = 500                   ! timesteps untill equilibrium  
  
  ! other declarations
  integer :: time, hist_runs = 0, mu_runs = 0           ! counts the number of times something is added to the histogram and no. enrties for chem pot.
  integer, parameter :: N = 4*N_cell**3                 ! the total number of particles in the system
  integer, parameter :: pairs_max = N*(N-1)/2           ! the maximum amount of entries of interacting particles
  real(8), parameter :: L_box = (N/rho)**(1.0_8/3.0_8)  ! length of the whole box
  real(8), parameter :: L_cell = L_box/N_cell           ! length of a unitcell
  real(8), parameter :: pi = 4*atan(1.0_8)              ! nice way to define pi
  real(8), dimension(3,N) :: pos, vel, acc, tot_dis=0   ! arrays with positions, velocities and accelerations
  real(8), dimension(time_end) :: T, U, kin_energy, mu=0, P=0,disp_sq   ! temperature, pressure, energy, kinetic energy and chemical potential mu in time 
  real(8) :: histogram(N_bin) = 0                       ! define the histogram with number of particles
  real(8) :: T_scaling_factor                           ! scaling factor for the temperature
  real(8) :: msq_vel                                    ! mean square velocity
  integer :: pairs(2,pairs_max)                         ! a array with the neirest pairs
  
  call initialize_system
  call plot_init(L_box)
  do time = 1, time_end
    if (modulo(time,50) == 0) call calculate_pairs      ! every 50 interations recalculate the nearest neighbour list
!    call plot_points(pos(:,:)) 
    call update_pos_vel_acc
    if ((mod(time,50) == 0) .and. time > eq_time) call widom(5000)  ! calculate the chemical potential
  end do
  call plot_close()
  call print_data
 contains

!----------------------------------------------------------------------------!
subroutine initialize_system
  time = 1
  call init_random_seed()                   ! random seed module
  call IC_pos                               ! positions
  call IC_vel                               ! IC  velocity
  call calculate_pairs                      ! calculate neighbour list
  call update_T                             ! rescales the temperate after system is initialized
  call update_acc                           ! calculates the accelerations (Force), LJ potential and the pressure
end subroutine
!----------------------------------------------------------------------------!

!----Put the atoms in a fcc lattice------------------------------------------!
subroutine IC_pos !sets the initial positions in N "boxes" of 4  
  integer :: x, y, z, pc 
  real(8) :: fcc(3,4), origin(3)
  fcc(:,1) = [0.0_8, 0.0_8, 0.0_8]
  fcc(:,2) = [0.0_8, 0.5_8, 0.5_8]
  fcc(:,3) = [0.5_8, 0.0_8, 0.5_8]
  fcc(:,4) = [0.5_8, 0.5_8, 0.0_8]
  fcc = fcc*L_cell
  pc = 0
  do x = 1, N_cell
    do y = 1, N_cell
      do z = 1, N_cell
        origin = [x-1._8, y-1._8, z-1._8]*L_cell
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,1)
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,2)
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,3)
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,4)
      end do
    end do
  end do
!   open (unit=45,file="fcc.dat")
!   write (45,"(3F15.5)") pos
end subroutine IC_pos
!----------------------------------------------------------------------------!

!----Box-Muller Gaussian random number generator-----------------------------!
subroutine random_gauss_routine(random_gauss)
  real(8) :: v1, v2
  real(8), intent(out) :: random_gauss
  call RANDOM_NUMBER(v1)
  call RANDOM_NUMBER(v2)
  random_gauss = sqrt(-2.0*log(v1))*sin(2.0_8*pi*v2)
end subroutine
!----------------------------------------------------------------------------!

!----Velocity initialization: Maxwell distribution---------------------------!
subroutine IC_vel
  real(8) :: avg_vel(3)
  integer :: i
  avg_vel = 0
  do i = 1,N        
    call random_gauss_routine(vel(1,i))
    call random_gauss_routine(vel(2,i))
    call random_gauss_routine(vel(3,i))
    avg_vel(:) = avg_vel(:) + vel(:,i)
  end do 
  avg_vel = avg_vel/N
  !subracting the average velocity so the center of mass does not move
  do i = 1, N
    vel(:,i) = vel(:,i) - avg_vel(:)
  end do  
end subroutine
!----------------------------------------------------------------------------!

!----The atoms that are closer to the atom i then 6 are put in pairs---------!
subroutine calculate_pairs
  integer :: i, j, counter
  real(8) :: r_sq, r(3)
  pairs = 0
  counter = 1
  do i = 1, N-1
    do j = i+1,N
      r = pos(:,i) - pos(:,j)
      r = r - Nint(r/L_box) * L_box
      r_sq = dot_product(r,r)
      if (r_sq < r_verlet**2) then
        pairs(1,counter) = i
        pairs(2,counter) = j
        counter = counter + 1
      end if
    end do
  end do
end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine update_acc
  integer :: pc, i, j
  real(8) :: r_sq, F(3), r(3)
  acc = 0
  U(time) = 0
  kin_energy(time) = 0.5*msq_vel
  do pc = 1, pairs_max
    i = pairs(1,pc)
    j = pairs(2,pc)
    if (i /= 0 .and. j /= 0) then ! only run when pair list entry .ne. 0
      r = pos(:,i) - pos(:,j)
      r = r - Nint(r/L_box) * L_box      
      r_sq = dot_product(r,r)
      if (r_sq < r_cut**2 ) then
        F = 4*(12/r_sq**7-6/r_sq**4)*r
        acc(:,i) = acc(:,i) + F
        acc(:,j) = acc(:,j) - F
        P(time) = P(time) + dot_product(F,r)
      end if
      U(time) = U(time) + 4*(1/r_sq**6-1/r_sq**3)
      call make_histogram(N_bin, sqrt(r_sq), 0._8, r_verlet, pc)
    end if
  end do
end subroutine
!----------------------------------------------------------------------------!

!----Verlet algorithm--------------------------------------------------------!
subroutine update_pos_vel_acc
   vel = vel + 0.5*acc*dt
   pos = modulo(pos + vel*dt, L_box)
   if (time >= eq_time) then 
    tot_dis = tot_dis + vel*dt
    disp_sq(time) = sum(tot_dis*tot_dis)/N
   end if
   call update_T
   call update_acc
   vel = vel + 0.5*acc*dt
end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine update_T
  msq_vel = sum(vel*vel)
  T(time) = msq_vel/(3._8*N)
  T_scaling_factor = sqrt(T_IC/T(time))
  if ((mod(time,20) == 0) .and. time < eq_time-100 ) vel = vel * T_scaling_factor ! scale the temperature every 20 step
end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!  
subroutine make_histogram(N_bin, x, xmin, xmax, pc)
  real(8) :: bin_size, x, xmin, xmax
  integer :: histogram_entry, N_bin, pc
  bin_size = (xmax - xmin) / N_bin
  if (modulo(time,hist_time) == 0 .and. time > eq_time) then
    histogram_entry = ceiling(x/bin_size) ! the bin in which x is put
    if (histogram_entry < N_bin) histogram(histogram_entry) = histogram(histogram_entry) + 1 ! Make sure it doensnt go outside the range of the histogram
    if (pc == 1) hist_runs = hist_runs + 1 ! if pc in subroutine update_acc is 1 -> count
  end if
end subroutine
!----------------------------------------------------------------------------!

!-------The Widom Insertion Method to calculate the chemical potential-------!
subroutine widom(Nghost)
  integer :: i, j, Nghost, k
  real(8) :: postest(3), r(3), r_sq, LJ(Nghost), rnd
  LJ = 0
  do i = 1, Nghost
    if(i==1) mu_runs = mu_runs + 1
    ! add test particle
    do j = 1, 3
      call RANDOM_NUMBER(rnd)
      postest(j) = L_box*rnd
    end do
    ! calculate the distances to other particles
    do k = 1, N
      r = postest(:) - pos(:,k)
      r = r - Nint(r/L_box) * L_box
      r_sq = dot_product(r,r)
      if (r_sq > 1) LJ(i) = LJ(i) + 4*(1/r_sq**6-1/r_sq**3)
    end do
  end do
  mu(mu_runs) = -T(time)*exp(avg(LJ)/T(time))
end subroutine  
!----------------------------------------------------------------------------!  

!----prints data to files----------------------------------------------------!
subroutine print_data
    real(8) :: array(5, size(T(eq_time:time_end)))
    
    P = 1._8 + avg(P(eq_time:time_end))/(3*N*T) - 16*PI*rho/(3*T*r_cut**3)
    histogram = histogram/hist_runs
    
    array(1,:) = T(eq_time:time_end)
    array(2,:) = kin_energy(eq_time:time_end)
    array(3,:) = U(eq_time:time_end)
    array(4,:) = P(eq_time:time_end)
    array(5,:) = disp_sq(eq_time:time_end)
    
    open (unit=10,file="array.dat")
    write (10,"(5F15.5)") array
    open (unit=11,file="mu.dat")
    write (11,"(1E15.5)") mu(1:mu_runs)
    open (unit=12,file="histogram.dat")
    write (12,"(1F15.5)") histogram
    open (unit=18,file="data.dat")
    write (18,*) rho, T_IC, eq_time, time_end, dt, r_cut, r_verlet, N_bin, N
    
    print *, "avg PkT:        ", avg(P(eq_time:time_end))
    print *, "avg T:        ", avg(T(eq_time:time_end))
    print *, "avg Kin:      ", avg(kin_energy(eq_time:time_end))
    print *, "avg mu:       ", avg(mu(1:mu_runs))  
end subroutine
!----------------------------------------------------------------------------!

real(8) function avg(x) result(average)
  real(8), intent(in) :: x(:)
  average = sum(x) / size(x)
end function avg

real(8) function var(x) result(variance)
  real(8), intent(in) :: x(:)
  variance = sum(x*x) / size(x) - avg(x)**2
end function var
end program argon
