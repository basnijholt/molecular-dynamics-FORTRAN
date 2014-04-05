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
  integer, parameter :: N_cell = 6                      ! number of unitcells in one direction  
  real(8), parameter :: T_IC = 1.095_8                  ! the temperature of the system in natural units
  real(8), parameter :: rho = .88_8                     ! density of Argon in natural units
  real(8), parameter :: r_verlet = 4.5_8                ! cut off radius for outer sphere
  real(8), parameter :: r_cut =  2.5_8                  ! cut off radius for inner sphere
  integer, parameter :: N_bin = 2000                    ! number of bins for pair correlation function
  integer, parameter :: time_end = 5000                 ! final time instant
  integer, parameter :: hist_time = 100                 ! timesteps between calculating the histogram
  integer, parameter :: eq_time = 1000                   ! timesteps untill equilibrium  
  
  ! other declarations
  integer :: time, i, hist_runs = 0, mu_runs = 0        ! counts the number of times something is added to the histogram and no. enrties for chem pot.
  integer, parameter :: N = 4*N_cell**3                 ! the total number of particles in the system
  integer, parameter :: pairs_max = N*(N-1)/2           ! the maximum amount of entries of interacting particles
  real(8), parameter :: L_box = (N/rho)**(1.0_8/3.0_8)  ! length of the whole box
  real(8), parameter :: L_cell = L_box/N_cell           ! length of a unitcell
  real(8), parameter :: pi = 4*atan(1.0_8)              ! nice way to define pi
  real(8), dimension(3,N) :: pos, vel, acc, tot_dis=0   ! arrays with positions, velocities and accelerations
  real(8), dimension(time_end) :: T, U, kin_energy, mu=0, PkT=0,disp_sq   ! temperature, pressure, energy, kinetic energy and chemical potential mu in time 
  real(8) :: C(2, time_end-eq_time)                     ! the auto correlation function array
  real(8) :: histogram(N_bin+1) = 0                     ! define the histogram with number of particles
  real(8) :: pair_correlation(2,N_bin)                  ! the pair correlation function vs distance
  real(8) :: T_scaling_factor                           ! scaling factor for the temperature
  real(8) :: msq_vel                                    ! mean square velocity
  integer :: pairs(2,pairs_max)                         ! a array with the neirest pairs
  call initialize_system
!   call plot_init(L_box)
  do time = 1, time_end
    if (modulo(time,50) == 0) call calculate_pairs                  ! every 50 interations recalculate the nearest neighbour list
!     call plot_points(pos(:,:)) 
    call update_pos_vel_acc
    if ((mod(time,50) == 0) .and. time > eq_time) call widom(5000)  ! calculate the chemical potential
  end do
  call plot_close()
  call make_pcf ! finnish the pair correlation function
  call auto_correlation(U(eq_time:time_end)-avg(U(eq_time:time_end))) ! subtracts the average
  call specific_heat
  call print_data(1) ! print T, mu, pcf, U, kin_energy, P and tot_dis to file
 contains

!----------------------------------------------------------------------------!
subroutine initialize_system
  time = 1
  call init_random_seed()                   ! random seed module
  call IC_pos                               ! positions
  call IC_vel                               ! IC  velocity
  call calculate_pairs                      ! calculate neighbour list
  call update_T
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
end subroutine IC_pos
!----------------------------------------------------------------------------!

!----Generate a random number------------------------------------------------!
real(8) function rand_number() result(q)
  call RANDOM_NUMBER(q)
end function rand_number
!----------------------------------------------------------------------------!

!----Box-Muller Gaussian random number generator-----------------------------!
subroutine random_gauss_routine(random_gauss)
  real(8) :: v1, v2
  real(8), intent(out) :: random_gauss
  v1 = rand_number()
  v2 = rand_number()
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
  integer :: pc, i, j !pc for pair counter, i and j index of paticle
  real(8) :: r_sq, F(3), r(3)
  acc = 0
  U(time) = 0
  kin_energy(time) = 0.5*msq_vel
  do pc = 1, pairs_max
    i = pairs(1,pc)
    j = pairs(2,pc)
    if (i /= 0 .and. j /= 0) then ! do not use the terms with entries which are 0 in the pairs
      r = pos(:,i) - pos(:,j)
      r = r - Nint(r/L_box) * L_box      
      r_sq = dot_product(r,r)
      ! calculate the F, LJ and the pressure
      if (r_sq < r_cut**2 ) then
        ! Newtons thrird law, the force of i on j, and the other way around.
        F = 4*(12/r_sq**7-6/r_sq**4)*r
        acc(:,i) = acc(:,i) + F
        acc(:,j) = acc(:,j) - F
        PkT(time) = PkT(time) + dot_product(F,r)
      end if
      U(time) = U(time) + 4*(1/r_sq**6-1/r_sq**3)
      call make_histogram(N_bin, sqrt(r_sq), 0._8, r_verlet, pc)
    end if
  end do
end subroutine
!----------------------------------------------------------------------------!

!----Verlet algorithm--------------------------------------------------------!
subroutine update_pos_vel_acc
  real(8), parameter :: dt = 0.004
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
subroutine make_histogram(N_bin, x, xmin, xmax, counter)
  real(8) :: bin_size, x, xmin, xmax
  integer :: histogram_entry, N_bin, counter
  bin_size = (xmax - xmin) / N_bin
  ! only run the histogram calculation every hist_time steps and only after the system is in equilibrium.
  if (modulo(time,hist_time) == 0 .and. time > eq_time) then
    histogram_entry = ceiling(x/bin_size) ! the bin in which x is put
    if (histogram_entry < N_bin) histogram(histogram_entry) = histogram(histogram_entry) + 1 ! Make sure it doensnt go outside the range of the histogram
    if (counter == 1) hist_runs = hist_runs + 1 ! count the number of times that the histogram is summed over
  end if
end subroutine
!----------------------------------------------------------------------------!

!----Part of the histogram is calculated in update_acc the rest here----!
subroutine make_pcf
  real(8), parameter :: bin_size = r_verlet / N_bin
  ! output the pair correlation function
  histogram = histogram/hist_runs ! to take the average over the histograms
  do i = 1, N_bin
    pair_correlation(1,i) = i*bin_size ! distance (r) in the first column
    pair_correlation(2,i) = 2*histogram(i)/(4*pi*(i*bin_size)**2*bin_size*N*rho)    ! normalize the histogram and change it into a pair correlation function
  end do
end subroutine
!----------------------------------------------------------------------------!

!-------The Widom Insertion Method to calculate the chemical potential-------!
subroutine widom(Nghost)
  integer :: i, j, Nghost, k
  real(8) :: postest(3), r(3), r_sq, LJ(Nghost)
  LJ = 0
  do i = 1, Nghost
    if(i==1) mu_runs = mu_runs + 1
    ! add test particle
    do j = 1, 3
      postest(j) = L_box*rand_number()
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

!----------------------------------------------------------------------------!
subroutine auto_correlation(A)
  integer :: time_delay, time_max
  real(8) :: A(:)
  C = 0
  do time_delay = 1, size(A)
    C(1,time_delay) = time_delay
    time_max = time_end - time_delay
    do time = eq_time, time_max
      C(2,time_delay) = C(2,time_delay) + A(time)*A(time + time_delay)/time_max
    end do
  end do
  call print_data(2) ! print auto correlation to file
end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine specific_heat
  real(8) :: Cv, x(eq_time:time_end)
  x = kin_energy(eq_time:time_end)
  Cv = 1/(2._8/(3._8*N) - Var(x)/avg(x)**2 )
  ! in later measurements the Cv is calculated in the python file, which incorperates the error!
  print *, "T and density:", T_IC, rho
  print *,"calculated Cv: ", Cv, "ideal Cv: high T low rho", 3.0_8/2.0_8*N, "ideal Cv: low T high rho", 3.0_8*N
  open (unit=26,file="cv.dat")
  write (26,*) cv
end subroutine
!----------------------------------------------------------------------------!

!----prints data to files----------------------------------------------------!
subroutine print_data(switch)
integer :: switch
real(8) :: diff_slope
  if (switch == 1) then
    PkT = 1._8 + avg(PkT(eq_time:time_end))/(3*N*T) - 16*PI*rho/(3*T*r_cut**3)
    ! print the temperature and pressure to file
    open (unit=10,file="t.dat")
    write (10,"(1F15.5)") T(eq_time:time_end)
    ! print the chemical potential
    open (unit=11,file="mu.dat")
    write (11,"(1E15.5)") mu(1:mu_runs)
    ! print the pair correlation function to file
    open (unit=12,file="pair_correlation.dat")
    write (12,"(2F15.5)") pair_correlation
    ! print the list with potential plus kinetic energy over time
    open (unit=14,file="u.dat")
    write (14,"(1F15.5)") U(eq_time:time_end)!+kin_energy(eq_time:time_end)
    ! print the list with kinetic energy over time
    open (unit=15,file="kin_energy.dat")
    write (15,"(1F15.5)") kin_energy(eq_time:time_end)
    ! print pressure to file
    open (unit=16,file="p.dat")
    write (16,"(1F15.5)") PkT(eq_time:time_end)
    diff_slope = (sum(tot_dis*tot_dis)/N) / size(disp_sq(eq_time:time_end))
    
!     do i = 1, size(disp_sq(eq_time:time_end))
!       ideal_diff(i) = i*diff_slope
!     end do  
!     open (unit=20,file="ideal_diff.dat")
!     write (20,"(1F15.5)") ideal_diff
    
    ! print the square displacement to file
    open (unit=17,file="disp_sq.dat")
    write (17,"(1E15.10)") disp_sq(eq_time:time_end)!-ideal_diff
    ! diffusion slope
    open (unit=19,file="diff_slope.dat")
    write (19,"(1E15.10)") diff_slope
    ! print density
    open (unit=18,file="rho.dat")
    write (18,*) rho
    open (unit=21,file="t0.dat")
    write (21,*) T_IC   
    open (unit=22,file="time_end.dat")
    write (22,*) time_end
    open (unit=27,file="eq_time.dat")
    write (27,*) eq_time

    print *, "avg PkT:        ", avg(PkT(eq_time:time_end))
    print *, "Var PkT:        ", var(PkT(eq_time:time_end))
    print *, "       "
    print *, "avg T:        ", avg(T(eq_time:time_end))
    print *, "Var T:        ", var(T(eq_time:time_end))
    print *, "       "
    print *, "avg Kin:      ", avg(kin_energy(eq_time:time_end))
    print *, "Var Kin:      ", var(kin_energy(eq_time:time_end))
    print *, "       "
    print *, "avg mu:       ", avg(mu(1:mu_runs))
    print *, "Var mu:       ", var(mu(1:mu_runs))
    print *, "       "
    print *, "diffusion D:  ", diff_slope / 6 ! divide by dt aswell!
  end if
  
  if (switch == 2) then
    ! print the auto correlation to file
    open (unit=13,file="auto_correlation.dat")
    write (13,"(2E20.8)") C
  end if
  
  
  
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
