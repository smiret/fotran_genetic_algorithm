!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!module genetic_module

!  implicit none

!  contains
  


!end module genetic_module



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Isotropic_Mech_GA

  !use genetic_module


  !Declarations
  implicit none
  
  double precision :: k1, u1, l1, rho1, kstardesired, ustardesired, lstardesired, rhostardesired
  double precision :: k_low, k_high, u_low, u_high, l_low, l_high, rho_low, rho_high
  double precision :: w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat, tolk, tolu
  double precision :: k2, u2, l2, v2, rho2, p, q
  double precision :: kstar, ustar, lstar, vstar, rhostar, cost
  double precision :: c1k, c1u, c2k, c2u
  double precision :: del_k1, del_k2, del_u1, del_u2, del_l1, del_l2, del_rho1, del_rho2, del_v2
  double precision, dimension(1000,1) :: k1vec, u1vec, l1vec, v1vec, rho1vec
  double precision, dimension(1000,1) :: k2vec, u2vec, l2vec, v2vec, rho2vec
  double precision, dimension(1000,1) :: kstarvec, ustarvec, lstarvec, vstarvec, rhostarvec
  double precision, dimension(1000,1) :: c1kvec, c1uvec, c2kvec, c2uvec
  double precision :: k1_range, u1_range, l1_range, rho1_range, v1_range, k1_low, k1_high, u1_low, u1_high
  double precision :: l1_low, l1_high, rho1_low, rho1_high, v1_high, v1_low
  double precision :: k2_range, u2_range, l2_range, rho2_range, v2_range, k2_low, k2_high, u2_low, u2_high
  double precision :: l2_low, l2_high, rho2_low, rho2_high, v2_high, v2_low
  double precision :: cost_grad_mag, loop_test, tolerance_cost, tolerance_grad
  double precision, dimension(9,1) :: cost_grad_vector, cost_grad
  double precision, dimension(1000,10) :: Pivec
  double precision, dimension(1000,1) :: rand1, cost_grad_mag_vec, rand2
  double precision, dimension(100,1) :: k1parents, u1parents, l1parents, rho1parents
  double precision, dimension(100,1) :: k2parents, u2parents, l2parents, rho2parents, v2parents
  double precision, dimension(100,1) :: k1children, u1children, l1children, rho1children
  double precision, dimension(100,1) :: k2children, u2children, l2children, rho2children, v2children
  integer :: i, j, generation_limit, gen_count, N, vec_size
  double precision, dimension(10000,1) :: cost_save, grad_cost_save
  character(len=50) :: tag1, tag2, tag3
  character :: file_input
  
  
  generation_limit = 10000
  vec_size = 1000
  !parents to keep & children to make
  N = 100
!!!!!!!!!!!!!!!!!  !Step 0: Define the initial parameters and I/O setup

  
  !I/O stuff

  !Open the input text file to read the material requirements

  

  read *, file_input

  Open(Unit = 32, File = file_input, Status = 'unknown')

  !Start the variable import line by line
  !Tag import - Skip the title line - 
  
  read(32,*)
  read(32,*) tag1, tag2, tag3
  read(32,*) !divider line

  print *, tag2

  !Open the output files

  !Filename1 = imported from txt file
  !Filename2 = imported from txt file

  Open(Unit = 64, File = 'Isotropic_GA'//tag1 , Status = 'unknown')
  Open(Unit = 176, File = 'Isotropic_GA'//tag2 , Status = 'unknown')
  Open(Unit = 200, File = 'Isotropic_GA'//tag3 , Status = 'unknown')

  
  !Material Properties - imported from txt file

  !Skip the title line
  read(32,*)
  read(32,*) kstardesired, ustardesired, lstardesired, rhostardesired
  read(32,*) !divider line
    
  !kstardesired = 98.0d0
  !ustardesired = 44.0d0
  !lstardesired = 6.20d0
  !rhostardesired = 7000.0d0


  !Range values for the material properties - imported from txt file

  !Skip the title line
  read(32,*)
  read(32,*) k1_low, k1_high, u1_low, u1_high, l1_low, l1_high, rho1_low, rho1_high
  read(32,*) !divider line

  !Skip the title line
  read(32,*)
  read(32,*) k2_low,k2_high,u2_low,u2_high,l2_low,l2_high,rho2_low,rho2_high,v2_low,v2_high
  read(32,*) !divider line

  !Cost Function weights - imported from txt file

  !Skip the title line
  read(32,*)
  read(32,*) w1,w2,w3,w4,w5hat,w6hat,w7hat,w8hat,p,q
  read(32,*) !divider line


  !Tolerances - import from txt file

  !Skip the title line
  read(32,*)
  read(32,*) tolk, tolu, tolerance_grad, tolerance_cost
  read(32,*) !divider line
  
  !tolk = 0.50d0
  !tolu = 0.50d0

  !tolerance_grad = 0.010d0
  !tolerance_cost = 0.010d0

  close(32)

  
  
!!!!!!!!!!!!!!!!!  !Step 1: Generate a number of random genetic strings that include the material parameters

  call init_random_seed()
  call random_number(rand1)
  call random_number(rand2)
  

  !Random number vectors for the matrix material
  k1vec = k1_low + (k1_high - k1_low)*rand1
  u1vec = u1_low + (u1_high - u1_low)*rand1
  l1vec = l1_low + (l1_high - l1_low)*rand1
  rho1vec = rho1_low + (rho1_high - rho1_low)*rand1

  k2vec = k2_low + (k2_high - k2_low)*rand2
  u2vec = u2_low + (u2_high - u2_low)*rand2
  l2vec = l2_low + (l2_high - l2_low)*rand2
  v2vec = v2_low + (v2_high - v2_low)*rand2
  rho2vec = rho_low + (rho2_high - rho2_low)*rand2


  k2_range = k2_high - k2_low 
  u2_range = u2_high - u2_low
  l2_range = l2_high - l2_low
  rho2_range = rho2_high - rho2_low
  v2_range = v2_high - v2_low

  

  !Compute the effective properties for each combination
  do i = 1,vec_size,1

     k1 = k1vec(i,1)
     u1 = u1vec(i,1)
     l1 = l1vec(i,1)
     rho1 = rho1vec(i,1)

     k2 = k2vec(i,1)
     u2 = u2vec(i,1)
     l2 = l2vec(i,1)
     rho2 = rho2vec(i,1)
     v2 = v2vec(i,1)



  
  call bounds_number(k1,k2,u1,u2,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)


  kstarvec(i,1) = kstar
  ustarvec(i,1) = ustar
  lstarvec(i,1) = lstar
  rhostarvec(i,1) = rhostar
  c1kvec(i,1) = c1k
  c1uvec(i,1) = c1u
  c2kvec(i,1) = c2k
  c2uvec(i,1) = c2u

end do



!!!!!!!!!!!!!!!!!  !Step 2: Compute the fitness of each genetic string using the cost function

do i = 1,vec_size,1

   kstar = kstarvec(i,1)
   ustar = ustarvec(i,1)
   lstar = lstarvec(i,1)
   rhostar = rhostarvec(i,1)
   c1k = c1kvec(i,1)
   c1u = c1uvec(i,1)
   c2k = c2kvec(i,1)
   c2u = c2uvec(i,1)

   

  call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
     w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
     kstardesired, ustardesired, lstardesired, rhostardesired, &
     cost)


  Pivec(i,1) = cost
  Pivec(i,2) = k1vec(i,1)
  Pivec(i,3) = u1vec(i,1)
  Pivec(i,4) = l1vec(i,1)
  Pivec(i,5) = rho1vec(i,1)
  Pivec(i,6) = k2vec(i,1)
  Pivec(i,7) = u2vec(i,1)
  Pivec(i,8) = l2vec(i,1)
  Pivec(i,9) = rho2vec(i,1)
  Pivec(i,10) = v2vec(i,1)

  
  
end do



!!!!!!!!!!!!!!!!!  !Step 3: Rank the genetic strings by the fitness using a sorting function

call sort1(Pivec)

k1parents(:,1) = Pivec(1:N,2)
u1parents(:,1) = Pivec(1:N,3)
l1parents(:,1) = Pivec(1:N,4)
rho1parents(:,1) = Pivec(1:N,5)

k2parents(:,1) = Pivec(1:N,6)
u2parents(:,1) = Pivec(1:N,7)
l2parents(:,1) = Pivec(1:N,8)
rho2parents(:,1) = Pivec(1:N,9)
v2parents(:,1) = Pivec(1:N,10)





!!!!!!!!!!!!!!!!!  !Step 4: Produce children from the best performing strings using them as bounds and a random number as their distance from each other


do j = 1,N-1,2

   !Matrix - use rand1

   k1children(j,1) = rand1(j,1)*k1parents(j,1) + (1-rand1(j,1))*k1parents(j+1,1)
   k1children(j+1,1) = rand1(j+1,1)*k1parents(j,1) + (1-rand1(j+1,1))*k1parents(j+1,1)

   u1children(j,1) = rand1(j,1)*u1parents(j,1) + (1-rand1(j,1))*u1parents(j+1,1)
   u1children(j+1,1) = rand1(j+1,1)*u1parents(j,1) + (1-rand1(j+1,1))*u1parents(j+1,1)

   l1children(j,1) = rand1(j,1)*l1parents(j,1) + (1-rand1(j,1))*l1parents(j+1,1)
   l1children(j+1,1) = rand1(j+1,1)*l1parents(j,1) + (1-rand1(j+1,1))*l1parents(j+1,1)

   rho1children(j,1) = rand1(j,1)*rho1parents(j,1) + (1-rand1(j,1))*rho1parents(j+1,1)
   rho1children(j+1,1) = rand1(j+1,1)*rho1parents(j,1) + (1-rand1(j+1,1))*rho1parents(j+1,1)


   !!Particles - use rand2

   k2children(j,1) = rand2(j,1)*k2parents(j,1) + (1.0d0-rand2(j,1))*k2parents(j+1,1)
   k2children(j+1,1) = rand2(j+1,1)*k2parents(j,1) + (1.0d0-rand2(j+1,1))*k2parents(j+1,1)

   u2children(j,1) = rand2(j,1)*u2parents(j,1) + (1-rand2(j,1))*u2parents(j+1,1)
   u2children(j+1,1) = rand2(j+1,1)*u2parents(j,1) + (1-rand2(j+1,1))*u2parents(j+1,1)

   l2children(j,1) = rand2(j,1)*l2parents(j,1) + (1-rand2(j,1))*l2parents(j+1,1)
   l2children(j+1,1) = rand2(j+1,1)*l2parents(j,1) + (1-rand2(j+1,1))*l2parents(j+1,1)

   rho2children(j,1) = rand2(j,1)*rho2parents(j,1) + (1-rand2(j,1))*rho2parents(j+1,1)
   rho2children(j+1,1) = rand2(j+1,1)*rho2parents(j,1) + (1-rand2(j+1,1))*rho2parents(j+1,1)

   v2children(j,1) = rand2(j,1)*v2parents(j,1) + (1-rand2(j,1))*v2parents(j+1,1)
   v2children(j+1,1) = rand2(j+1,1)*v2parents(j,1) + (1-rand2(j+1,1))*v2parents(j+1,1)
   

end do



!!!!!!!!!!!!!!!!!  !Step 4.5: Produce mutations of the other strings using random numbers


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIX THE VECTORS HERE


!New Design Variables

!Matrix - rand1

k1vec(1:N,1) = k1parents(1:N,1)
k1vec(N+1:2*N,1) = k1children(1:N,1)
k1vec(2*N+1:vec_size,1) = k1_low + (k1_high - k1_low)*rand1(2*N+1:vec_size,1)

u1vec(1:N,1) = u1parents(1:N,1)
u1vec(N+1:2*N,1) = u1children(1:N,1)
u1vec(2*N+1:vec_size,1) = u1_low + (u1_high - u1_low)*rand1(2*N+1:vec_size,1)

l1vec(1:N,1) = l1parents(1:N,1)
l1vec(N+1:2*N,1) = l1children(1:N,1)
l1vec(2*N+1:vec_size,1) = l1_low + (l1_high - l1_low)*rand1(2*N+1:vec_size,1)

rho1vec(1:N,1) = rho1parents(1:N,1)
rho1vec(N+1:2*N,1) = rho1children(1:N,1)
rho1vec(2*N+1:vec_size,1) = rho1_low + (rho1_high - rho1_low)*rand1(2*N+1:vec_size,1)



!Particles - rand2

k2vec(1:N,1) = k2parents(1:N,1)
k2vec(N+1:2*N,1) = k2children(1:N,1)
k2vec(2*N+1:vec_size,1) = k2_low + (k2_high - k2_low)*rand2(2*N+1:vec_size,1)

u2vec(1:N,1) = u2parents(1:N,1)
u2vec(N+1:2*N,1) = u2children(1:N,1)
u2vec(2*N+1:vec_size,1) = u2_low + (u2_high - u2_low)*rand2(2*N+1:vec_size,1)

l2vec(1:N,1) = l2parents(1:N,1)
l2vec(N+1:2*N,1) = l2children(1:N,1)
l2vec(2*N+1:vec_size,1) = l2_low + (l2_high - l2_low)*rand2(2*N+1:vec_size,1)

rho2vec(1:N,1) = rho2parents(1:N,1)
rho2vec(N+1:2*N,1) = rho2children(1:N,1)
rho2vec(2*N+1:vec_size,1) = rho2_low + (rho2_high - rho2_low)*rand2(2*N+1:vec_size,1)
 
v2vec(1:N,1) = v2parents(1:N,1)
v2vec(N+1:2*N,1) = v2children(1:N,1)
v2vec(2*N+1:vec_size,1) = v2_low + (v2_high - v2_low)*rand2(2*N+1:vec_size,1)




!!!!!!!!!!!!!!!!!  !Step 5: Enforce the design contraints on the new generations of genetic strings

do i = 1,vec_size,1

   !Matrix - material 1

   if (k1vec(i,1) > k1_high) then
      k1vec(i,1) = k1_high

   end if

   if (k1vec(i,1) < k1_low) then
      k1vec(i,1) = k1_low

   end if

   if (u1vec(i,1) > u1_high) then
      u1vec(i,1) = u1_high

   end if

   if (u1vec(i,1) < u1_low) then
      u1vec(i,1) = u1_low

   end if

   if (l1vec(i,1) > l1_high) then
      l1vec(i,1) = l1_high

   end if

   if (l1vec(i,1) < l1_low) then
      l1vec(i,1) = l1_low

   end if

 

   if (rho1vec(i,1) > rho1_high) then
      rho1vec(i,1) = rho1_high

   end if

   if (rho1vec(i,1) < rho1_low) then
      rho1vec(i,1) = rho1_low

   end if

   !Particles - material 2
   
   if (k2vec(i,1) > k2_high) then
      k2vec(i,1) = k2_high

   end if

   if (k2vec(i,1) < k2_low) then
      k2vec(i,1) = k2_low

   end if

   if (u2vec(i,1) > u2_high) then
      u2vec(i,1) = u2_high

   end if

   if (u2vec(i,1) < u2_low) then
      u2vec(i,1) = u2_low

   end if

   if (l2vec(i,1) > l2_high) then
      l2vec(i,1) = l2_high

   end if

   if (l2vec(i,1) < l2_low) then
      l2vec(i,1) = l2_low

   end if

   if (rho2vec(i,1) > rho2_high) then
      rho2vec(i,1) = rho2_high

   end if

   if (rho2vec(i,1) < rho2_low) then
      rho2vec(i,1) = rho2_low

   end if

   if (v2vec(i,1) > v2_high) then
      v2vec(i,1) = v2_high

   end if

   if (v2vec(i,1) < v2_low) then
      v2vec(i,1) = v2_low

   end if

   
   

end do





!!!!!!!!!!!!!!!!!  !Step 6: Repeat the process with the new gene pool and include mutations equal to number of excluded strings

loop_test = 0.0d0

gen_count = 0


do while (loop_test < 1.0d0) !DEBUG HERE!!!!!!!!!!! MEMORY SEGFAULT

   print*, 'Cost minimum'

   print *, minval(Pivec(:,1))
   print *, 'Generation'
   print *, gen_count

   call sort1(Pivec)

  
   call random_number(rand1)
   call random_number(rand2)


   !Compute the effective properties for each combination - NOT FOR PARENTS AND CHILDREN
  do i = 2*N+1,vec_size,1

     k1 = k1vec(i,1)
     u1 = u1vec(i,1)
     l1 = l1vec(i,1)
     rho1 = rho1vec(i,1)

     k2 = k2vec(i,1)
     u2 = u2vec(i,1)
     l2 = l2vec(i,1)
     rho2 = rho2vec(i,1)
     v2 = v2vec(i,1)


  call bounds_number(k1,k2,u1,u2,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)



  kstarvec(i,1) = kstar
  ustarvec(i,1) = ustar
  lstarvec(i,1) = lstar
  rhostarvec(i,1) = rhostar
  c1kvec(i,1) = c1k
  c1uvec(i,1) = c1u
  c2kvec(i,1) = c2k
  c2uvec(i,1) = c2u

end do


   !Calculate the value of the cost function - DONT REDO IT FOR PARENTS AND CHILDREN


   do i = 2*N+1,vec_size,1
            
      
      kstar = kstarvec(i,1)
      ustar = ustarvec(i,1)
      lstar = lstarvec(i,1)
      rhostar = rhostarvec(i,1)
      
      c1k = c1kvec(i,1)
      c1u = c1uvec(i,1)
      c2k = c2kvec(i,1)
      c2u = c2uvec(i,1)

      call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost)
      
      

      Pivec(i,1) = cost
      Pivec(i,2) = k1vec(i,1)
      Pivec(i,3) = u1vec(i,1)
      Pivec(i,4) = l1vec(i,1)
      Pivec(i,5) = rho1vec(i,1)
      Pivec(i,6) = k2vec(i,1)
      Pivec(i,7) = u2vec(i,1)
      Pivec(i,8) = l2vec(i,1)
      Pivec(i,9) = rho2vec(i,1)
      Pivec(i,10) = v2vec(i,1)

      
   end do

   
  
   

   !Calculate the gradient of the cost function and its magnitude

 do i = 1,vec_size,1

     k1 = k1vec(i,1)
     u1 = u1vec(i,1)
     l1 = l1vec(i,1)
     rho1 = rho1vec(i,1)

     
     
     k2 = k2vec(i,1)
     u2 = u2vec(i,1)
     l2 = l2vec(i,1)
     rho2 = rho2vec(i,1)
     v2 = v2vec(i,1)

    

     del_k1 = k1*10**(-2.0d0)
     del_k2 = k2*10**(-2.0d0)
     del_u1 = u1*10**(-2.0d0)
     del_u2 = u2*10**(-2.0d0)
     del_l1 = l1*10**(-2.0d0)
     del_l2 = l2*10**(-2.0d0)
     del_rho1 = rho1*10**(-2.0d0)
     del_rho2 = rho2*10**(-2.0d0)
     del_v2 = v2*10**(-2.0d0)

     !Function creates segmentation fault


   call cost_gradient(k1,k2,u1,u2,l1,l2,rho1,rho2,v2, &
        del_k1, del_k2, del_u1, del_u2, del_l1, del_l2, del_rho1, del_rho2, del_v2, &
        w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
        kstardesired, ustardesired, lstardesired, rhostardesired, &
        cost_grad_vector, cost_grad_mag)

   !!!

   cost_grad_mag_vec(i,1) = cost_grad_mag

   


end do


print *, 'Gradient'
print *, minval(cost_grad_mag_vec)


cost_save(gen_count,1) =  minval(Pivec(:,1))
grad_cost_save(gen_count,1) = minval(cost_grad_mag_vec)


!Check if the conditions are satisfied




   if (minval(Pivec(:,1)) < tolerance_cost) then
      loop_test = 2.0d0
      
   end if



   if (gen_count > generation_limit) then
      loop_test = 2.0d0

   end if


   gen_count = gen_count + 1



   !Update the material properties

  ! print *, 'Line 750'
!!!!!!!!!!!!!!!!! Rank the genetic strings by the fitness using a sorting function
   
   call sort1(Pivec)


   k1parents(:,1) = Pivec(1:N,2)
   u1parents(:,1) = Pivec(1:N,3)
   l1parents(:,1) = Pivec(1:N,4)
   rho1parents(:,1) = Pivec(1:N,5)
   
   k2parents(:,1) = Pivec(1:N,6)
   u2parents(:,1) = Pivec(1:N,7)
   l2parents(:,1) = Pivec(1:N,8)
   rho2parents(:,1) = Pivec(1:N,9)
   v2parents(:,1) = Pivec(1:N,10)


!!!!!!!!!!!!!!!!!  Produce children from the best performing strings using them as bounds and a random number as their distance from each other


   do j = 1,N-1,2


      !Matrix - use rand1
      
      k1children(j,1) = rand1(j,1)*k1parents(j,1) + (1-rand1(j,1))*k1parents(j+1,1)
      k1children(j+1,1) = rand1(j+1,1)*k1parents(j,1) + (1-rand1(j+1,1))*k1parents(j+1,1)

      u1children(j,1) = rand1(j,1)*u1parents(j,1) + (1-rand1(j,1))*u1parents(j+1,1)
      u1children(j+1,1) = rand1(j+1,1)*u1parents(j,1) + (1-rand1(j+1,1))*u1parents(j+1,1)

      l1children(j,1) = rand1(j,1)*l1parents(j,1) + (1-rand1(j,1))*l1parents(j+1,1)
      l1children(j+1,1) = rand1(j+1,1)*l1parents(j,1) + (1-rand1(j+1,1))*l1parents(j+1,1)

      rho1children(j,1) = rand1(j,1)*rho1parents(j,1) + (1-rand1(j,1))*rho1parents(j+1,1)
      rho1children(j+1,1) = rand1(j+1,1)*rho1parents(j,1) + (1-rand1(j+1,1))*rho1parents(j+1,1)

      !!Particles - use rand2

      k2children(j,1) = rand2(j,1)*k2parents(j,1) + (1-rand2(j,1))*k2parents(j+1,1)
      k2children(j+1,1) = rand2(j+1,1)*k2parents(j,1) + (1-rand2(j+1,1))*k2parents(j+1,1)

      u2children(j,1) = rand2(j,1)*u2parents(j,1) + (1-rand2(j,1))*u2parents(j+1,1)
      u2children(j+1,1) = rand2(j+1,1)*u2parents(j,1) + (1-rand2(j+1,1))*u2parents(j+1,1)

      l2children(j,1) = rand2(j,1)*l2parents(j,1) + (1-rand2(j,1))*l2parents(j+1,1)
      l2children(j+1,1) = rand2(j+1,1)*l2parents(j,1) + (1-rand2(j+1,1))*l2parents(j+1,1)

      rho2children(j,1) = rand2(j,1)*rho2parents(j,1) + (1-rand2(j,1))*rho2parents(j+1,1)
      rho2children(j+1,1) = rand2(j+1,1)*rho2parents(j,1) + (1-rand2(j+1,1))*rho2parents(j+1,1)

      v2children(j,1) = rand2(j,1)*v2parents(j,1) + (1-rand2(j,1))*v2parents(j+1,1)
      v2children(j+1,1) = rand2(j+1,1)*v2parents(j,1) + (1-rand2(j+1,1))*v2parents(j+1,1)

   

   end do


!!!!!!!!!!!!!!!!!  Produce mutations of the other strings using random numbers


   !New Design Variables

   !Matrix - rand1

   k1vec(1:N,1) = k1parents(1:N,1)
   k1vec(N+1:2*N,1) = k1children(1:N,1)
   k1vec(2*N+1:vec_size,1) = k1_low + (k1_high - k1_low)*rand1(2*N+1:vec_size,1)

   u1vec(1:N,1) = u1parents(1:N,1)
   u1vec(N+1:2*N,1) = u1children(1:N,1)
   u1vec(2*N+1:vec_size,1) = u1_low + (u1_high - u1_low)*rand1(2*N+1:vec_size,1)

   l1vec(1:N,1) = l1parents(1:N,1)
   l1vec(N+1:2*N,1) = l1children(1:N,1)
   l1vec(2*N+1:vec_size,1) = l1_low + (l1_high - l1_low)*rand1(2*N+1:vec_size,1)

   rho1vec(1:N,1) = rho1parents(1:N,1)
   rho1vec(N+1:2*N,1) = rho1children(1:N,1)
   rho1vec(2*N+1:vec_size,1) = rho1_low + (rho1_high - rho1_low)*rand1(2*N+1:vec_size,1)



   !Particles - rand2

   k2vec(1:N,1) = k2parents(1:N,1)
   k2vec(N+1:2*N,1) = k2children(1:N,1)
   k2vec(2*N+1:vec_size,1) = k2_low + (k2_high - k2_low)*rand2(2*N+1:vec_size,1)

   u2vec(1:N,1) = u2parents(1:N,1)
   u2vec(N+1:2*N,1) = u2children(1:N,1)
   u2vec(2*N+1:vec_size,1) = u2_low + (u2_high - u2_low)*rand2(2*N+1:vec_size,1)

   l2vec(1:N,1) = l2parents(1:N,1)
   l2vec(N+1:2*N,1) = l2children(1:N,1)
   l2vec(2*N+1:vec_size,1) = l2_low + (l2_high - l2_low)*rand2(2*N+1:vec_size,1)

   rho2vec(1:N,1) = rho2parents(1:N,1)
   rho2vec(N+1:2*N,1) = rho2children(1:N,1)
   rho2vec(2*N+1:vec_size,1) = rho2_low + (rho2_high - rho2_low)*rand2(2*N+1:vec_size,1)
 
   v2vec(1:N,1) = v2parents(1:N,1)
   v2vec(N+1:2*N,1) = v2children(1:N,1)
   v2vec(2*N+1:vec_size,1) = v2_low + (v2_high - v2_low)*rand2(2*N+1:vec_size,1)

   
!!!!!!!!!!!!!!!!!  Enforce the design contraints on the new generations of genetic strings

   do i = 1,vec_size,1

      !Matrix - material 1

   if (k1vec(i,1) > k1_high) then
      k1vec(i,1) = k1_high

   end if

   if (k1vec(i,1) < k1_low) then
      k1vec(i,1) = k1_low

   end if

   if (u1vec(i,1) > u1_high) then
      u1vec(i,1) = u1_high

   end if

   if (u1vec(i,1) < u1_low) then
      u1vec(i,1) = u1_low

   end if

   if (l1vec(i,1) > l1_high) then
      l1vec(i,1) = l1_high

   end if

   if (l1vec(i,1) < l1_low) then
      l1vec(i,1) = l1_low

   end if

   if (rho1vec(i,1) > rho1_high) then
      rho1vec(i,1) = rho1_high

   end if

   if (rho1vec(i,1) < rho1_low) then
      rho1vec(i,1) = rho1_low

   end if

   !Particles - material 2
   
   if (k2vec(i,1) > k2_high) then
      k2vec(i,1) = k2_high

   end if

   if (k2vec(i,1) < k2_low) then
      k2vec(i,1) = k2_low

   end if

   if (u2vec(i,1) > u2_high) then
      u2vec(i,1) = u2_high

   end if

   if (u2vec(i,1) < u2_low) then
      u2vec(i,1) = u2_low

   end if

   if (l2vec(i,1) > l2_high) then
      l2vec(i,1) = l2_high

   end if

   if (l2vec(i,1) < l2_low) then
      l2vec(i,1) = l2_low

   end if

   if (rho2vec(i,1) > rho2_high) then
      rho2vec(i,1) = rho2_high

   end if

   if (rho2vec(i,1) < rho2_low) then
      rho2vec(i,1) = rho2_low

   end if

   if (v2vec(i,1) > v2_high) then
      v2vec(i,1) = v2_high

   end if

   if (v2vec(i,1) < v2_low) then
      v2vec(i,1) = v2_low

   end if
      
  
end do



end do



!I/O to write output file
write(64, *) 'Study Details'
write(64, *) !break line

write(64, *) 'Desired kstardesired ustardesired lstardesired rhostardesired', &
     kstardesired, ustardesired, lstardesired, rhostardesired, &
     'Ranges for k1 k1_low k1_high u1_low u1_high l1_low l1_high rho1_low rho1_high', &
     k1_low, k1_high, u1_low, u1_high, l1_low, l1_high, rho1_low, rho1_high, &
     'Ranges for k2 k2_low k2_high u2_low u2_high l2_low l2_high rho2_low rho2_high v2_low v2 high', &
     k2_low,k2_high,u2_low,u2_high,l2_low,l2_high,rho2_low,rho2_high,v2_low,v2_high, &
     'Weights w1 w2 w3 w4 w5hat w6hat w7hat w8hat p q', &
     w1,w2,w3,w4,w5hat,w6hat,w7hat,w8hat,p,q, &
     'Tolerances tolk tolu tolerance_grad tolerance_cost', &
     tolk, tolu, tolerance_grad, tolerance_cost

write(64, *) !break line
     
     



write(64,*)  'Matrix Bulk Modulus', k1vec, 'Matrix Shear Modulus', u1vec, &
       'Matrix Conductivity', l1vec, 'Matrix Density', rho1vec, &
       'Dopant Bulk Modulus', k2vec, 'Dopant Shear Modulus', u2vec, &
       'Dopant Conductivity', l2vec, 'Dopant Density', rho2vec, 'Dopant Volume', v2vec
       !'Cost Minima', cost_save, 'Gradient Minima', grad_cost_save

  !100 FORMAT(A,1X,F16.5, A,1X,F16.5, A,1X,F16.5, A,1X,F16.5, A,1X,F16.5)


write(176, *) 'Study Details'
write(176, *) !break line

write(176, *) 'Desired kstardesired ustardesired lstardesired rhostardesired', &
     kstardesired, ustardesired, lstardesired, rhostardesired, &
     'Ranges for k1 k1_low k1_high u1_low u1_high l1_low l1_high rho1_low rho1_high', &
     k1_low, k1_high, u1_low, u1_high, l1_low, l1_high, rho1_low, rho1_high, &
     'Ranges for k2 k2_low k2_high u2_low u2_high l2_low l2_high rho2_low rho2_high v2_low v2 high', &
     k2_low,k2_high,u2_low,u2_high,l2_low,l2_high,rho2_low,rho2_high,v2_low,v2_high, &
     'Weights w1 w2 w3 w4 w5hat w6hat w7hat w8hat p q', &
     w1,w2,w3,w4,w5hat,w6hat,w7hat,w8hat,p,q, &
     'Tolerances tolk tolu tolerance_grad tolerance_cost', &
     tolk, tolu, tolerance_grad, tolerance_cost

write(176, *) !break line

write(200, *) 'Study Details'
write(200, *) !break line

write(200, *) 'Desired kstardesired ustardesired lstardesired rhostardesired', &
     kstardesired, ustardesired, lstardesired, rhostardesired, &
     'Ranges for k1 k1_low k1_high u1_low u1_high l1_low l1_high rho1_low rho1_high', &
     k1_low, k1_high, u1_low, u1_high, l1_low, l1_high, rho1_low, rho1_high, &
     'Ranges for k2 k2_low k2_high u2_low u2_high l2_low l2_high rho2_low rho2_high v2_low v2 high', &
     k2_low,k2_high,u2_low,u2_high,l2_low,l2_high,rho2_low,rho2_high,v2_low,v2_high, &
     'Weights w1 w2 w3 w4 w5hat w6hat w7hat w8hat p q', &
     w1,w2,w3,w4,w5hat,w6hat,w7hat,w8hat,p,q, &
     'Tolerances tolk tolu tolerance_grad tolerance_cost', &
     tolk, tolu, tolerance_grad, tolerance_cost

write(200, *) !break line
  
  

write(176,*) 'Cost Minima'
write(200,*) 'Gradient Minima'

do i = 1,10000,1

   write(176,*) cost_save(i,1), ',', i, ','
   write(200,*) grad_cost_save(i,1), ',', i, ','

end do



  
  close(64)
  close(176)
  close(200)
 

end program Isotropic_Mech_GA




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Subroutines

subroutine bounds_number(k1,k2,u1,u2,l1,l2,rho1,rho2,v2, &
     kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)

  double precision, intent(in) :: k1, k2, u1, u2, l1, l2, rho1, rho2, v2
  double precision, intent(out) :: kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u
  double precision :: khsl, khsu, uhsl, uhsu, lhsl, lhsu, theta, v1

  !call init_random_seed()
  !call random_number(theta)
  theta = 0.50d0
  !theta = 0.50d0

  !theta = 0.50d0
  v1 = 1.0d0 - v2


  !HS Bounds for k, u and l
  khsl = k1 + (v2/((1.0d0/(k2-k1)) + (3.0d0*(1.0d0-v2))/(3.0d0*k1 + 4.0d0*u1)))
  
  khsu = k2 + ((1.0d0-v2)/((1.0d0/(k1-k2)) + (3.0d0*v2)/(3.0d0*k2 + 4.0d0*u2)))

  uhsl = u1 + (v2/((1.0d0/(u2-u1))+((6.0d0*(1.0d0-v2)*(k1+2.0d0*u1))/	&
       (5.0d0*u1*(3.0d0*k1+4.0d0*u1)))))
  
  uhsu = u2 + ((1.0d0-v2)/((1.0d0/(u1-u2))+((6.0d0*v2*(k2*2.0d0*u2))/ &
       (5.0d0*u2*(3.0d0*k2+4.0d0*u2)))))

  lhsl = l1 + (v2 / ((1.0d0/(l2 - l1)) + ((1 - v2)/l1))) !old formula !((l1*l2) + 2.0d0*l1*(l1*v2 + l2*v1))/(2.0d0*l1 + l1*v1 + l2*v2)

  lhsu = l2 + ((1 - v2) / ((1.0d0/(l1 - l2)) + ((v2)/l2))) !old formula !(l1*l2 + 2.0d0*l2*(l1*v2 + l2*v1))/(2*l2 + l2*v2 + l1*v1)


  
  !Effective Values for k and u
  kstar = theta*khsu + (1.0d0-theta)*khsl
  ustar = theta*uhsu + (1.0d0-theta)*uhsl
  lstar = theta*lhsu + (1.0d0-theta)*lhsl
  rhostar = (rho1*v1 + rho2*v2)

  
  
  !Stress concentration tensors
  c1k = (1.0d0/v2)*(k2/kstar)*((kstar-k1)/(k2-k1))
  c1u = (1.0d0/v2)*(u2/ustar)*((ustar-u1)/(u2-u1))
  v1 = 1.0d0 - v2
  c2k = (1.0d0/v1)*(1.0d0-v2*c1k)
  c2u = (1.0d0/v1)*(1.0d0-v2*c1u)
  
end subroutine bounds_number


subroutine cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
     w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
     kstardesired, ustardesired, lstardesired, rhostardesired, &
     cost)

  double precision, intent(in) :: kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u
  double precision, intent(in) :: kstardesired, ustardesired, lstardesired, rhostardesired
  double precision, intent(in) :: w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu
  double precision, intent(out) :: cost
  double precision :: w5, w6, w7, w8


  if (ABS((c1k-1.0d0)/c1k) <= tolk) then
        w5 = 0.0d0
     else
        w5 = w5hat 
     end if
     
     if (ABS((c1u-1.0d0)/c1u) <= tolu) then
        w6 = 0.0d0
     else
        w6 = w6hat
     end if
     
     if (ABS((c2k-1.0d0)/c2k) <= tolk) then
        w7 = 0.0d0
     else
        w7 = w7hat
     end if

     if (ABS((c2u-1.0d0)/c2u) <= tolu) then
        w8 = 0.0d0
     else
        w8 = w8hat
     end if
     
     cost = w1*abs((kstar/kstardesired)-1.0d0)**p  	&
          + w2*abs((ustar/ustardesired)-1.0d0)**p  	&
          + w3*abs((lstar/lstardesired)-1.0d0)**p  	&
          + w4*(abs(rhostar/rhostardesired - 1.0d0))**p &
          + w5*(abs((((c1k-1.0d0)/c1k)/tolk)-1.0d0))**q &
          + w6*(abs((((c1u-1.0d0)/c1u)/tolu)-1.0d0))**q &
          + w7*(abs((((c2k-1.0d0)/c2k)/tolk)-1.0d0))**q	&
          + w8*(abs((((c2u-1.0d0)/c2u)/tolu)-1.0d0))**q



end subroutine cost_number




   
          




!Watch this sorting algorithm - it's replacing values in the genetic string (maybe)
subroutine sort1 (Pi)

  double precision,intent(inout), dimension(1000,10) :: Pi
  double precision :: a2
  double precision, dimension(10,1) :: a1, a3


     
  do j = 2,1000,1

     !Take the genetic string at point j
     a1(:,1) = Pi(j,:)

     !Take the cost function value of that genetic string
     a2 = a1(1,1)


     !Compare the value of the cost function for the given string to its neighbor
     !Loop backwards from the jth value
     do i = j,2,-1
        
        if (Pi(i,1) < Pi(i-1,1)) then

           a3(:,1) = Pi(i,:)

           Pi(i,:) = Pi(i-1,:)
           Pi(i-1,:) = a3(:,1)

     end if

  end do

end do



 
end subroutine sort1



subroutine cost_gradient(k1,k2,u1,u2,l1,l2,rho1,rho2,v2, &
     del_k1, del_k2, del_u1, del_u2, del_l1, del_l2, del_rho1, del_rho2, del_v2, &
     w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
     kstardesired, ustardesired, lstardesired, rhostardesired, &
     cost_grad_vector, cost_grad_mag)

  !Do everything as floats for now, not vector valued ... keep the vectors in the main program
  double precision, intent(in) :: k1, k2, u1, u2, l1, l2, rho1, rho2, v2
  double precision, intent(in) :: del_k1, del_k2, del_u1, del_u2, del_l1, del_l2
  double precision, intent(in) :: del_rho1, del_rho2, del_v2
  double precision, intent(in) ::  w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu
  double precision, intent(in) :: kstardesired, ustardesired, lstardesired, rhostardesired
  double precision, dimension(9,1), intent(out) :: cost_grad_vector
  double precision, intent(out) :: cost_grad_mag
  double precision :: k1_minus, k1_plus, k2_minus, k2_plus, u1_minus, u1_plus, u2_minus, u2_plus, l1_minus, l1_plus
  double precision :: l2_minus, l2_plus, rho1_minus, rho1_plus, rho2_minus, rho2_plus, v2_minus, v2_plus
  double precision :: cost_grad_k1, cost_grad_k2, cost_grad_u1, cost_grad_u2, cost_grad_l1, cost_grad_l2
  double precision :: cost_grad_rho1, cost_grad_rho2, cost_grad_v2
  double precision :: kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u
  double precision :: cost_k1_minus, cost_k1_plus, cost_k2_minus, cost_k2_plus, cost_u1_minus, cost_u1_plus
  double precision :: cost_u2_minus, cost_u2_plus, cost_l1_minus, cost_l1_plus, cost_l2_minus, cost_l2_plus
  double precision :: cost_rho1_minus, cost_rho1_plus, cost_rho2_minus, cost_rho2_plus, cost_v2_minus, cost_v2_plus
  !double precision :: temp1, temp2
  integer :: xx, yy, rounds


  
  rounds  = 5


  !Compute the effective properties for the given combination to yield the derivative

 

  !Derivative distances
  k1_minus = k1 - del_k1
  k1_plus = k1 + del_k1
  k2_minus = k2 - del_k2
  k2_plus = k2 + del_k2
  u1_minus = u1 - del_u1
  u1_plus = u1 + del_u1
  u2_minus = u2 - del_u2
  u2_plus = u2 + del_u2
  l1_minus = l1 - del_l1
  l1_plus = l1 + del_l1
  l2_minus = l2 - del_l2
  l2_plus = l2 + del_l2
  rho1_minus = rho1 - del_rho1
  rho1_plus = rho1 + del_rho1
  rho2_minus = rho2 - del_rho2
  rho2_plus = rho2 + del_rho2
  v2_minus = v2 - del_v2
  v2_plus = v2 + del_v2

  

  !Run the calculation loops for each derivative term

   
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For k1 - First component

  

   
  !Minus
  call bounds_number(k1_minus,k2,u1,u2,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)

  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_k1_minus)


    
   !Plus
   call bounds_number(k1_plus,k2,u1,u2,l1,l2,rho1,rho2,v2, &
        kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)

   

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_k1_plus)


   
   !Rounding

   
  

  cost_grad_k1 = (cost_k1_plus - cost_k1_minus)/(2.0d0*del_k1)

 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !print *, 'Component 2'
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For k2 - Second component

  
  !Minus
  call bounds_number(k1,k2_minus,u1,u2,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_k2_minus)

   

   !Plus
   call bounds_number(k1,k2_plus,u1,u2,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_k2_plus)

   

   cost_grad_k2 = (cost_k2_plus - cost_k2_minus)/(2.0d0*del_k2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For u1 - Third component

  
  !Minus
  call bounds_number(k1,k2,u1_minus,u2,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_u1_minus)

   

   !Plus
   call bounds_number(k1,k2,u1_plus,u2,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_u1_plus)

   !Rounding

   !xx = anint(cost_u1_minus*10.0**rounds)
   !cost_u1_minus = xx / 10.0**rounds

   !yy = anint(cost_u1_plus*10.0**rounds)
   !cost_u1_plus = yy / 10.0**rounds

   cost_grad_u1 = (cost_u1_plus - cost_u1_minus)/(2.0d0*del_u1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For u2 - Fourth component

  
  !Minus
  call bounds_number(k1,k2,u1,u2_minus,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_u2_minus)

   !Plus
   call bounds_number(k1,k2,u1,u2_plus,l1,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_u2_plus)

    !Rounding

   !xx = anint(cost_u2_minus*10.0**rounds)
   !cost_u2_minus = xx / 10.0**rounds

   !yy = anint(cost_u2_plus*10.0**rounds)
   !cost_u2_plus = yy / 10.0**rounds


   cost_grad_u2 = (cost_u2_plus - cost_u2_minus)/(2.0d0*del_u2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For l1 - Fifth component

  
  !Minus
  call bounds_number(k1,k2,u1,u2,l1_minus,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_l1_minus)

   !Plus
   call bounds_number(k1,k2,u1,u2,l1_plus,l2,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_l1_plus)

  

   cost_grad_l1 = (cost_l1_plus - cost_l1_minus)/(2.0d0*del_l1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For l2 - Sixth component

  
  !Minus
  call bounds_number(k1,k2,u1,u2,l1,l2_minus,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_l2_minus)

   !Plus
   call bounds_number(k1,k2,u1,u2,l1,l2_plus,rho1,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_l2_plus)

  

   cost_grad_l2 = (cost_l2_plus - cost_l2_minus)/(2.0d0*del_l2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For rho1 - Seventh component

  
  !Minus
  call bounds_number(k1,k2,u1,u2,l1,l2,rho1_minus,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_rho1_minus)

   !Plus
   call bounds_number(k1,k2,u1,u2,l1,l2,rho1_plus,rho2,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_rho1_plus)

   

   cost_grad_rho1 = (cost_rho1_plus - cost_rho1_minus)/(2.0d0*del_rho1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    

   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For rho2 - Eight component

  
  !Minus
  call bounds_number(k1,k2,u1,u2,l1,l2,rho1,rho2_minus,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_rho2_minus)

   !Plus
   call bounds_number(k1,k2,u1,u2,l1,l2,rho1,rho2_plus,v2, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_rho2_plus)


   cost_grad_rho2 = (cost_rho2_plus - cost_rho2_minus)/(2.0d0*del_rho2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !For v2 - Ninth component

  
   !Minus



   
  call bounds_number(k1,k2,u1,u2,l1,l2,rho1,rho2,v2_minus, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_v2_minus)
   
   

   !Plus
   call bounds_number(k1,k2,u1,u2,l1,l2,rho1,rho2,v2_plus, &
       kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u)  

   call cost_number(kstar, ustar, lstar, rhostar, c1k, c1u, c2k, c2u, &
           w1, w2, w3, w4, w5hat, w6hat, w7hat, w8hat,p,q, tolk, tolu, &
           kstardesired, ustardesired, lstardesired, rhostardesired, &
           cost_v2_plus)

   
  

   

   cost_grad_v2 = (cost_v2_plus - cost_v2_minus)/(2.0d0*del_v2)

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  


   !Put the gradient into one vector and calculate its magnitude


   cost_grad_vector(:,1) = (/ cost_grad_k1, cost_grad_k2, cost_grad_u1, cost_grad_u2, &
        cost_grad_l1, cost_grad_l2, cost_grad_rho1, cost_grad_rho2, cost_grad_v2 /)

   cost_grad_mag = ( cost_grad_k1**2.0d0 + cost_grad_k2**2.0d0 + cost_grad_u1**2.0d0 &
        + cost_grad_u2**2.0d0 + cost_grad_l1**2.0d0 + cost_grad_l2**2.0d0 + cost_grad_rho1**2.0d0 &
        + cost_grad_rho2**2.0d0)**0.50d0 !+ cost_grad_v2**2.0d0 )**0.50d0

  



end subroutine cost_gradient








subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t
          
  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed


