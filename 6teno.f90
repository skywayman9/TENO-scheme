!!!    This program solves Double Mach reflection problem of 2D Euler equations TENO scheme
!**********************************************************************
! Written by Sainath Ch, Email: s.chamarthi@al.t.u-tokyo.ac.jp
!**********************************************************************


	program doublemach5

	implicit none


    integer 				:: i, j
	integer					:: N, rk_step
	integer					:: NTMAX=100000

	double precision		:: t_end, time , dt

	double precision		:: CFL =0.4d0, gamma =1.4d0

	integer, parameter 		:: NX = 768, NY = 256, ghostp = 5, n_eqn = 4


	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, xmin, xmax, ymin, ymax

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn,0:2) ! 0:2 is for the three Runge-Kutta time steps


	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)

    integer, parameter		:: file_save=200

	integer 				:: time_ini,time_end,time_calc
	double precision 		:: start, finish

	common /grid/ dx, dy


	call cpu_time(start)
	write(*,*) 'Program start...'

	call system_clock(count = time_ini)


	xmin = 0.0d0
	xmax = 3.0d0

	ymin = 0.0d0
	ymax = 1.0d0

	t_end = 0.20
! Generate simple grid

	
    dx = (xmax - xmin)/NX
	dy = (ymax - ymin)/NY

	do i = -ghostp, NX + ghostp

		  x(i) = xmin + (i-0.5d0)*dx
	enddo

	do j = -ghostp, NY + ghostp

	  	y(j) = ymin + (j-0.5d0)*dy
	
	enddo	

	call initialconditions(x,y,density, u_vel, v_vel,pressure, sound, gamma,NX,NY,ghostp)
	call timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)


	time = 0.0d0
	N=1
	
	! call output(density,u_vel,v_vel,pressure)
	call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

	  write(*,*)'*********************************************'
      write(*,*)'   time step N        time             '
      write(*,*)'*********************************************'

      ! Computations starts here

	do while(time.lt.t_end)
			call timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)
			
			time = time + dt

			write(*,*) N ,time

			do rk_step = 0,2

				call boundaryconditions(density, u_vel, v_vel, pressure, x, y, time, gamma,NX,NY,ghostp)

				call FX(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp)

				call GY(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp)
				
				
				call rungekutta(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step, dt,NX,NY,ghostp,n_eqn)

			enddo

			N=N+1
			

			if(MOD(N,file_save) .eq. 0) then

				call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)
				! call tecplot(N/file_save,density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

			endif	

			if (abs(time-t_end) .le. 1.0d-06) then
				
				call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

				write(*,*)'*********************************************'
           		write(*,*)'   Number of time steps = ',N
          	    write(*,*)'*********************************************'

          	    exit
          	endif

    enddo
    	write(*,*)'*********************************************'
        write(*,*)'   Number of time steps = ',N
        write(*,*)'*********************************************'
    	
    		call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)
    		call tecplot(N/file_save,density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

			call system_clock(count = time_end)
    
   			time_calc = time_end-time_ini
    
    	write(*,'(A20,I10,A)')   'Calculation time ',time_calc,' [CPU ticks]'

    		call cpu_time(finish)
    
    	write(*,*) " Total CPU time to solution = ", finish-start, " seconds"

    	write(*,*) 'Program ends...'


	end program doublemach5


	!***********************************************************************
	!*****                       Initial conditions                    *****
	!***********************************************************************


	subroutine initialconditions(x,y,density, u_vel, v_vel,pressure, sound, gamma,NX,NY,ghostp)

	implicit none

	integer 				:: i, j

	integer			 		:: NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, xmin, xmax, ymin, ymax, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0)

	do i = 1, NX

		do j =  1, NY

			if(y(j) .lt. (x(i) - (1.0d0/6.0d0))*sqrt(3.0d0) )then

				density(i,j) 	= 1.4d0 
				u_vel(i,j)		= 0.0d0
				v_vel(i,j)		= 0.0d0
				pressure(i,j)	= 1.0d0

			else

				density(i,j)	= 8.0d0
				u_vel(i,j)		= 8.25*(cos(pi/6))
				v_vel(i,j)		= -8.25*(sin(pi/6))
				pressure(i,j)	= 116.5d0

			endif
			
			sound(i,j)			= (gamma*pressure(i,j)/density(i,j))**0.5
		enddo	
	enddo	


	end subroutine initialconditions

	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

	integer 				:: i, j

	integer			 		:: NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)


	open(unit=25, file="soln.txt",action="write",status="replace")
    do j = 1, NY
        do i = 1, NX
	 
	  write(25,'(6F25.8)') x(i),y(j),density(i,j),pressure(i,j),u_vel(i,j),v_vel(i,j)
	 enddo
	enddo
	close(25)


 	end subroutine output

 	subroutine tecplot(file,density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)
	implicit none


    integer 				:: i, j, file,l

	integer			 		:: NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	character(len=8) 		:: number*4, file_name

	write(number,'(i4.4)') file
	file_name="Rslt"//number
	open(unit=1,file=file_name//'.plt')
	
	write(1,*) 'TITLE="',file_name,'"'
    write(1,*) 'VARIABLES = "x","y","rho","vx","vy","Pre"'
	write(1,*) "ZONE I=",NX," J=",NY," F=POINT"


    do j = 1, NY
        do i = 1, NX

          write(1,'(6F25.8)') x(I), y(J), density(I,J), u_vel(I,J), v_vel(I,J), pressure(I,J)

	  enddo
    enddo
	
    close(1)



    end subroutine tecplot


 	!***********************************************************************
	!*****                       Compute time step                	   *****
	!***********************************************************************

 	subroutine timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)
 	implicit none


 	integer 				:: i, j

	integer			 		:: NX, NY, ghostp
	double precision 		:: dx, dy,gamma

	double precision 		:: dt, time, t_end, dtnew, CFL

 	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
 	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: x_velocity, y_velocity

	common /grid/ dx, dy

	dt = 1.0d10

		do i = 0, NX

			do j = 0, NY
				sound(i,j) =  (gamma*pressure(i,j)/density(i,j))**0.5d0
				x_velocity =  ABS(u_vel(i,j)) + sound(i,j)
				y_velocity =  ABS(v_vel(i,j)) + sound(i,j)

				dtnew = min(dx/x_velocity, dy/y_velocity)

				if(dtnew .lt. dt) dt = dtnew
			enddo
		enddo

		dt = CFL*dt

		if ((time+dt) .gt. t_end ) then

			dt = t_end - time
		endif	


 	end subroutine timestep


 	!***********************************************************************
	!*****                       Boundary conditions                   *****
	!***********************************************************************


 	subroutine boundaryconditions(density, u_vel, v_vel, pressure, x, y, time, gamma,NX,NY,ghostp)
    implicit none

	integer 				:: i, j, NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, x1, y1, gamma

	double precision		:: time

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0)

	double precision		:: shockspeed=10.0d0 		!  The shock is propagating with a speed of Mach 10

	common /mesh/ dx, dy

	

	! Set top boundary condition

	do j = 1,ghostp

		do i = 1, NX

			x1 = x(i) + 0.5d0*dx
			y1 = y(NY) + 0.5d0*dy


			if(y1-sqrt(3.0d0)*(x1-1.0d0/6.0d0).gt.-shockspeed*time/0.5) then ! Time dependent boundary conditions

				density(i,NY+j)		= 8.0d0
				u_vel(i,NY+j)		= 8.25*(cos(pi/6))
				v_vel(i,NY+j)		= -8.25*(sin(pi/6))
				pressure(i,NY+j)	= 116.5d0

			else
				density(i,NY+j) 	= 1.4d0 
				u_vel(i,NY+j)		= 0.0d0
				v_vel(i,NY+j)		= 0.0d0
				pressure(i,NY+j)	= 1.0d0
			endif


		enddo
	enddo

	! Set bottom boundary condition

	do j = 1,ghostp

		do i = 1, NX

			x1 = x(i) + 0.5d0*dx
			y1 = y(0) + 0.5d0*dy


			if(x1 .lt. 1.0d0/6.0d0) then ! Post-shock conditions

			density(i,-j+1)	= 8.0d0
			u_vel(i,-j+1)		= 8.25*(cos(pi/6))
			v_vel(i,-j+1)		= -8.25*(sin(pi/6))
			pressure(i,-j+1)	= 116.5d0

	    else
	      	density(i,-j+1)	=  density(i,j-1+1)
	      	u_vel(i,-j+1)		=  u_vel(i,j-1+1)
	      	v_vel(i,-j+1)		= -v_vel(i,j-1+1)		! Normal velocity is reversed. Wall boundary condition
	      	pressure(i,-j+1)	=  pressure(i,j-1+1)
	    endif

	  enddo
	enddo	


	! Set left and right boundary conditions. Gradient boundary conditions
	do i = 1,ghostp

		do j = 1, NY

		! left -1=0, -2=1, -3=2

		! left 0=1, -1=2 and -2=3

		density(-i+1,j)	= density(i,j)
		u_vel(-i+1,j)		= u_vel (i,j)
		v_vel(-i+1,j)		= v_vel (i,j)
		pressure(-i+1,j)	= pressure(i,j)


		! right NX+1=NX, NX+2=NX-1, NX+3=NX-2
		
		density(NX+i,j)		= density(NX-i+1,j)
		u_vel(NX+i,j)		= u_vel(NX-i+1,j)
		v_vel(NX+i,j)		= v_vel(NX-i+1,j)
		pressure(NX+i,j)	= pressure(NX-i+1,j)

		enddo
	enddo
     
      end

 	!***********************************************************************
	!*****                       Time step, TVD- Runge Kutta           *****
	!***********************************************************************

	subroutine rungekutta(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step, dt,NX,NY,ghostp,n_eqn)
	implicit none

	integer 				:: i, j, k, rk_step, rk

	integer 				:: NX, NY, ghostp,n_eqn

	double precision 		:: dt, uv, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn,0:2) ! 0:2 is for the three Runge-Kutta time steps

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)


	if(rk_step .EQ.0) then

	do k=1, n_eqn

	    do j = 1, NY
	
		  do i = 1, NX
	
		    cons(i,j,k,1) = cons(i,j,k,0) + dt*residual(i,j,k)
	
		  enddo
	
		enddo
	enddo

	elseif(rk_step .EQ.1) then

	do k=1, n_eqn
	    
	    do j = 1, NY
		
		  do i = 1, NX
		
		    cons(i,j,k,2) = (3.0d0/4.0d0)*cons(i,j,k,0) + (1.0/4.0d0)*(cons(i,j,k,1) + dt*residual(i,j,k))
		
		  enddo
		
		enddo
	enddo

	else

	do k=1, n_eqn
	    
	    do j = 1, NY
		
		  do i = 1, NX
		
		    cons(i,j,k,0) = (1.0d0/3.0d0)*cons(i,j,k,0) + (2.0d0/3.0d0)*(cons(i,j,k,2) + dt*residual(i,j,k))
		
		  enddo
		
		enddo
	
	enddo

	endif
	
	rk = MOD(rk_step +1, 3)

			do i = 1, NX
				do j = 1, NY


			    density(i,j)		= cons(i,j,1,rk)
		       	u_vel(i,j)			= cons(i,j,2,rk)/density(i,j)
		        v_vel(i,j)			= cons(i,j,3,rk)/density(i,j)
		        uv 					= u_vel(i,j)**2 + v_vel(i,j)**2
			    pressure(i,j)		= (gamma-1.0)*(cons(i,j,4,rk)-0.5*cons(i,j,1,rk)*uv)

				enddo
			enddo

	end subroutine rungekutta



 !***********************************************************************
	!*****                      Flux in X-direction                    *****
	!***********************************************************************


 	subroutine FX(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp)
 	implicit none


 	integer     		:: i, j, rk_step, M, NX, NY, ghostp,k
 	integer, parameter	::  n_eqn =4

 	double precision    ::	 dx, dy, gamma



	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: sound(-ghostp:NX+ghostp, -ghostp:NY+ghostp),enthalpy(-ghostp:NX+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4,0:2),constime(-ghostp:NX+ghostp,4), primitive(-ghostp:NX+ghostp,4) ! 0:2 is for the three Runge-Kutta time steps

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4)

	double precision		:: lefteigen(-ghostp:NX+ghostp,4,4) , righteigen(-ghostp:NX+ghostp,4,4) 		! 2D 4 by 4 matrix 4 eigen values and .....

	double precision		:: consl(-ghostp:NX+ghostp,4),consr(-ghostp:NX+ghostp,4)

	double precision		:: fluxl(-ghostp:NX+ghostp,4),fluxr(-ghostp:NX+ghostp,4)

	double precision		:: priml(-ghostp:NX+ghostp,4),primr(-ghostp:NX+ghostp,4)

	double precision		:: flux(-ghostp:NX+ghostp,4), flux_half(-ghostp:NX+ghostp,4)

	double precision		:: fright(-ghostp:NX+ghostp,4), fleft(-ghostp:NX+ghostp,4)

	double precision		:: lambda1(-ghostp:NX+ghostp),lambda2(-ghostp:NX+ghostp),lambda3(-ghostp:NX+ghostp),lambda4(-ghostp:NX+ghostp)

	double precision		:: den_average(-ghostp:NX+ghostp)
	double precision   		:: enthalpy_average, uv_average, sound_average, u_average, v_average,sound_left,sound_right

	
	double precision   		:: T0, T1, T2, T3, B1, B2,inverse_sound,Da, ML,MR

	double precision		:: lambda(4)

	!!!! Riemann solver

	double precision		:: densityleft(-ghostp:NX+ghostp), pressureleft(-ghostp:NX+ghostp),u_velleft(-ghostp:NX+ghostp), v_velleft(-ghostp:NX+ghostp)

	double precision		:: densityright(-ghostp:NX+ghostp), pressureright(-ghostp:NX+ghostp),u_velright(-ghostp:NX+ghostp), v_velright(-ghostp:NX+ghostp)
	double precision		:: enthalpyleft(-ghostp:NX+ghostp), enthalpyright(-ghostp:NX+ghostp),energyleft(-ghostp:NX+ghostp),energyright(-ghostp:NX+ghostp)


	double precision 		:: left(-ghostp:NX+ghostp,4),right(-ghostp:NX+ghostp,4),mx,my



    

	common /grid/ dx, dy


	mx=1.0d0;my=0.0d0
	

	do j = 1, NY



		do i = -ghostp, NX + ghostp

				cons(i,j,1,rk_step ) = density(i,j)
			    cons(i,j,2,rk_step ) = density(i,j)*u_vel(i,j)
			    cons(i,j,3,rk_step ) = density(i,j)*v_vel(i,j)
			    cons(i,j,4,rk_step ) = pressure(i,j)/(gamma-1.0d0) + 0.5d0*density(i,j)*(u_vel(i,j)**2.0d0 + v_vel(i,j)**2.0d0)

			    constime(i,1) = density(i,j)
			    constime(i,2) = density(i,j)*u_vel(i,j)
			    constime(i,3) = density(i,j)*v_vel(i,j)
			    constime(i,4) = pressure(i,j)/(gamma-1.0d0) + 0.5d0*density(i,j)*(u_vel(i,j)**2.0d0 + v_vel(i,j)**2.0d0)

			    primitive(i,1) = density(i,j)
			    primitive(i,2) = u_vel(i,j)
			    primitive(i,3) = v_vel(i,j)
			    primitive(i,4) = pressure(i,j)

			    flux(i,1) = constime(i,2)									! rho * u
			    flux(i,2) = constime(i,2) * u_vel(i,j) + pressure(i,j)		! rho * u * u + p
			    flux(i,3) = constime(i,2) * v_vel(i,j)						! rho * u * v
			    flux(i,4) = u_vel(i,j) * (constime(i,4) + pressure(i,j))

			      sound(i,j) = sqrt(gamma*pressure(i,j)/density(i,j))

			   	lambda(1) = ABS(u_vel(i,j)) - sound(i,j)
			    lambda(2) = ABS(u_vel(i,j))
			    lambda(3) = lambda(2)
			    lambda(4) = ABS(u_vel(i,j)) + sound(i,j)


			left(i,:)=0.5d0*(flux(i+0,:) + (lambda(4))*(constime(i+0,:)))
			right(i,:)=0.5d0*(flux(i+0,:) - (lambda(4))*(constime(i+0,:)))

		enddo
	  
	  call WCNS_flux(NX, constime, flux, fluxl,fluxr,mx,my,ghostp,n_eqn,left,right)


	  do i =-1,NX+1


			flux_half(i,:) =  (fluxr(i,:) +fluxl(i,:))
      
		
		enddo


	  do M = 1, n_eqn
	    do i = 1, NX
		  residual(i,j,M) = -(flux_half(i,M) - flux_half(i-1,M))/dx
		enddo
	  enddo
	
	enddo
 	

 	end subroutine FX

 	!***********************************************************************
	!*****                      Flux in Y-direction                    *****
	!***********************************************************************

	subroutine GY(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp)
	implicit none

 	integer     		:: i, j, rk_step, M, NX, NY, ghostp,k
 	integer, parameter	::  n_eqn =4

 	double precision    ::	 dx, dy, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: sound(-ghostp:NX+ghostp, -ghostp:NY+ghostp),enthalpy(-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4,0:2),constime(-ghostp:NY+ghostp,4),primitive(-ghostp:NY+ghostp,4) ! 0:2 is for the three Runge-Kutta time steps

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4)

	double precision		:: lefteigen(-ghostp:NY+ghostp,4,4) , righteigen(-ghostp:NY+ghostp,4,4) 		! 2D 4 by 4 matrix 4 eigen values and .....

	double precision		:: consl(-ghostp:NY+ghostp,4),consr(-ghostp:NY+ghostp,4)


	double precision		:: fluxl(-ghostp:NY+ghostp,4),fluxr(-ghostp:NY+ghostp,4)


	double precision		:: priml(-ghostp:NY+ghostp,4),primr(-ghostp:NY+ghostp,4)

	double precision		:: flux(-ghostp:NY+ghostp,4), flux_half(-ghostp:NY+ghostp,4)

	double precision		:: fright(-ghostp:NY+ghostp,4), fleft(-ghostp:NY+ghostp,4)

	double precision		:: lambda1(-ghostp:NY+ghostp),lambda2(-ghostp:NY+ghostp),lambda3(-ghostp:NY+ghostp),lambda4(-ghostp:NY+ghostp)

	double precision		:: den_average(-ghostp:NY+ghostp)
	double precision   		:: enthalpy_average, uv_average, sound_average, u_average, v_average, sound_right,sound_left

	
	double precision   		:: T0, T1, T2, T3, B1, B2,inverse_sound,Da,ML,MR

	double precision  :: lambda(4)

	!!!! Riemann solver
	double precision 		:: left(-ghostp:NY+ghostp,4),right(-ghostp:NY+ghostp,4),mx,my



	common /grid/ dx, dy


	mx=0.0d0;my=1.0d0
	

	do j = 1, NX


	 		do i = -ghostp, NY + ghostp

				cons(j,i,1,rk_step ) = density(j,i)
			    cons(j,i,2,rk_step ) = density(j,i)*u_vel(j,i)
			    cons(j,i,3,rk_step ) = density(j,i)*v_vel(j,i)
			    cons(j,i,4,rk_step ) = pressure(j,i)/(gamma-1.0d0) + 0.5d0*density(j,i)*(u_vel(j,i)**2.0d0 + v_vel(j,i)**2.0d0)

			   	constime(i,1) = density(j,i)
			    constime(i,2) = density(j,i)*u_vel(j,i)
			    constime(i,3) = density(j,i)*v_vel(j,i)
			    constime(i,4) = pressure(j,i)/(gamma-1.0d0) + 0.5d0*density(j,i)*(u_vel(j,i)**2.0d0 + v_vel(j,i)**2.0d0)

			   	primitive(i,1) = density(j,i)
			    primitive(i,2) = u_vel(j,i)
			    primitive(i,3) = v_vel(j,i)
			    primitive(i,4) = pressure(j,i)

				flux(i,1) = constime(i,3)
			    flux(i,2) = constime(i,2) * v_vel(j,i)
			    flux(i,3) = constime(i,3) * v_vel(j,i) + pressure(j,i)
			    flux(i,4) = v_vel(j,i) * (constime(i,4) + pressure(j,i))

			   	sound(j,i) = sqrt(gamma*pressure(j,i)/density(j,i))

			    lambda(1) = ABS(v_vel(j,i)) - sound(j,i)
			    lambda(2) = ABS(v_vel(j,i))
			    lambda(3) = lambda(2)
			    lambda(4) = ABS(v_vel(j,i)) + sound(j,i)

			left(i,:)=0.5d0*(flux(i+0,:) + (lambda(4))*(constime(i+0,:)))
			right(i,:)=0.5d0*(flux(i+0,:) - (lambda(4))*(constime(i+0,:)))


			enddo

	  call WCNS_flux(NY, constime, flux, fluxl,fluxr,mx,my,ghostp,n_eqn,left,right)

	  do i =-1, NY+1
	  		
			flux_half(i,:) =  (fluxr(i,:) +fluxl(i,:))
		enddo




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	  do M = 1, n_eqn
	    do i = 1, NY
	      residual(j,i,M) = residual(j,i,M) - (flux_half(i,M) - flux_half(i-1,M))/dy
		enddo
	  enddo
	enddo


	end subroutine GY
	



	!***********************************************************************
	!*****                      WCNS interpolation	                   *****
	!***********************************************************************
	subroutine WCNS_flux (NS, un, flux, ulnew,urnew,mx,my,ghostp,n_eqn,left,right)

	implicit none
	integer				:: ix,NS, ghostp, n_eqn,i,j,k,total
	double precision	:: un(-ghostp:NS+ghostp,n_eqn),flux(-ghostp:NS+ghostp,n_eqn)


	double precision	:: ulnew(-ghostp:NS+ghostp,n_eqn),urnew(-ghostp:NS+ghostp,n_eqn)

	double precision	:: left(-ghostp:NS+ghostp,n_eqn),right(-ghostp:NS+ghostp,n_eqn)


	double precision 	:: mx,my, lx, ly
	double precision 	:: interp_poly5(0:3), alpha(0:3), beta(0:3), omega(0:3)

			double precision	:: v(-3:4,n_eqn),charstencil(-3:3)

    double precision 	:: lefteigen(4,4), righteigen(4,4)
    double precision, parameter	:: epsilon = 1.0d-40, gamma = 1.4d0, Constant =1.0d00



	double precision    :: ul(n_eqn),ur(n_eqn)

  	double precision    :: sqrt_rho, divisor,rho,u_vel,v_vel,p,c,sqrt_rho_L,sqrt_rho_R

  	double precision, PARAMETER ::  B2 = 4.0d0/3.0d0, MP5=4.0d0, EPSM=1.0d-40
	double precision :: VOR,VMP,DJM1,DJ,DJP1,DM4JPH,DM4JMH,VUL,VAV,VMD,VLC,VMIN,VMAX, MINMOD2,MINMOD4,l2norm

	double precision 	:: du_dx(-4:NS+4,n_eqn),du2_dx2(-4:NS+4,n_eqn)

	double precision, parameter	:: power = 6.0d0,CT=1.0d-5


    double precision 	::  tau, teno(0:3),temp


        double precision  :: vl,vr,pl,pr,hl,hr,enth,uv_average,c2,qn,ql,q2

        double precision :: u_right, v_right,pre_right,u_left,v_left,pre_left,uv,rol,ror


	lx = -my ; ly = mx;

		



	do ix =-1,NS+1


		rol=(un(ix,1))**0.5
        ror=(un(ix+1,1))**0.5
        rho=rol*ror

        u_vel 	= (rol*(un(ix,2)/un(ix,1))+ror*(un(ix+1,2)/un(ix+1,1)))/(rol+ror)


        v_vel 	= (rol*(un(ix,3)/un(ix,1))+ror*(un(ix+1,3)/un(ix+1,1)))/(rol+ror)

		uv_average 		= 	(u_vel**2.0d0 + v_vel**2.0d0)

        pl=(gamma-1.d0)*(un(ix,4)  -0.5*un(ix,1)  *(uv_average))
        pr=(gamma-1.d0)*(un(ix+1,4)-0.5*un(ix+1,1)*(uv_average))
         p=(rol*pl+ror*pr)/(rol+ror)

        ! if(p.lt. 0.0) stop

        ! if((gamma*p/rho) .lt. 0.0d0) stop

        hl=(un(ix,4)  + pl)/un(ix,1)
        hr=(un(ix+1,4)+ pr)/un(ix+1,1)

		enth=(hl*rol+hr*ror)/(rol+ror)

		
		c 	= 	dsqrt((gamma - 1.0d0) * (enth - 0.5d0*uv_average))	
		! c       = sqrt(gamma*p/rho)

		c2 	= 	c**2.0d0

		if(p.lt.0.0d0 .or. rho.lt.0.0d0 .or. c .lt. 0.0d0)then

		write(*,*)p,rho,c

		write(*,*)ix

		write(*,*) (gamma-1.d0)*(un(ix,4)  -0.5*un(ix,1)  * ((un(ix,2)/un(ix,1))**2.0d0 + (un(ix,3)/un(ix,1))**2.0d0 ) ),pl,pr
		! stop
		endif

		qn=u_vel*mx+v_vel*my
		ql=u_vel*lx+v_vel*ly
		q2=u_vel*u_vel+v_vel*v_vel


	! lx = -my ; ly = mx;
        righteigen(1,1)=	1.0d0;	 	righteigen(1,2)= 0.0d0;	righteigen(1,3)=	1.0d0;		righteigen(1,4)=	1.0d0;
        righteigen(2,1)=u_vel-c*mx;	 	righteigen(2,2)=lx;		righteigen(2,3)=	u_vel;		righteigen(2,4)= u_vel+c*mx;
        righteigen(3,1)=v_vel-c*my;	 	righteigen(3,2)=ly;		righteigen(3,3)=	v_vel;		righteigen(3,4)= v_vel+c*my;
        righteigen(4,1)=enth-qn*c;		righteigen(4,2)=ql;     righteigen(4,3)=   0.5*q2;		righteigen(4,4)= enth+qn*c;

        
        lefteigen(1,1)=0.5d0*(((gamma-1.0d0)*0.5d0*q2/c2)+(qn/c));	
        lefteigen(1,2)=-0.5d0*(((gamma-1.0d0)*u_vel/c2)+(mx/c));
        lefteigen(1,3)=-0.5d0*(((gamma-1.0d0)*v_vel/c2)+(my/c));
        lefteigen(1,4)= (gamma-1.0d0)/(2.0d0*c2);

        lefteigen(2,1)=-ql;	
        lefteigen(2,2)=lx;	
        lefteigen(2,3)=ly; 
        lefteigen(2,4)=0.0d0;

        lefteigen(3,1)= 1.0d0- ((gamma-1.0d0)*q2/(2.0d0*c2));	
        lefteigen(3,2)= (gamma-1.0d0)*u_vel/c2;	
        lefteigen(3,3)= (gamma-1.0d0)*v_vel/c2;	
        lefteigen(3,4)= -((gamma-1.0d0)/(c2));	


        lefteigen(4,1)=0.5d0*(((gamma-1.0d0)*0.5d0*q2/c2)-(qn/c));	
        lefteigen(4,2)=-0.5d0*(((gamma-1.0d0)*u_vel/c2)-(mx/c));
        lefteigen(4,3)=-0.5d0*(((gamma-1.0d0)*v_vel/c2)-(my/c));
        lefteigen(4,4)= (gamma-1.0d0)/(2.0d0*c2);


        do i=-3,4
		! 
			v(i,:)   = matmul(lefteigen,left(i+ix,:))

		enddo

		do i=1,n_eqn


			charstencil=v(-3:3,i)

			interp_poly5(0)=charstencil(0)*(2.0/6.0d0)+charstencil(1)*(5.0/6.0d0)+charstencil(2)*(-1.0/6d0)
            interp_poly5(1)=charstencil(-1)*(-1.0/6.0d0)+charstencil(0)*(5.0/6.0d0)+charstencil(1)*(2.0/6.0d0)
            interp_poly5(2)=charstencil(-2)*(2.0/6.0d0)+charstencil(-1)*(-7.0/6.0d0)+charstencil(0)*(11.0/6.0d0)
            interp_poly5(3)=1.0d0/12.0d0*( 3.0d0*charstencil(0)+ 13.0d0*charstencil(0+1)-5.0d0*charstencil(0+2)+charstencil(0+3))

     
            beta(0)=13.d0/12.*(charstencil(0)-2.*charstencil(1)+charstencil(2))**2&
                    +(1.0/4.0d0)*(3.*charstencil(0)-4.*charstencil(1)+charstencil(2))**2
            beta(1)=13.d0/12.*(charstencil(-1)-2.*charstencil(0)+charstencil(1))**2&
                    +(1.0/4.0d0)*(charstencil(-1)-charstencil(1))**2
            beta(2)=13.d0/12.*(charstencil(-2)-2.*charstencil(-1)+charstencil(0))**2&
                    +(1.0/4.0d0)*(charstencil(-2)-4.*charstencil(-1)+3.*charstencil(0))**2
			
			beta(3)= (1.0d0/240.d0)*(+2107.0d0*charstencil(0)**2-9402.0d0*charstencil(0)*charstencil(0+1)+7402.0d0*charstencil(0)*charstencil(0+2)-1854.0d0*charstencil(0)*charstencil(0+3)+&
											 11003.0d0*charstencil(0+1)*charstencil(0+1)-17246.0d0*charstencil(0+1)*charstencil(0+2)+4642.0d0*charstencil(0+3)*charstencil(0+1)+7043.0d0*charstencil(0+2)*charstencil(0+2)&
											 -3882.0d0*charstencil(0+2)*charstencil(0+3)+547.0d0*charstencil(0+3)*charstencil(0+3))


		  	temp = (1.0d0/120960.0d0)*(271779.0d0*charstencil(-2)**2.d0+&
		  				charstencil(-2)*(-2380800.0d0*charstencil(-1)+4086352.00d0*charstencil(0)-3462252.0d0*charstencil(1)+1458762.0d0*charstencil(2)-245620.0d0*charstencil(3))+&
		  				charstencil(-1)*(+5653317.0d0*charstencil(-1)-20427884.0d0*charstencil(0)+1705032.0d0*charstencil(1)-7727988.0d0*charstencil(2)+1325006.0d0*charstencil(3))+&
		  				charstencil(0) *(19510972.0d0*charstencil(0)-35817664.0d0*charstencil(1)+15929912.0d0*charstencil(2)-2792660.0d0*charstencil(3))+&
		  				charstencil(1) *(17195652.0d0*charstencil(1)-15880404.0d0*charstencil(2)+2863984.0d0*charstencil(3))+&
		  				charstencil(2) *(3824847.0d0*charstencil(2)-1429976.0d0*charstencil(3)+139633.0d0*charstencil(3)**2.0d0))


		 	tau = ABS(temp - 1.0/6.d0*(beta(0)+beta(2)+4.0d0*beta(1))) 

! !WENO-Z by Borges
            alpha(0) = (Constant + (abs(tau) / (epsilon + beta(0))))**power
            alpha(1) = (Constant + (abs(tau) / (epsilon + beta(1))))**power
            alpha(2) = (Constant + (abs(tau) / (epsilon + beta(2))))**power
            alpha(3) = (Constant + (abs(tau) / (epsilon + beta(3))))**power

	            teno(0)  =alpha(0)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
	            teno(1)  =alpha(1)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
	            teno(2)  =alpha(2)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
	            teno(3)  =alpha(3)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))

	            do j=0,3
	            if(teno(j) .le. CT) then
	            teno(j) =0.0d0 
	            else
	            teno(j) =1.0d0
	            endif
	            enddo
            
           	alpha(0) = (0.30d0)*teno(0)
            alpha(1) = (0.45d0)*teno(1)
            alpha(2) = (0.05d0)*teno(2)
            alpha(3) = (0.20d0)*teno(3)
            
            omega(0)=alpha(0)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            omega(1)=alpha(1)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            omega(2)=alpha(2)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            omega(3)=alpha(3)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))


            ul(i)=omega(0)*interp_poly5(0)+omega(1)*interp_poly5(1)+omega(2)*interp_poly5(2)+omega(3)*interp_poly5(3)
					
        enddo

        do i=-3,4
		! 
			v(i,:)   = matmul(lefteigen,right(i+ix,:))

		enddo

		do i=1,n_eqn


			charstencil=v(-3:3,i)
            
            interp_poly5(0)=1.0d0/12.0d0*( 3.0d0*charstencil(0)+ 13.0d0*charstencil(0-1)-5.0d0*charstencil(0-2)+charstencil(0-3))
            interp_poly5(1)=charstencil(2)*(2.0/6.0d0)+charstencil(+1)*(-7.0/6.0d0)+charstencil(0)*(11.0/6.0d0)
            interp_poly5(2)=charstencil(1)*(-1.0/6.0d0)+charstencil(0)*(5.0/6.0d0)+charstencil(-1)*(2.0/6.0d0)
			interp_poly5(3)=charstencil(0)*(2.0/6.0d0)+charstencil(-1)*(5.0/6.0d0)+charstencil(-2)*(-1.0/6d0)
           

			beta(0)		= (1.0d0/240.d0)*(+2107.0d0*charstencil(0)**2-9402.0d0*charstencil(0)*charstencil(0-1)+7402.0d0*charstencil(0)*charstencil(0-2)-1854.0d0*charstencil(0)*charstencil(0-3)+&
											 11003.0d0*charstencil(0-1)*charstencil(0-1)-17246.0d0*charstencil(0-1)*charstencil(0-2)+4642.0d0*charstencil(0-3)*charstencil(0-1)+7043.0d0*charstencil(0-2)*charstencil(0-2)&
											 -3882.0d0*charstencil(0-2)*charstencil(0-3)+547.0d0*charstencil(0-3)*charstencil(0-3))

            beta(1)=13.d0/12.*(charstencil(0)-2.*charstencil(1)+charstencil(2))**2&
                    +(1.0/4.0d0)*(3.*charstencil(0)-4.*charstencil(1)+charstencil(2))**2

            beta(2)=13.d0/12.*(charstencil(-1)-2.*charstencil(0)+charstencil(1))**2&
                    +(1.0/4.0d0)*(charstencil(-1)-charstencil(1))**2

            beta(3)=13.d0/12.*(charstencil(-2)-2.*charstencil(-1)+charstencil(0))**2&
                    +(1.0/4.0d0)*(charstencil(-2)-4.*charstencil(-1)+3.*charstencil(0))**2

		  	temp = (1.0d0/120960.0d0)*(271779.0d0*charstencil(2)**2.d0+&
		  				charstencil(2)*(-2380800.0d0*charstencil(1)+4086352.00d0*charstencil(0)-3462252.0d0*charstencil(-1)+1458762.0d0*charstencil(-2)-245620.0d0*charstencil(-3))+&
		  				charstencil(1)*(+5653317.0d0*charstencil(1)-20427884.0d0*charstencil(0)+1705032.0d0*charstencil(-1)-7727988.0d0*charstencil(-2)+1325006.0d0*charstencil(-3))+&
		  				charstencil(0) *(19510972.0d0*charstencil(0)-35817664.0d0*charstencil(-1)+15929912.0d0*charstencil(-2)-2792660.0d0*charstencil(3))+&
		  				charstencil(-1) *(17195652.0d0*charstencil(-1)-15880404.0d0*charstencil(-2)+2863984.0d0*charstencil(-3))+&
		  				charstencil(-2) *(3824847.0d0*charstencil(-2)-1429976.0d0*charstencil(-3)+139633.0d0*charstencil(-3)**2.0d0))


		 	tau = ABS(temp- 1.0/6.d0*(beta(1)+beta(3)+4.0d0*beta(2))) 

! omega(2)=alpha(2)/(alpha(0)+alpha(1)+alpha(2))

!WENO-Z by Borges
            alpha(0) = (Constant + (abs(tau) / (epsilon + beta(0))))**power
            alpha(1) = (Constant + (abs(tau) / (epsilon + beta(1))))**power
            alpha(2) = (Constant + (abs(tau) / (epsilon + beta(2))))**power
            alpha(3) = (Constant + (abs(tau) / (epsilon + beta(3))))**power

            teno(0)  =alpha(0)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            teno(1)  =alpha(1)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            teno(2)  =alpha(2)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            teno(3)  =alpha(3)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))

            do j=0,3
            if(teno(j) .le. CT) then
            teno(j) =0.0d0 
            else
            teno(j) =1.0d0
            endif
            enddo
            
                       
            alpha(0) = (0.20d0)*teno(0)
           	alpha(1) = (0.05d0)*teno(1)
            alpha(2) = (0.45d0)*teno(2)
            alpha(3) = (0.30d0)*teno(3)


            
            omega(0)=alpha(0)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            omega(1)=alpha(1)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            omega(2)=alpha(2)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))
            omega(3)=alpha(3)/(alpha(0)+alpha(1)+alpha(2)+alpha(3))

            ur(i)=omega(0)*interp_poly5(0)+omega(1)*interp_poly5(1)+omega(2)*interp_poly5(2)+omega(3)*interp_poly5(3)

		enddo

        urnew(ix-1,:)=matmul(righteigen,ur)
        ulnew(ix,:)=matmul(righteigen,ul)
	enddo
!*****************************!*****************************!*****************************!*****************************!*****************************
	
	end subroutine WCNS_flux

  
double precision function minmod2(X,Y)
double precision, INTENT(IN) :: X,Y
minmod2 = 0.5d0*(sign(1.0d0,X)+sign(1.0d0,Y))*min(ABS(X),ABS(Y))
end function minmod2


double precision function minmod4(W,X,Y,Z)
double precision, INTENT(IN) :: W,X,Y,Z

minmod4 = 0.125d0*(sign(1.0d0,W)+sign(1.0d0,X))* &
          ABS( (sign(1.0d0,W)+sign(1.0d0,Y))*(sign(1.0d0,W)+SIGN(1.0d0,Z)) )* &
          min(ABS(W),ABS(X),ABS(Y),ABS(Z))

end function minmod4

