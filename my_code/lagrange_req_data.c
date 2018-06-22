/* advance in time particles and assign them to the right processors */
void move_particles(){

  int ipart,i,j;
  int npart_here,npart_there,all_npart_there,all_npart;  
  point_particle part;
  int *displs,*rcounts;
  my_double invtau;
  vector Dt_u;
  vector vec, v_old;
  vector omega,lift_coeff;
#ifdef LAGRANGE_ORIENTATION
  my_double matA[3][3],matS[3][3],matW[3][3];
  my_double scalOSO, f_alpha, alpha,norm , gyro;
  my_double vecF[3],vecFold[3],vecP[3],vecTMP[3],vecA[3];
 #ifdef LAGRANGE_ORIENTATION_SECONDORIENTATION
  my_double vecFN[3],vecFNold[3],vecN[3];
 #endif
/* #ifdef LAGRANGE_ORIENTATION_DIFFUSION
  my_double vec_xi[3];
  my_double two_d_r;
  my_double matI[3][3]; */
 #endif 
/*  #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP
  my_double shear_rate,jump_time_duration ,velocity_amplitude;
  #endif  */
#endif
  my_double reactivity;

#ifdef LAGRANGE_ORIENTATION
/* #ifdef LAGRANGE_ORIENTATION_DIFFUSION */
  /* Define the identity matrix */
/*	     matI[0][0] = matI[1][1] = matI[2][2] = 1.0;
	     matI[0][1] = matI[0][2] = matI[1][2] = 0.0;
	     matI[1][0] = matI[2][0] = matI[2][1] = 0.0;
 #endif */
#endif

  //fprintf(stderr,"me %d I am here, npart %d time %g\n",me, npart,time_now);

/* Begin loop on particles */
 for (ipart=0;ipart<npart;ipart++) {

   /* We take care of of computing scalar derivatives at particle position */
#ifdef LB_TEMPERATURE
   if(itime==0 && resume==0){ 
     (tracer+ipart)->t_old =  (tracer+ipart)->t;
   }else{
     (tracer+ipart)->dt_t  = ((tracer+ipart)->t - (tracer+ipart)->t_old )/property.time_dt; 
     (tracer+ipart)->t_old =  (tracer+ipart)->t;
   }
#endif
#ifdef LB_SCALAR
   if(itime==0 && resume==0){ 
     (tracer+ipart)->s_old =  (tracer+ipart)->s;
   }else{
     (tracer+ipart)->dt_s  = ((tracer+ipart)->s - (tracer+ipart)->s_old )/property.time_dt; 
     (tracer+ipart)->s_old =  (tracer+ipart)->s;
   }
#endif

   /* This is a trick for Eulerian probes i.e. fixed probes */
   if((tracer+ipart)->tau_drag < 0.0){  /* tau_drag < 0 here conventionally indicate an Eulerian probe */

   (tracer+ipart)->vx = (tracer+ipart)->ux;
   (tracer+ipart)->vy = (tracer+ipart)->uy;
   (tracer+ipart)->vz = (tracer+ipart)->uz;
   }

   /* if we just have tracers */
   if((tracer+ipart)->tau_drag == 0.0){

   /* copy fluid velocity into particle velocity NOTE that this is true only for tracers */
   (tracer+ipart)->vx = (tracer+ipart)->ux;
   (tracer+ipart)->vy = (tracer+ipart)->uy;
   (tracer+ipart)->vz = (tracer+ipart)->uz;

#ifdef LAGRANGE_ORIENTATION 
 #ifdef LAGRANGE_ORIENTATION_ACTIVE 
   /* if the particle is alive there is an extra velocity to add */

  #if defined(LAGRANGE_ORIENTATION_ACTIVE_JUMP)
   /* We perform jumps like a Copepod */

   #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP_TEMPERATURE
    #ifdef LB_TEMPERATURE
   /* the particle can react to the temperature value and decide if to jump */
   /* put temperature value in in the shear_rate */
   shear_rate = (tracer+ipart)->t;
    #endif
   #else
   /* the particle can react to flow gradients and decide if to jump */
   /* compute gamma_dot = sqrt( 2*S \ddot S) */
   (tracer+ipart)->shear_rate =  shear_rate = sqrt(2.)*sqrt(
   ((tracer+ipart)->dx_ux) *  ((tracer+ipart)->dx_ux) +
   ((tracer+ipart)->dy_uy) *  ((tracer+ipart)->dz_uy) +
   ((tracer+ipart)->dz_uz) *  ((tracer+ipart)->dz_uz) +
   0.5*( ((tracer+ipart)->dy_ux) + ((tracer+ipart)->dx_uy) )*( ((tracer+ipart)->dy_ux) + ((tracer+ipart)->dx_uy) ) +  
   0.5*( ((tracer+ipart)->dz_ux) + ((tracer+ipart)->dx_uz) )*( ((tracer+ipart)->dz_ux) + ((tracer+ipart)->dx_uz) ) +  
   0.5*( ((tracer+ipart)->dy_uz) + ((tracer+ipart)->dz_uy) )*( ((tracer+ipart)->dy_uz) + ((tracer+ipart)->dz_uy) ) );
   #endif

   jump_time_duration = - ((tracer+ipart)->jump_time)*log(0.01);

    if((tracer+ipart)->time_from_jump > jump_time_duration){
      /* set to zero the velocity amplitude */
      velocity_amplitude = 0.0;
      /* but we are in the condition to jump */

    #ifdef LAforGRANGE_ORIENTATION_ACTIVE_JUMP_HIGHPASS
	/* avoids low shear rate values */
      if( shear_rate  < (tracer+ipart)->critical_shear_rate )
    #else 
	/* avoids high shear rate values */	
      if( shear_rate  > (tracer+ipart)->critical_shear_rate )
    #endif

	{ /* begins one of the  if above */	

	/* reset the time from jump */ 	
	(tracer+ipart)->time_from_jump = 0.0;
        /* set initial velocity amplitude */
        velocity_amplitude = (tracer+ipart)->swim_velocity; 

     /* old version 
        // generarate a new random vector    
        vec = random_vector();
        // retrive last particle velocity 
        v_old.x = (tracer+ipart)->vx_old;
        v_old.y = (tracer+ipart)->vy_old;
        v_old.z = (tracer+ipart)->vz_old;
        // choose vec direction in the same emisphere of the particle velocity 
        if( scalar_product(vec,v_old) < 0.0 ) vec = vector_scale( -1.0 , vec);
        (tracer+ipart)->px = vec.x;
        (tracer+ipart)->py = vec.y;
        (tracer+ipart)->pz = vec.z;
     */
	/* new version : take for jump the present orientation */  
        (tracer+ipart)->px_jump = (tracer+ipart)->px;
        (tracer+ipart)->py_jump = (tracer+ipart)->py;
        (tracer+ipart)->pz_jump = (tracer+ipart)->pz;	
      } /* end if on shear rate */

    }else{
      /* hence (tracer+ipart)->time_from_jump < jump_time_duration  : we are already jumping */
      velocity_amplitude = (tracer+ipart)->swim_velocity * exp( -((tracer+ipart)->time_from_jump) / ((tracer+ipart)->jump_time) ); 
      
      if((tracer+ipart)->jump_time == 0.0) velocity_amplitude = 0.0;  /* this is for the tracer */
    }

      /* add jump part to the particle velocity */
      /* old version 
      (tracer+ipart)->vx += velocity_amplitude*((tracer+ipart)->px);
      (tracer+ipart)->vy += velocity_amplitude*((tracer+ipart)->py);
      (tracer+ipart)->vz += velocity_amplitude*((tracer+ipart)->pz);
      */
      /* new version : we jump with jeffrey */
      (tracer+ipart)->vx += velocity_amplitude*((tracer+ipart)->px_jump);
      (tracer+ipart)->vy += velocity_amplitude*((tracer+ipart)->py_jump);
      (tracer+ipart)->vz += velocity_amplitude*((tracer+ipart)->pz_jump);

      /* increase time from jump */
      (tracer+ipart)->time_from_jump += property.time_dt;

      //#endif /* LAGRANGE_ORIENTATION_ACTIVE_JUMP */
  #elif defined(LAGRANGE_ORIENTATION_ACTIVE_BALLISTIC)
      //#ifdef LAGRANGE_ORIENTATION_ACTIVE_BALLISTIC 
      /* Does not take the fluid velocity */
      /* it just swim with direction p evolved by jeffrey,random, etc. (if they are activated) */
   (tracer+ipart)->vx = ((tracer+ipart)->swim_velocity)*((tracer+ipart)->px);
   (tracer+ipart)->vy = ((tracer+ipart)->swim_velocity)*((tracer+ipart)->py);
   (tracer+ipart)->vz = ((tracer+ipart)->swim_velocity)*((tracer+ipart)->pz);
   //#endif /* LAGRANGE_ORIENTATION_ACTIVE_BALLISTIC */
  #elif defined(LAGRANGE_ORIENTATION_ACTIVE_TEMPERATURE)
   /* Kemotactic behaviour, take into account the temperature time increments to adjust swim velocity */
   /* add swimming velocity to fluid velocity */
   if( (tracer+ipart)->dt_t>0.0 ){ 
     reactivity = 1.0; 
   }else{ 
     reactivity = 0.0;
   }
   (tracer+ipart)->vx += reactivity * ((tracer+ipart)->swim_velocity)*((tracer+ipart)->px);
   (tracer+ipart)->vy += reactivity * ((tracer+ipart)->swim_velocity)*((tracer+ipart)->py);
   (tracer+ipart)->vz += reactivity * ((tracer+ipart)->swim_velocity)*((tracer+ipart)->pz);
  #else
   /* Default behaviour for LAGRANGE_ORIENTATION_ACTIVE : fluid velocity + swim with direction p evolved by jeffrey,random, etc. (if they are activated) */
   (tracer+ipart)->vx += ((tracer+ipart)->swim_velocity)*((tracer+ipart)->px);
   (tracer+ipart)->vy += ((tracer+ipart)->swim_velocity)*((tracer+ipart)->py);
   (tracer+ipart)->vz += ((tracer+ipart)->swim_velocity)*((tracer+ipart)->pz);
   // fprintf(stderr,"%g %g %g %g \n",(tracer+ipart)->name , (tracer+ipart)->vx , (tracer+ipart)->swim_velocity , (tracer+ipart)->px );
   //  fprintf(stderr,"%g %g %g %g \n",(tracer+ipart)->name , (tracer+ipart)->vy , (tracer+ipart)->swim_velocity , (tracer+ipart)->py );
  #endif

 #endif /* LAGRANGE_ORIENTATION_ACTIVE */
#endif /* LAGRANGE_ORIENTATION */

  if(itime==0 && resume==0){ 
    (tracer+ipart)->vx_old = (tracer+ipart)->vx;
    (tracer+ipart)->vy_old = (tracer+ipart)->vy;
    (tracer+ipart)->vz_old = (tracer+ipart)->vz;
  }
  /* Compute tracer acceleration : if v=u as here, then a = D_t u */
  (tracer+ipart)->ax = ((tracer+ipart)->vx - (tracer+ipart)->vx_old )/property.time_dt;
  (tracer+ipart)->ay = ((tracer+ipart)->vy - (tracer+ipart)->vy_old )/property.time_dt;
  (tracer+ipart)->az = ((tracer+ipart)->vz - (tracer+ipart)->vz_old )/property.time_dt;


   if(itime==0 && resume==0){
   /* Explicit Euler 1st order */
   (tracer+ipart)->x += property.time_dt*(tracer+ipart)->vx;
   (tracer+ipart)->y += property.time_dt*(tracer+ipart)->vy;
   (tracer+ipart)->z += property.time_dt*(tracer+ipart)->vz;
   }else{
   /* Adams-Bashforth 2nd order */
   (tracer+ipart)->x += property.time_dt*0.5*(3.0*(tracer+ipart)->vx - (tracer+ipart)->vx_old);
   (tracer+ipart)->y += property.time_dt*0.5*(3.0*(tracer+ipart)->vy - (tracer+ipart)->vy_old);
   (tracer+ipart)->z += property.time_dt*0.5*(3.0*(tracer+ipart)->vz - (tracer+ipart)->vz_old);
   }
   /* copy particle velocity in old */
   (tracer+ipart)->vx_old = (tracer+ipart)->vx; 
   (tracer+ipart)->vy_old = (tracer+ipart)->vy; 
   (tracer+ipart)->vz_old = (tracer+ipart)->vz;  

   /* copy fluid velocity in old */
   (tracer+ipart)->ux_old = (tracer+ipart)->ux; 
   (tracer+ipart)->uy_old = (tracer+ipart)->uy; 
   (tracer+ipart)->uz_old = (tracer+ipart)->uz; 

   }/* end if on fluid tracer */

   /* With drag force */ 
   //   if((tracer+ipart)->tau_drag != 0.0){
   if((tracer+ipart)->tau_drag > 0.0){    /* why >0 , because ==0 is a tracer and <0 is an eulerian probe */
  
   invtau = 1.0 / (tracer+ipart)->tau_drag;
   (tracer+ipart)->ax = ((tracer+ipart)->ux - (tracer+ipart)->vx)*invtau;
   (tracer+ipart)->ay = ((tracer+ipart)->uy - (tracer+ipart)->vy)*invtau;
   (tracer+ipart)->az = ((tracer+ipart)->uz - (tracer+ipart)->vz)*invtau;


#ifdef LAGRANGE_GRAVITY /* note that works only if LB_FORCING_GRAVITY is defined */
   /*  add: -g to acceleration */ 
 #ifdef LAGRANGE_GRAVITY_VARIABLE
     (tracer+ipart)->ax -= (tracer+ipart)->gravity_coeff*property.gravity_x;
     (tracer+ipart)->ay -= (tracer+ipart)->gravity_coeff*property.gravity_y;
     (tracer+ipart)->az -= (tracer+ipart)->gravity_coeff*property.gravity_z; 
 #else   
     (tracer+ipart)->ax -= property.gravity_x;
     (tracer+ipart)->ay -= property.gravity_y;
     (tracer+ipart)->az -= property.gravity_z; 
 #endif
#endif

#ifdef LAGRANGE_ADDEDMASS
  /* With Added mass */ 
   if((tracer+ipart)->beta_coeff != 0.0){

#ifdef LAGRANGE_GRAVITY
  /*  add also: -\beta*g to acceleration */ 
 #ifdef LAGRANGE_GRAVITY_VARIABLE
     (tracer+ipart)->ax -= (  - (tracer+ipart)->beta_coeff )*(tracer+ipart)->gravity_coeff*property.gravity_x;
     (tracer+ipart)->ay -= (  - (tracer+ipart)->beta_coeff )*(tracer+ipart)->gravity_coeff*property.gravity_y;
     (tracer+ipart)->az -= (  - (tracer+ipart)->beta_coeff )*(tracer+ipart)->gravity_coeff*property.gravity_z; 
 #else    
     (tracer+ipart)->ax -= (  - (tracer+ipart)->beta_coeff )*property.gravity_x;
     (tracer+ipart)->ay -= (  - (tracer+ipart)->beta_coeff )*property.gravity_y;
     (tracer+ipart)->az -= (  - (tracer+ipart)->beta_coeff )*property.gravity_z; 
 #endif    
#endif

  if(itime==0 && resume==0){ 
    (tracer+ipart)->ux_old = (tracer+ipart)->ux;
    (tracer+ipart)->uy_old = (tracer+ipart)->uy;
    (tracer+ipart)->uz_old = (tracer+ipart)->uz;
  }

   /* Here I will write the computation of the fluid material derivative */
   Dt_u.x = ((tracer+ipart)->ux - (tracer+ipart)->ux_old )/property.time_dt
          + ((tracer+ipart)->ux - (tracer+ipart)->vx)*(tracer+ipart)->dx_ux 
          + ((tracer+ipart)->uy - (tracer+ipart)->vy)*(tracer+ipart)->dy_ux 
          + ((tracer+ipart)->uz - (tracer+ipart)->vz)*(tracer+ipart)->dz_ux;

   Dt_u.y = ((tracer+ipart)->uy - (tracer+ipart)->uy_old )/property.time_dt
          + ((tracer+ipart)->ux - (tracer+ipart)->vx)*(tracer+ipart)->dx_uy 
          + ((tracer+ipart)->uy - (tracer+ipart)->vy)*(tracer+ipart)->dy_uy 
          + ((tracer+ipart)->uz - (tracer+ipart)->vz)*(tracer+ipart)->dz_uy;

   Dt_u.z = ((tracer+ipart)->uz - (tracer+ipart)->uz_old )/property.time_dt         
	  + ((tracer+ipart)->ux - (tracer+ipart)->vx)*(tracer+ipart)->dx_uz 
          + ((tracer+ipart)->uy - (tracer+ipart)->vy)*(tracer+ipart)->dy_uz 
          + ((tracer+ipart)->uz - (tracer+ipart)->vz)*(tracer+ipart)->dz_uz;

   (tracer+ipart)->ax += Dt_u.x*(tracer+ipart)->beta_coeff;
   (tracer+ipart)->ay += Dt_u.y*(tracer+ipart)->beta_coeff;
   (tracer+ipart)->az += Dt_u.z*(tracer+ipart)->beta_coeff;


 #ifdef LAGRANGE_ADDEDMASS_LIFT
   /* Here we add the Lift force */

   /* compute vorticity vector : omega = nabla x u */
   omega.x =  (tracer+ipart)->dy_uz  - (tracer+ipart)->dz_uy;
   omega.y =  (tracer+ipart)->dz_ux  - (tracer+ipart)->dx_uz;
   omega.z =  (tracer+ipart)->dx_uy  - (tracer+ipart)->dy_ux;
   
   /* lift force computation assuming lift coefficient CL = 1/2 */
   /*  - beta/3 (u - v) x omega */
   lift_coeff = ( (tracer+ipart)->beta_coeff )/3.0;

   (tracer+ipart)->ax += -lift_coeff*( ((tracer+ipart)->uy - (tracer+ipart)->vy))*omega.z - ((tracer+ipart)->uz - (tracer+ipart)->vz))*omega.y );
   (tracer+ipart)->ay += -lift_coeff*( ((tracer+ipart)->uz - (tracer+ipart)->vz))*omega.x - ((tracer+ipart)->ux - (tracer+ipart)->vx))*omega.z );
   (tracer+ipart)->az += -lift_coeff*( ((tracer+ipart)->ux - (tracer+ipart)->vx))*omega.y - ((tracer+ipart)->uy - (tracer+ipart)->vy))*omega.x );

 #endif /* end of lift */

   }/* end of if on addedd mass */
#endif

   if(itime==0 && resume==0){
     (tracer+ipart)->vx += property.time_dt*(tracer+ipart)->ax;
     (tracer+ipart)->vy += property.time_dt*(tracer+ipart)->ay;
     (tracer+ipart)->vz += property.time_dt*(tracer+ipart)->az;
   }else{
     (tracer+ipart)->vx += property.time_dt*0.5*(3.0*(tracer+ipart)->ax - (tracer+ipart)->ax_old);
     (tracer+ipart)->vy += property.time_dt*0.5*(3.0*(tracer+ipart)->ay - (tracer+ipart)->ay_old);
     (tracer+ipart)->vz += property.time_dt*0.5*(3.0*(tracer+ipart)->az - (tracer+ipart)->az_old);
   }

   (tracer+ipart)->ax_old = (tracer+ipart)->ax; 
   (tracer+ipart)->ay_old = (tracer+ipart)->ay; 
   (tracer+ipart)->az_old = (tracer+ipart)->az; 

   if(itime==0 && resume==0){
     (tracer+ipart)->x += property.time_dt*(tracer+ipart)->vx; 
     (tracer+ipart)->y += property.time_dt*(tracer+ipart)->vy;
     (tracer+ipart)->z += property.time_dt*(tracer+ipart)->vz;
   }else{
     (tracer+ipart)->x += property.time_dt*0.5*(3.0*(tracer+ipart)->vx - (tracer+ipart)->vx_old);
     (tracer+ipart)->y += property.time_dt*0.5*(3.0*(tracer+ipart)->vy - (tracer+ipart)->vy_old);
     (tracer+ipart)->z += property.time_dt*0.5*(3.0*(tracer+ipart)->vz - (tracer+ipart)->vz_old);
   }
   (tracer+ipart)->vx_old = (tracer+ipart)->vx; 
   (tracer+ipart)->vy_old = (tracer+ipart)->vy; 
   (tracer+ipart)->vz_old = (tracer+ipart)->vz; 

   (tracer+ipart)->ux_old = (tracer+ipart)->ux; 
   (tracer+ipart)->uy_old = (tracer+ipart)->uy; 
   (tracer+ipart)->uz_old = (tracer+ipart)->uz; 

   }/* end of if on tau_drag different from zero */

#ifdef LAGRANGE_ORIENTATION

	      /* assign P vector */
     vecP[0] = (tracer+ipart)->px;
     vecP[1] = (tracer+ipart)->py;
     vecP[2] = (tracer+ipart)->pz;
              /* assign the last dP /dt  vector */
     vecFold[0] = (tracer+ipart)->dt_px;
     vecFold[1] = (tracer+ipart)->dt_py;
     vecFold[2] = (tracer+ipart)->dt_pz;

 #ifdef LAGRANGE_ORIENTATION_SECONDORIENTATION
	      /* assign N normal vector */
     vecN[0] = (tracer+ipart)->nx;
     vecN[1] = (tracer+ipart)->ny;
     vecN[2] = (tracer+ipart)->nz;
              /* assign the last dN /dt  vector */
     vecFNold[0] = (tracer+ipart)->dt_nx;
     vecFNold[1] = (tracer+ipart)->dt_ny;
     vecFNold[2] = (tracer+ipart)->dt_nz;
 #endif

 #ifdef LAGRANGE_ORIENTATION_JEFFREY
 /* Here we implement Jeffrey equation */

              /* aspect ratio factor */
   alpha = (tracer+ipart)->aspect_ratio;
   f_alpha = (alpha*alpha-1.0)/(1.0+alpha*alpha);


	      /* velocity gradient matrix */
	      matA[0][0]=(tracer+ipart)->dx_ux ; matA[0][1]=(tracer+ipart)->dy_ux; matA[0][2]=(tracer+ipart)->dz_ux;
	      matA[1][0]=(tracer+ipart)->dx_uy ; matA[1][1]=(tracer+ipart)->dy_uy; matA[1][2]=(tracer+ipart)->dz_uy;
	      matA[2][0]=(tracer+ipart)->dx_uz ; matA[2][1]=(tracer+ipart)->dy_uz; matA[2][2]=(tracer+ipart)->dz_uz;	 

	      /* Compute Sij */
	      /*make simmetric*/
	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matS[i][j] = 0.5*(matA[i][j]+matA[j][i]);
		}

	      /* Compute Wij */
	      /*make simmetric*/
	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matW[i][j] = 0.5*(matA[i][j]-matA[j][i]);
		}

	      /* multiply S by the aspect ratio stretching factor */
	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matS[i][j] *= f_alpha;
		}

  #ifdef LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS
	       /* gravitational gyrotaxis : the stretched S matrix has an extra term -1/(2*v0) * g_i p_j      */	      
	      if((tracer+ipart)->gyrotaxis_velocity !=0){
		gyro = - 0.5 / (tracer+ipart)->gyrotaxis_velocity;
   #ifdef LAGRANGE_GRAVITY	      
		/* the full term */	      
		vecA[0] = gyro * (- property.gravity_x  + (tracer+ipart)->ax);
		vecA[1] = gyro * (- property.gravity_y  + (tracer+ipart)->ay);
		vecA[2] = gyro * (- property.gravity_z  + (tracer+ipart)->az);
	      
   #else
		/* just for a test: fixed vector along z. Like g_x=0, g_y=0, g_z=-1.0 */
		vecA[0] = gyro * 0.0;  
		vecA[1] = gyro * 0.0;
		vecA[2] = gyro * 1.0;
		/* only acceleration */
		/*
		  vecA[0] = gyro * (tracer+ipart)->ax;
		  vecA[1] = gyro * (tracer+ipart)->ay;
		  vecA[2] = gyro * (tracer+ipart)->az;
		*/
   #endif	      
		matS[0][0]+= vecA[0]*vecP[0]; matA[0][1] += vecA[0]*vecP[1]; matA[0][2] += vecA[0]*vecP[2];
		matS[1][0]+= vecA[1]*vecP[0]; matA[1][1] += vecA[1]*vecP[1]; matA[1][2] += vecA[1]*vecP[2];
		matS[2][0]+= vecA[2]*vecP[0]; matA[2][1] += vecA[2]*vecP[1]; matA[2][2] += vecA[2]*vecP[2];	      
	      }/* end if on gyrotaxis_velocity !=0 , to avoid nan */
  #endif /* LAGRANGE_ORIENTATION_JEFFREY_GYROTAXIS */

	      /* Now we compute RHS of the Jeffrey equation */

	      /* first product TMP[i] = S[i][j]*P[j] */
	      vecTMP[0]=vecTMP[1]=vecTMP[2]=0;
	      for (i=0; i<3; i++)
		for (j=0; j<3; j++){
		  vecTMP[i] += matS[i][j]*vecP[j];
		}
	      /* then the product of the two vectors P[i]*TMP[i] to form a scalar */
	      scalOSO = 0.0;	      
	      for (i=0; i<3; i++){
		  scalOSO += vecP[i]*vecTMP[i];
		}

	    /* here I add on all the contributions */
	      for (i=0; i<3; i++){
		  vecF[i] = 0.0; 
		for (j=0; j<3; j++){	
		  vecF[i] += matW[i][j]*vecP[j] + ( matS[i][j]*vecP[j]);
		}
		  vecF[i] -=  vecP[i]*scalOSO;
	      }

 #ifdef LAGRANGE_ORIENTATION_SECONDORIENTATION
 	      /* Now we compute RHS of the Jeffrey equation for N (it is the vector normal to P) */
	      /* first product TMP[i] = S[i][j]*N[j] */
	      vecTMP[0]=vecTMP[1]=vecTMP[2]=0;
	      for (i=0; i<3; i++)
		for (j=0; j<3; j++){
                  vecTMP[i] += matS[i][j]*vecN[j];  /* <----- NOTE differently from before here we have the N vector */
		}
	      /* then the product of the two vectors P[i]*TMP[i] to form a scalar */
	      scalOSO = 0.0;	      
	      for (i=0; i<3; i++){
		  scalOSO += vecP[i]*vecTMP[i];
		}

	    /* here I add on all the contributions */
	      for (i=0; i<3; i++){
		  vecFN[i] = 0.0; 
		for (j=0; j<3; j++){	
		  vecFN[i] += matW[i][j]*vecN[j];
		}
		  vecFN[i] -=  vecP[i]*scalOSO;
	      }

 #endif
	     	      
   /* if restart Euler 1st order G = G0 + (DT)*F  */
 if(itime==0 && resume==0){
		  for (i=0; i<3; i++){	
		    vecP[i] =  vecP[i] + property.time_dt*vecF[i];	  
		    /* copy the old term */
		    vecFold[i]  = vecF[i];
		}
 }else{
   /* AB 2nd order G = G0 + (DT/2)*(3*F - Fold)  */
	      for (i=0; i<3; i++){		 
		  vecP[i] =  vecP[i] + 0.5*property.time_dt*( 3.*vecF[i] -  vecFold[i] );
		  /* copy the old term */
		  vecFold[i]  = vecF[i];
		}
 }

 #ifdef LAGRANGE_ORIENTATION_SECONDORIENTATION

   /* if restart Euler 1st order G = G0 + (DT)*F  */
 if(itime==0 && resume==0){
		  for (i=0; i<3; i++){	
		    vecN[i] =  vecN[i] + property.time_dt*vecFN[i];	  
		    /* copy the old term */
		    vecFNold[i]  = vecFN[i];
		}
 }else{
   /* AB 2nd order G = G0 + (DT/2)*(3*F - Fold)  */
	      for (i=0; i<3; i++){		 
		  vecN[i] =  vecN[i] + 0.5*property.time_dt*( 3.*vecFN[i] -  vecFNold[i] );
		  /* copy the old term */
		  vecFNold[i]  = vecFN[i];
		}
 }

 #endif

 #endif /*LAGRANGE_ORIENTATION_JEFFREY */
 #ifdef LAGRANGE_ORIENTATION_DIFFUSION
	     vec_xi[0] =  random_gauss(0.0,1.0);
	     vec_xi[1] =  random_gauss(0.0,1.0);
	     vec_xi[2] =  random_gauss(0.0,1.0); 	     

	     vecTMP[0]=vecTMP[1]=vecTMP[2]=0;
	      for (i=0; i<3; i++)
		for (j=0; j<3; j++){
		  vecTMP[i] += ( matI[i][j]- vecP[i]*vecP[j] )*vec_xi[j];
		}

	      two_d_r = 2.0*(tracer+ipart)->rotational_diffusion;
	      /* stochastic part of the rotation equation */
	      for (i=0; i<3; i++)
 	      vecP[i] +=  sqrt(two_d_r*property.time_dt)*vecTMP[i];
 #endif /* LAGRANGE_ORIENTATION_DIFFUSION */

 #ifdef LAGRANGE_ORIENTATION_RANDOM
 /* randomly oriented vector */
 vec = random_vector();
#ifdef GRID_POP_D2Q9
 /* generate a random vector in the x,y plane */
 vec = random_vector_2d();
#endif
 vecP[0] = vec.x;
 vecP[1] = vec.y;
 vecP[2] = vec.z;

 /* compute acceleration */
 vecFold[0] = (vecP[0] - (tracer+ipart)->px)/property.time_dt;
 vecFold[1] = (vecP[1] - (tracer+ipart)->py)/property.time_dt;
 vecFold[2] = (vecP[2] - (tracer+ipart)->pz)/property.time_dt;
 #endif /* LAGRANGE_ORIENTATION_RANDOM */

	      /* normalize P vector */
	      norm=0.0;
	      for (i=0; i<3; i++) norm += vecP[i]*vecP[i];
	      for (i=0; i<3; i++) vecP[i]/=sqrt(norm);

	      /* assign P vector */
     (tracer+ipart)->px = vecP[0];
     (tracer+ipart)->py = vecP[1];
     (tracer+ipart)->pz = vecP[2];

              /* assign the just computed dP /dt  vector */     
     (tracer+ipart)->dt_px = vecFold[0];
     (tracer+ipart)->dt_py = vecFold[1];
     (tracer+ipart)->dt_pz = vecFold[2];

  #ifdef LAGRANGE_ORIENTATION_SECONDORIENTATION

	      /* normalize N vector */
	      norm=0.0;
	      for (i=0; i<3; i++) norm += vecN[i]*vecN[i];
	      for (i=0; i<3; i++) vecN[i]/=sqrt(norm);

	      /* assign N vector */
     (tracer+ipart)->nx = vecN[0];
     (tracer+ipart)->ny = vecN[1];
     (tracer+ipart)->nz = vecN[2];

              /* assign the just computed dN /dt  vector */     
     (tracer+ipart)->dt_nx = vecFNold[0];
     (tracer+ipart)->dt_ny = vecFNold[1];
     (tracer+ipart)->dt_nz = vecFNold[2];
  #endif
     
#endif /* end of lagrange orientation */

   /* In case of BC we use elastic bouncing rule for the particle */
  
#ifdef LB_FLUID_BC
 #ifdef LB_FLUID_BC_Y

   if( (tracer+ipart)->y < 0.0 ){
     (tracer+ipart)->y *= -1.0; 
     //(tracer+ipart)->vx *= -1.0;
     if( (tracer+ipart)->vy < 0.0 ) (tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
  #ifdef LAGRANGE_ORIENTATION_BC
     if( (tracer+ipart)->py < 0.0 ) (tracer+ipart)->py *= -1.0;
  #endif 
   }

   if( (tracer+ipart)->y >= property.SY ){
     (tracer+ipart)->y = property.SY- ( (tracer+ipart)->y-property.SY ); 
     //(tracer+ipart)->vx *= -1.0;
     if( (tracer+ipart)->vy > 0.0 ) (tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
  #ifdef LAGRANGE_ORIENTATION_BC
     if( (tracer+ipart)->py > 0.0 ) (tracer+ipart)->py *= -1.0;
  #endif 
   }
 #endif

 #ifdef LB_FLUID_BC_X

   if( (tracer+ipart)->x < 0.0 ){
     (tracer+ipart)->x *= -1.0; 
     (tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
   }

   if( (tracer+ipart)->x >= property.SX ){
     (tracer+ipart)->x = property.SX- ( (tracer+ipart)->x-property.SX ); 
     (tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     //(tracer+ipart)->vz *= -1.0;
   }
 #endif

 #ifdef LB_FLUID_BC_Z

   if( (tracer+ipart)->z < 0.0 ){
     (tracer+ipart)->z *= -1.0; 
     //(tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     (tracer+ipart)->vz *= -1.0;
   }

   if( (tracer+ipart)->z >= property.SZ ){
     (tracer+ipart)->z = property.SZ- ( (tracer+ipart)->z-property.SZ ); 
     //(tracer+ipart)->vx *= -1.0;
     //(tracer+ipart)->vy *= -1.0;
     (tracer+ipart)->vz *= -1.0;
   }
 #endif

#endif /* endof  LB_FLUID_BC */
#ifdef LAGRANGE_GRADIENT
 #ifdef LAGRANGE_POLYMER
  evolve_lagrangian_polymer_conformation_tensor(ipart);
 #endif
#endif   

}/* end of loop on particles */

  sendrecv_particles();

}/* end of move_particles */

