
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

    #ifdef LAGRANGE_ORIENTATION_ACTIVE_JUMP_HIGHPASS
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
  
#ifndef LAGRANGE_ORIENTATION_DRAG
   /* Stokes drag acceleration on a sphere */
   invtau = 1.0 / (tracer+ipart)->tau_drag; 
   (tracer+ipart)->ax = ((tracer+ipart)->ux - (tracer+ipart)->vx)*invtau;
   (tracer+ipart)->ay = ((tracer+ipart)->uy - (tracer+ipart)->vy)*invtau;
   (tracer+ipart)->az = ((tracer+ipart)->uz - (tracer+ipart)->vz)*invtau;
#else
   /* Stokes drag force on a NON spherical (axisymmetric ellipsoidal) particle */
   /* same notation as in PRL 119, 254501 (2017) see also its supplementary materials */
   /* The aspect ratio is a_par / a_perp */ 
    alpha = (tracer+ipart)->aspect_ratio; 
   /* Note that the meaning of the parameter tau_drag in our code is  ( \rho_p 2 a_perp^2 ) / (9 \nu \rho_f )  */
   /* The particle relaxation time however is ~ ( \rho_p 2 a_perp a_par ) / (9 \nu \rho_f ) = tau_drag * alpha */
    invtau = 1.0 / ( (tracer+ipart)->tau_drag * alpha );  
   /*compute prefactors*/
   if (alpha == 1){ 
     beta = c_perp = c_par  = 1.0;
   }else{
    if (alpha >  1)  beta = log(alpha + sqrt(alpha*alpha - 1.0)) / ( alpha*sqrt(alpha*alpha - 1.0) );
    if (alpha <  1)  beta = acos(alpha) / ( alpha*sqrt(1.0 - alpha*alpha) );
    c_perp = 8.0*(alpha*alpha - 1.) / ( 3.0*alpha*((2.0*alpha*alpha - 3.0)*beta + 1.0)  );
    c_par  = 4.0*(alpha*alpha - 1.) / ( 3.0*alpha*((2.0*alpha*alpha - 1.0)*beta - 1.0)  );
   }
   //fprintf(stderr,"%e %e %e %e\n",alpha,beta,c_perp,c_par);
   /* assign P vector */
   vecP[0] = (tracer+ipart)->px;
   vecP[1] = (tracer+ipart)->py;
   vecP[2] = (tracer+ipart)->pz;
   /* build drag tensor*/
   matM[0][0] =     c_perp + (c_par - c_perp)*vecP[0]*vecP[0];
   matM[1][1] =     c_perp + (c_par - c_perp)*vecP[1]*vecP[1];
   matM[2][2] =     c_perp + (c_par - c_perp)*vecP[2]*vecP[2];
   matM[0][1] = matM[1][0] = (c_par - c_perp)*vecP[0]*vecP[1];
   matM[0][2] = matM[2][0] = (c_par - c_perp)*vecP[0]*vecP[2];
   matM[1][2] = matM[2][1] = (c_par - c_perp)*vecP[1]*vecP[2];
   /* velocity difference */
   uvx = (tracer+ipart)->ux - (tracer+ipart)->vx;
   uvy = (tracer+ipart)->uy - (tracer+ipart)->vy;
   uvz = (tracer+ipart)->uz - (tracer+ipart)->vz;
   /* acceleration */
   (tracer+ipart)->ax = (matM[0][0]*uvx +  matM[0][1]*uvy +  matM[0][2]*uvz)*invtau;
   (tracer+ipart)->ay = (matM[1][0]*uvx +  matM[1][1]*uvy +  matM[1][2]*uvz)*invtau;
   (tracer+ipart)->az = (matM[2][0]*uvx +  matM[2][1]*uvy +  matM[2][2]*uvz)*invtau;    
#endif





/* --------------------------------------------------------------------------------------------------------------*/
/* inserting the program for translation and rotation using quaternions */ 
/* --------------------------------------------------------------------------------------------------------------*/




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



 #endif /*LAGRANGE_ORIENTATION_JEFFREY */





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

     
#endif /* end of lagrange orientation */

   /* In case of BC we use elastic bouncing rule for the particle */
  
#ifdef LB_FLUID_BC

#endif /* endof  LB_FLUID_BC */
 

}/* end of loop on particles */

  sendrecv_particles();

}/* end of move_particles */


