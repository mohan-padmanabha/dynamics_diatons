
      /* translation */
    my_double q0, q1, q2, q3;
    my_double matQ[2][2];
    my_double aspr, mp, a, rhop, mu, pi;
    my_double k_xx, k_yy, k_zz;
    my_double k_tranx, k_trany, k_tranz;
    my_double fx, fy, fz;
    my_double v12x, v12y, v12z;
    


/* orinetation */
          my_double I_xx, I_yy, I_zz;
          my_double q0, q1, q2, q3;
          my_double matD[3][3], matA[3][3], matS[3][3], matDT[3][3], matW[3][3], matB[3][3] ;
          my_double alpha, beta, gamma;
          my_double T_x, T_y, T_z;
          my_double do_x, do_y, do_z;
          my_double o12x, o12y, o12z;
          my_double obx, oby, obz;
          my_double dq0dt, dq1dt, dq2dt, dq3dt;
          my_double qsq, invq, p, r;
          int i, j, k;

   

   pi = 3.1415;
    
    



    /* assigning to the local variable */
    q0 = (tracer+ipart)->q[0];
    q1 = (tracer+ipart)->q[1];
    q2 = (tracer+ipart)->q[2];
    q3 = (tracer+ipart)->q[3];





/* to initilize the values of the quaternions */


      quaternion euler2quaternion(euler e)
      {
        quaternion q;
        q.q[0] = cos(0.5*e.theta)*cos(0.5*(e.phi+e.psi));
        q.q[1] = sin(0.5*e.theta)*cos(0.5*(e.phi-e.psi));
        q.q[2] = sin(0.5*e.theta)*sin(0.5*(e.phi-e.psi));
        q.q[3] = cos(0.5*e.theta)*sin(0.5*(e.phi+e.psi));
        return q;
      }

     /* rotation matrix of quaternions */
    matQ[0][0] = 1-2(q2*q2+q3*q3);  matQ[1][0] = 2*(q1*q2-q0*q3);   matQ[2][0] = 2*(q1*q3+q0*q2); 
    matQ[0][1] = 2*(q1*q2+q0*q3);   matQ[1][1] = 1-2(q1*q1+q3*q3);  matQ[2][1] = 2*(q2*q3-q0*q1);
    matQ[0][2] = 2*(q1*q3-q0*q2);   matQ[1][2] = 2*(q2*q3+q0*q1);   matQ[2][2] = 1-2(q1*q1+q2*q2);  
 
     

 // aspect ratio is major axis / minor axis
 // already present in the code 
//    (tracer+ipart)->aspect_ratio = b/a ;
  aspr = (tracer+ipart)->aspect_ratio ; 
// mp is the mass of the particle (zhang et al 2001)
    mp = 4/3*pi*pow(a,3)*aspr*rhop ;

    K_xx = (16*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr*aspr - 3)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));

    K_yy = K_xx;

    K_zz = (8*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr*aspr - 1)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));


matK[0][0]= kxx;
matK[1][1]= kyy;
matK[2][2]= kzz;


/* translation dyadic or resistance tensor */

          /* For transformation of co-ordinates from Lab frome to body frame */
          /* general form  */
               for (k=1; k<3; k++){
                     p = 0;
                 for (j=0; j<3; j++){
                     r = 0;
                   for (i=0; i<3; i++){
                     r += matQ[j][i]*matK[i][j];
                    }
                     matKT[k][j] += r * matA[k][j];
                      }
                         }
                       
            
                       

//    K_tranx = (matQ[0][0] * K_xx)* matQ[0][0]+ (matQ[0][1] * K_xx)* matQ[0][1]+(matQ[0][2] * K_xx)* matQ[0][2];

//    K_trany = (matQ[1][0] * K_yy )* matQ[1][0]+ ( matQ[1][1] * K_yy)* matQ[1][1]+(matQ[1][2] * K_yy)* matQ[1][2];

//    K_tranz = (matQ[2][0] * K_zz)* matQ[2][0]+ ( matQ[2][1] * K_zz)* matQ[2][1]+(matQ[2][2] * K_zz)* matQ[2][2];

/* hydrodynamic drag force (brenner 1964 && zhang et al 2001) */

//    fx = mu*K_tranx * ((tracer+ipart)->ux - (tracer+ipart)->vx);
//    fy = mu*K_trany * ((tracer+ipart)->uy - (tracer+ipart)->vy);
//    fz = mu*K_tranz * ((tracer+ipart)->uz - (tracer+ipart)->vz);   

    fx = property.nu * matKT[0][0] * ((tracer+ipart)->ux - (tracer+ipart)->vx);
    fy = property.nu * matKT[1][1] * ((tracer+ipart)->uy - (tracer+ipart)->vy);
    fz = property.nu * matKT[2][2] * ((tracer+ipart)->uz - (tracer+ipart)->vz); 



    /* Here compute v particle at time n-1/2 */
    v12x =  0.5 * ((tracer+ipart)->vx + (tracer+ipart)->vx0);
    v12y =  0.5 * ((tracer+ipart)->vy + (tracer+ipart)->vy0);
    v12z =  0.5 * ((tracer+ipart)->vz + (tracer+ipart)->vz0);

    /* here we assign the present time value to previous variable */
    (tracer+ipart)->vx0 = (tracer+ipart)->vx;
    (tracer+ipart)->vy0 = (tracer+ipart)->vy;
    (tracer+ipart)->vz0 = (tracer+ipart)->vz;

    /* Here compute v particle at time n+1/2 */

    (tracer+ipart)->vx = v12x + (fx * property.time_dt) / mp;
    (tracer+ipart)->vy = v12y + (fy * property.time_dt) / mp;
    (tracer+ipart)->vz = v12z + (fz * property.time_dt) / mp;

    /* Here compute v particle at time n */

      /* computation of the position */
    (tracer+ipart)->x += ((tracer+ipart)->vx * property.time_dt);
    (tracer+ipart)->y += ((tracer+ipart)->vy * property.time_dt);
    (tracer+ipart)->z += ((tracer+ipart)->vz * property.time_dt);
 
    /* to find the value of v  by extrapolation at time n +1 for next iteration */ 
    (tracer+ipart)->vx = 0.5 * (v12x + 3*(tracer+ipart)->vx);
    (tracer+ipart)->vy = 0.5 * (v12y + 3*(tracer+ipart)->vy);
    (tracer+ipart)->vz = 0.5 * (v12z + 3*(tracer+ipart)->vz);   




    /* multiplication matrix of quaternions to advance in time */
    B[0][0] =  q0;    B[0][1] = -q1;    B[0][2] = -q2;    B[0][3] = -q3;    
    B[1][0] =  q1;    B[1][1] =  q0;    B[1][2] = -q3;    B[1][3] =  q2;
    B[2][0] =  q2;    B[2][1] =  q3;    B[2][2] =  q0;    B[2][3] = -q1;
    B[3][0] =  q3;    B[3][1] = -q2;    B[3][2] =  q1;    B[3][3] =  q0;



              /* moment of inertia (diagnal in body frame)*/
              I_xx = ((1+(aspr*aspr))*(a*a)*mp)/(5) ;
              I_yy = I_xx ;
              I_zz = (2*a*a*mp)/(5) ;


	      /* velocity gradient matrix */
	      matD[0][0]=(tracer+ipart)->dx_ux ; matD[0][1]=(tracer+ipart)->dy_ux; matD[0][2]=(tracer+ipart)->dz_ux;
	      matD[1][0]=(tracer+ipart)->dx_uy ; matD[1][1]=(tracer+ipart)->dy_uy; matD[1][2]=(tracer+ipart)->dz_uy;
	      matD[2][0]=(tracer+ipart)->dx_uz ; matD[2][1]=(tracer+ipart)->dy_uz; matD[2][2]=(tracer+ipart)->dz_uz;	 

          /* For transformation of co-ordinates from Lab frome to body frame for tensor */
               for (k=1; k<3; k++){
                     p = 0;
                 for (j=0; j<3; j++){
                     r = 0;
                   for (i=0; i<3; i++){
                     r += matA[j][i]*matD[i][j];
                    }
                     matDT[k][j] += r * matA[k][j];
                      }
                         }
      

    	      /* Compute Sij for body frame*/
	      /*make simmetric*/

	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matS[i][j] = 0.5*(matDT[i][j]+matDT[j][i]);
		}

	      /* Compute Wij for body frame*/
	      /*make simmetric*/

	      for (i=0; i<3; i++)
                for (j=0; j<3; j++){
                  matW[i][j] = 0.5*(matDT[i][j]-matDT[j][i]);
		}

              /* computation of constants alpha beta and gamma to compute torque */
    alpha = ((aspr*aspr)/((aspr*aspr) - 1)) + aspr/(2*pow((aspr*aspr)-1),3/2))* ln ((aspr - sqrt((aspr*aspr) - 1))/(aspr + sqrt((aspr*aspr)-1)));
    beta = alpha ;
    gamma = - (2/((aspr*aspr) - 1)) - aspr/(pow((aspr*aspr)-1),3/2))* ln ((aspr - sqrt((aspr*aspr) - 1))/(aspr + sqrt((aspr*aspr)-1))) ;


        /* computation of torque using the values computed */
    T_x = (16*pi*property.nu*pow(a,3)*aspr)/(3*(beta + (aspr*aspr)*gamma)) * ((1-(aspr*aspr)) * matS[2][1] + (1+ (aspr*aspr)) * ( matW[1][2] - (tracer+ipart)->ox)) ;

    T_y = (16*pi*property.nu*pow(a,3)*aspr)/(3*((aspr*aspr)*gamma + alpha))*(((aspr*aspr) - 1) * matS[0][2] + ((aspr*aspr)+1) * ( matW[0][2] - (tracer+ipart)->oy)) ;

    T_z = (32*pi*property.nu*pow(a,3)*aspr)/(3*(alpha + gamma))*( matW[1][0] - (tracer+ipart)->oz) ;
    

    /* compute angular velocity */
 //   obx = lbx/(tracer+ipart)->Ixx;
 //  oby = lby/(tracer+ipart)->Iyy;
 //   obz = lbz/(tracer+ipart)->Izz;

    do_x = (T_x + (tracer+ipart)->oy*(tracer+ipart)->oz*(I_yy - I_zz))/I_xx;
    do_y = (T_y + (tracer+ipart)->oz* (tracer+ipart)->ox *(I_zz - I_xx))/I_yy;
    
    /* I_xx = I_yy hence the equation is simplified */
    do_z = (T_z)/I_zz;


    /* taking average to determine the middle value */
    o12x = 0.5((tracer+ipart)->ox0+(tracer+ipart)->ox);
    o12y = 0.5((tracer+ipart)->oy0+(tracer+ipart)->oy);
    o12z = 0.5((tracer+ipart)->oz0+(tracer+ipart)->oz);


   /*assigning the old values*/
    (tracer+ipart)->ox0=(tracer+ipart)->ox;
    (tracer+ipart)->oy0=(tracer+ipart)->oy;
    (tracer+ipart)->oz0=(tracer+ipart)->oz;

    /* estimating the value of omega at time t */
    (tracer+ipart)->ox = o12x + (do_x * (property.time_dt*0.5));
    (tracer+ipart)->oy = o12y + (do_y * (property.time_dt*0.5));
    (tracer+ipart)->oz = o12z + (do_z * (property.time_dt*0.5));   





    /* transformation: From angular momentum in the lab frame to body frame */
    obx = A[0][0] *(tracer+ipart)->ox + A[0][1] *(tracer+ipart)->oy + A[0][2] *(tracer+ipart)->oz;
    oby = A[1][0] *(tracer+ipart)->ox + A[1][1] *(tracer+ipart)->oy + A[1][2] *(tracer+ipart)->oz;
    obz = A[2][0] *(tracer+ipart)->ox + A[2][1] *(tracer+ipart)->oy + A[2][2] *(tracer+ipart)->oz;

   /* Now compute dq/dt = 0.5* B*(0,obx,oby,obz) , this is just a guess of dq/dt(t) */
    dq0dt = 0.5 *  (B[0][1] * obx + B[0][2] * oby + B[0][3] * obz );
    dq1dt = 0.5 *  (B[1][1] * obx + B[1][2] * oby + B[1][3] * obz );
    dq2dt = 0.5 *  (B[2][1] * obx + B[2][2] * oby + B[2][3] * obz );
    dq3dt = 0.5 *  (B[3][1] * obx + B[3][2] * oby + B[3][3] * obz );

    
    /* advancement in time  q(t+0.5*time_dt) = q(t) + 0.5*time_dt*dqdt(t) */
    q0 += 0.5*property.time_dt*dq0dt;
    q1 += 0.5*property.time_dt*dq1dt;
    q2 += 0.5*property.time_dt*dq2dt;
    q3 += 0.5*property.time_dt*dq3dt;

    
    /* quaternion normalization */
    
    qsq =q0*q0 + q1*q1 + q2*q2 + q3*q3;
    
    invq = 1.0/ sqrt(qsq);

    q0 *= invq;
    q1 *= invq;
    q2 *= invq;
    q3 *= invq;

    /* estimating the value of omega at time t+dt/2  */
    (tracer+ipart)->ox = o12x + do_x *property.time_dt;
    (tracer+ipart)->oy = o12y + do_y *property.time_dt;
    (tracer+ipart)->oz = o12z + do_z *property.time_dt;


/* Rotation matrix of quaternions to transform  to body frame */
    A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    A[0][1] = 2*(q1*q2+q0*q3);
    A[0][2] = 2*(q1*q3-q0*q2);
    
    A[1][0] = 2*(q1*q2-q0*q3);
    A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A[1][2] = 2*(q2*q3+q0*q1);
    
    A[2][0] = 2*(q1*q3+q0*q2);
    A[2][1] = 2*(q2*q3-q0*q1);
    A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;


    /* multiplication matrix of quaternions to advance in time */
    B[0][0] =  q0;    B[0][1] = -q1;    B[0][2] = -q2;    B[0][3] = -q3;    
    B[1][0] =  q1;    B[1][1] =  q0;    B[1][2] = -q3;    B[1][3] =  q2;
    B[2][0] =  q2;    B[2][1] =  q3;    B[2][2] =  q0;    B[2][3] = -q1;
    B[3][0] =  q3;    B[3][1] = -q2;    B[3][2] =  q1;    B[3][3] =  q0;
    
    /* transformation: From angular momentum in the solid frame to angular momentum in the body frame */
    obx = A[0][0] * (tracer+ipart)->ox + A[0][1] * (tracer+ipart)->oy + A[0][2] * (tracer+ipart)->oz;
    oby = A[1][0] * (tracer+ipart)->ox + A[1][1] * (tracer+ipart)->oy + A[1][2] * (tracer+ipart)->oz;
    obz = A[2][0] * (tracer+ipart)->ox + A[2][1] * (tracer+ipart)->oy + A[2][2] * (tracer+ipart)->oz;



   /* Now compute dq/dt = 0.5* B*(0,obx,oby,obz) , this is just a guess of dq/dt(t+0.5*time_dt) */
    dq0dt = 0.5 *  (B[0][1] * obx + B[0][2] * oby + B[0][3] * obz );
    dq1dt = 0.5 *  (B[1][1] * obx + B[1][2] * oby + B[1][3] * obz );
    dq2dt = 0.5 *  (B[2][1] * obx + B[2][2] * oby + B[2][3] * obz );
    dq3dt = 0.5 *  (B[3][1] * obx + B[3][2] * oby + B[3][3] * obz );

    q0 += 0.5*property.time_dt*dq0dt;
    q1 += 0.5*property.time_dt*dq1dt;
    q2 += 0.5*property.time_dt*dq2dt;
    q3 += 0.5*property.time_dt*dq3dt;


    /* quaternion normalization */
    
    qsq =q0*q0 + q1*q1 + q2*q2 + q3*q3;
    
    invq = 1.0/ sqrt(qsq);

    q0 *= invq;
    q1 *= invq;
    q2 *= invq;
    q3 *= invq;

   /* assigning the values back to variables */  
    (tracer+ipart)->q[0] = q0;
    (tracer+ipart)->q[1] = q1;
    (tracer+ipart)->q[2] = q2;
    (tracer+ipart)->q[3] = q3;

/* extrapolation to predict the value of ox at time t+dt */

    (tracer+ipart)->ox = 0.5 * (o12x + 3*(tracer+ipart)->ox);
    (tracer+ipart)->oy = 0.5 * (o12y + 3*(tracer+ipart)->oy);
    (tracer+ipart)->oz = 0.5 * (o12z + 3*(tracer+ipart)->oz);

