 

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

   /* assigning the value to local variable for multiple use */
    q0 =  (tracer+ipart)->q[0];
    q1 =  (tracer+ipart)->q[1];
    q2 =  (tracer+ipart)->q[2];
    q3 =  (tracer+ipart)->q[3];

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



