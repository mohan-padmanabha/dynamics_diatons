  #ifdef LAGRANGE_ORIENTATION_TRANSLATION
      /* translation */
    my_double q0, q1, q2, q3;
    my_double matR[3][3], matK[3][3], matKT[3][3];
    my_double aspr, mp, a, rhop, mu, pi;
    my_double k_tranx, k_trany, k_tranz;
    my_double fx, fy, fz;
    my_double v12x, v12y, v12z;


 #endif   


/* orinetation */
    my_double I_xx, I_yy, I_zz;
    my_double matD[3][3], matSB[3][3], matDT[3][3], matWB[3][3], B[4][4] ;
    my_double alpha0, beta0, gamma0, theta;
    my_double T_x, T_y, T_z;
    my_double do_x, do_y, do_z;
    my_double o12x, o12y, o12z;
    my_double obx, oby, obz;
    my_double dq0dt, dq1dt, dq2dt, dq3dt;
    my_double qsq, invq, p, r;
    int l, m, k;
 #endif

  #ifdef LAGRANGE_ORIENTATION_TRANSLATION
  my_double qt0, qt1, qt2, qt3;
  my_double ox,oy,oz;
  my_double ox_old,oy_old,oz_old;
  #endif

/* initial conditions */
#ifdef LAGRANGE_ORIENTATION
/*
(tracer+i)->px = 0.0;
(tracer+i)->py = 1.0;
(tracer+i)->pz = 0.0;
*/
/* uniform random oriented vector */
/*
theta = acos(2.0*myrand()-1.0);
phi = two_pi*myrand();
(tracer+i)->px = sin(theta)*sin(phi);
(tracer+i)->py = sin(theta)*cos(phi);
(tracer+i)->pz = cos(theta);

 vec = random_vector();
*/


  #ifdef LAGRANGE_ORIENTATION_TRANSLATION

   pi = 3.1415;
   theta = 0. ;

/* conversion from axis angles to quaternions */
(tracer+ipart)->qt1 = (tracer+ipart)->px * sin(theta/2);
(tracer+ipart)->qt2 = (tracer+ipart)->py * sin(theta/2);
(tracer+ipart)->qt3 = (tracer+ipart)->pz * sin(theta/2);
(tracer+ipart)->qt0 = cos(theta/2);



    /* assigning to the local variable */
    q0 = (tracer+ipart)->qt0;
    q1 = (tracer+ipart)->qt1;
    q2 = (tracer+ipart)->qt2;
    q3 = (tracer+ipart)->qt3;


     /* rotation matrix of quaternions */
    matR[0][0] = (q0*q0)+(q1*q1)-(q2*q2)-(q3*q3);  matR[1][0] = 2*(q1*q2+q0*q3);   matR[2][0] = 2*(q1*q3-q0*q2); 
    matR[0][1] = 2*(q1*q2-q0*q3);   matR[1][1] = (q0*q0)-(q1*q1)+(q2*q2)-(q3*q3);  matR[2][1] = 2*(q2*q3+q0*q1);
    matR[0][2] = 2*(q1*q3+q0*q2);   matR[1][2] = 2*(q2*q3-q0*q1);   matR[2][2] = (q0*q0)-(q1*q1)-(q2*q2)+(q3*q3);   
 
     

 // aspect ratio is major axis / minor axis
 // already present in the code 
//    (tracer+ipart)->aspect_ratio = b/a ;
  aspr = (tracer+ipart)->aspect_ratio ; 
// mp is the mass of the particle (zhang et al 2001)
    mp = 4/3*pi*pow(a,3)*aspr*rhop ;



    matK[0][0] = (16*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr - 3)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));

    matK[1][1] = matK[0][0];

    matK[2][2] = (8*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr - 1)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));

/* translation dyadic or resistance tensor */

          /* For transformation of co-ordinates from Lab frome to body frame */
          /* general form  */
               for (m=0; m<3; m++){
                     p = 0;
                 for (l=0; l<3; l++){
                     r = 0;
                   for (k=0; k<3; k++){
                     r += matR[m][k]*matK[k][l];                    
                    }
                     p += r * matR[m][l];
                     matKT[m][l] = p;
                    //  printf("%e %e %e %e %e\n",matKT[m][l], vx, ux, x, y);
 
                      }  
                     
                  }                       

/* hydrodynamic drag force (brenner 1964 && zhang et al 2001) */
 
    fx = (matKT[0][0]*((tracer+ipart)->ux - (tracer+ipart)->vx) +  matKT[0][1]*((tracer+ipart)->uy - (tracer+ipart)->vy) +  matKT[0][2]*((tracer+ipart)->uz - (tracer+ipart)->vz))*property.nu; 
    fy = (matKT[1][0]*((tracer+ipart)->ux - (tracer+ipart)->vx) +  matKT[1][1]*((tracer+ipart)->uy - (tracer+ipart)->vy) +  matKT[1][2]*((tracer+ipart)->uz - (tracer+ipart)->vz))*property.nu;
    fz = (matKT[2][0]*((tracer+ipart)->ux - (tracer+ipart)->vx) +  matKT[2][1]*((tracer+ipart)->uy - (tracer+ipart)->vy) +  matKT[2][2]*((tracer+ipart)->uz - (tracer+ipart)->vz))*property.nu;

 //   fx = property.nu * matKT[0][0] * ((tracer+ipart)->ux - (tracer+ipart)->vx);
 //   fy = property.nu * matKT[1][1] * ((tracer+ipart)->uy - (tracer+ipart)->vy);
 //   fz = property.nu * matKT[2][2] * ((tracer+ipart)->uz - (tracer+ipart)->vz); 



    /* Here compute v particle at time n-1/2 */
    v12x =  0.5 * ((tracer+ipart)->vx + (tracer+ipart)->vx_old);
    v12y =  0.5 * ((tracer+ipart)->vy + (tracer+ipart)->vy_old);
    v12z =  0.5 * ((tracer+ipart)->vz + (tracer+ipart)->vz_old);

    /* here we assign the present time value to previous variable */
    (tracer+ipart)->vx_old = (tracer+ipart)->vx;
    (tracer+ipart)->vy_old = (tracer+ipart)->vy;
    (tracer+ipart)->vz_old = (tracer+ipart)->vz;

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

     /* ----------------------------------------------------------------------------------------------------*/
     /* orientation computation */

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

      for (m=0; m<3; m++){
           p = 0;
           for (l=0; l<3; l++){
               r = 0;
               for (k=0; k<3; k++){
                   r += matR[l][k]*matD[k][l];
               }
               p += r * matR[m][l];
               matDT[m][l] = p;
            }
             
       }

     /* Compute Sij for body frame*/
    /*make simmetric*/

	for (k=0; k<3; k++)
            for (l=0; l<3; l++){
                matSB[k][l] = 0.5*(matDT[k][l]+matDT[l][k]);
               }

    /* Compute Wij for body frame*/
    /*make simmetric*/


	for (k=0; k<3; k++)
            for (l=0; l<3; l++){
                matWB[k][l] = 0.5*(matDT[k][l]+matDT[l][k]);
               }

    /* computation of constants alpha0 beta0 and gamma0 to compute torque */

    alpha0 = ((aspr*aspr)/((aspr*aspr) - 1)) + aspr/(2*pow(((aspr*aspr)-1),3/2)* log((aspr - sqrt((aspr*aspr) - 1))/(aspr + sqrt((aspr*aspr)-1))));
    beta0 = alpha0 ;
    gamma0 = - (2/((aspr*aspr) - 1)) - aspr/(pow(((aspr*aspr)-1),3/2)* log((aspr - sqrt((aspr*aspr) - 1))/(aspr + sqrt((aspr*aspr)-1))));

    /* computation of torque using the values computed */

    T_x = (16*pi*property.nu*pow(a,3)*aspr)/(3*(beta0 + (aspr*aspr)*gamma0)) * ((1-(aspr*aspr)) * matSB[2][1] + (1+ (aspr*aspr)) * ( matWB[1][2] - (tracer+ipart)->ox)) ;

    T_y = (16*pi*property.nu*pow(a,3)*aspr)/(3*((aspr*aspr)*gamma0 + alpha0))*(((aspr*aspr) - 1) * matSB[0][2] + ((aspr*aspr)+1) * ( matWB[0][2] - (tracer+ipart)->oy)) ;

    T_z = (32*pi*property.nu*pow(a,3)*aspr)/(3*(alpha0 + gamma0))*( matWB[1][0] - (tracer+ipart)->oz) ;
    

    /* compute angular velocity */
    //   obx = lbx/(tracer+ipart)->Ixx;
    //  oby = lby/(tracer+ipart)->Iyy;
    //   obz = lbz/(tracer+ipart)->Izz;

    do_x = (T_x + (tracer+ipart)->oy*(tracer+ipart)->oz*(I_yy - I_zz))/I_xx;
    do_y = (T_y + (tracer+ipart)->oz* (tracer+ipart)->ox *(I_zz - I_xx))/I_yy;

    /* I_xx = I_yy hence the equation is simplified */
    do_z = (T_z)/I_zz;


    /* taking average to determine the middle value */
    o12x = 0.5*((tracer+ipart)->ox_old+(tracer+ipart)->ox);
    o12y = 0.5*((tracer+ipart)->oy_old+(tracer+ipart)->oy);
    o12z = 0.5*((tracer+ipart)->oz_old+(tracer+ipart)->oz);


   /*assigning the old values*/
    (tracer+ipart)->ox_old=(tracer+ipart)->ox;
    (tracer+ipart)->oy_old=(tracer+ipart)->oy;
    (tracer+ipart)->oz_old=(tracer+ipart)->oz;

    /* estimating the value of omega at time t */
    (tracer+ipart)->ox = o12x + (do_x * (property.time_dt*0.5));
    (tracer+ipart)->oy = o12y + (do_y * (property.time_dt*0.5));
    (tracer+ipart)->oz = o12z + (do_z * (property.time_dt*0.5));   

    /* transformation: From angular momentum in the lab frame to body frame */
    obx = matR[0][0] *(tracer+ipart)->ox + matR[0][1] *(tracer+ipart)->oy + matR[0][2] *(tracer+ipart)->oz;
    oby = matR[1][0] *(tracer+ipart)->ox + matR[1][1] *(tracer+ipart)->oy + matR[1][2] *(tracer+ipart)->oz;
    obz = matR[2][0] *(tracer+ipart)->ox + matR[2][1] *(tracer+ipart)->oy + matR[2][2] *(tracer+ipart)->oz;

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


     /* rotation matrix of quaternions */
    matR[0][0] = (q0*q0)+(q1*q1)-(q2*q2)-(q3*q3);  matR[1][0] = 2*(q1*q2+q0*q3);   matR[2][0] = 2*(q1*q3-q0*q2); 
    matR[0][1] = 2*(q1*q2-q0*q3);   matR[1][1] = (q0*q0)-(q1*q1)+(q2*q2)-(q3*q3);  matR[2][1] = 2*(q2*q3+q0*q1);
    matR[0][2] = 2*(q1*q3+q0*q2);   matR[1][2] = 2*(q2*q3-q0*q1);   matR[2][2] = (q0*q0)-(q1*q1)-(q2*q2)+(q3*q3);    


    /* multiplication matrix of quaternions to advance in time */
    B[0][0] =  q0;    B[0][1] = -q1;    B[0][2] = -q2;    B[0][3] = -q3;    
    B[1][0] =  q1;    B[1][1] =  q0;    B[1][2] = -q3;    B[1][3] =  q2;
    B[2][0] =  q2;    B[2][1] =  q3;    B[2][2] =  q0;    B[2][3] = -q1;
    B[3][0] =  q3;    B[3][1] = -q2;    B[3][2] =  q1;    B[3][3] =  q0;
    
    /* transformation: From angular momentum in the solid frame to angular momentum in the body frame */
    obx = matR[0][0] * (tracer+ipart)->ox + matR[0][1] * (tracer+ipart)->oy + matR[0][2] * (tracer+ipart)->oz;
    oby = matR[1][0] * (tracer+ipart)->ox + matR[1][1] * (tracer+ipart)->oy + matR[1][2] * (tracer+ipart)->oz;
    obz = matR[2][0] * (tracer+ipart)->ox + matR[2][1] * (tracer+ipart)->oy + matR[2][2] * (tracer+ipart)->oz;



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

/* extrapolation to predict the value of ox at time t+dt */

    (tracer+ipart)->ox = 0.5 * (o12x + 3*(tracer+ipart)->ox);
    (tracer+ipart)->oy = 0.5 * (o12y + 3*(tracer+ipart)->oy);
    (tracer+ipart)->oz = 0.5 * (o12z + 3*(tracer+ipart)->oz);


/* conversion from quaternions to axis angles */
theta = 2 * acos(q0) ;
(tracer+ipart)->px = q1 / sqrt(1-q0*q0) ;
(tracer+ipart)->py = q2 / sqrt(1-q0*q0) ;
(tracer+ipart)->pz = q3 / sqrt(1-q0*q0) ;


   /* assigning the values back to variables */  
    (tracer+ipart)->qt0 = q0;
    (tracer+ipart)->qt1 = q1;
    (tracer+ipart)->qt2 = q2;
    (tracer+ipart)->qt3 = q3;





