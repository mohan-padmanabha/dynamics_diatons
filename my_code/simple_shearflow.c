#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int main {

      /* translation */
    my_double q0, q1, q2, q3;
    my_double matR[3][3], matK[3][3], matKT[3][3];
    my_double aspr, mp, a, rhop, mu, pi;
    my_double k_tranx, k_trany, k_tranz;
    my_double fx, fy, fz;
    my_double v12x, v12y, v12z;

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
    my_double qt0, qt1, qt2, qt3;
    my_double ox,oy,oz;
    my_double ox_old,oy_old,oz_old;


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
   uz = 0;
   du_dx = 2;
   

/* conversion from axis angles to quaternions */
 q1 =  px * sin(theta/2);
 q2 =  py * sin(theta/2);
 q3 =  pz * sin(theta/2);
 q0 = cos(theta/2);


     /* rotation matrix of quaternions */
    matR[0][0] = 1-2*(q2*q2+q3*q3);  matR[1][0] = 2*(q1*q2-q0*q3);   matR[2][0] = 2*(q1*q3+q0*q2); 
    matR[0][1] = 2*(q1*q2+q0*q3);   matR[1][1] = 1-2*(q1*q1+q3*q3);  matR[2][1] = 2*(q2*q3-q0*q1);
    matR[0][2] = 2*(q1*q3-q0*q2);   matR[1][2] = 2*(q2*q3+q0*q1);   matR[2][2] = 1-2*(q1*q1+q2*q2);  
 
     

 // aspect ratio is major axis / minor axis
 // already present in the code 
//     aspect_ratio = b/a ;
// mp is the mass of the particle (zhang et al 2001)
    mp = 4/3*pi*pow(a,3)*aspr*rhop ;

    matK[0][0] = (16*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr*aspr - 3)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));

    matK[1][1] = matK[0][0];

    matK[2][2] = (8*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr*aspr - 1)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));



/* translation dyadic or resistance tensor */

          /* For transformation of co-ordinates from Lab frome to body frame */
          /* general form  */
               for (m=1; m<3; m++){
                     p = 0;
                 for (l=0; l<3; l++){
                     r = 0;
                   for (k=0; k<3; k++){
                     r += matR[l][k]*matK[k][l];
                    }
                     matKT[m][l] += r * matR[m][l];
                      }
                         }
                       

/* hydrodynamic drag force (brenner 1964 && zhang et al 2001) */
 

    fx = property.nu * matKT[0][0] * ( ux -  vx);
    fy = property.nu * matKT[1][1] * ( uy -  vy);
    fz = property.nu * matKT[2][2] * ( uz -  vz); 



    /* Here compute v particle at time n-1/2 */
    v12x =  0.5 * ( vx +  vx_old);
    v12y =  0.5 * ( vy +  vy_old);
    v12z =  0.5 * ( vz +  vz_old);

    /* here we assign the present time value to previous variable */
     vx_old =  vx;
     vy_old =  vy;
     vz_old =  vz;

    /* Here compute v particle at time n+1/2 */

     vx = v12x + (fx * property.time_dt) / mp;
     vy = v12y + (fy * property.time_dt) / mp;
     vz = v12z + (fz * property.time_dt) / mp;

    /* Here compute v particle at time n */

      /* computation of the position */
     x += ( vx * property.time_dt);
     y += ( vy * property.time_dt);
     z += ( vz * property.time_dt);
 
    /* to find the value of v  by extrapolation at time n +1 for next iteration */ 
     vx = 0.5 * (v12x + 3* vx);
     vy = 0.5 * (v12y + 3* vy);
     vz = 0.5 * (v12z + 3* vz);   

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
    matD[0][0]= dx_ux ; matD[0][1]= dy_ux; matD[0][2]= dz_ux;
    matD[1][0]= dx_uy ; matD[1][1]= dy_uy; matD[1][2]= dz_uy;
    matD[2][0]= dx_uz ; matD[2][1]= dy_uz; matD[2][2]= dz_uz;	 

    /* For transformation of co-ordinates from Lab frome to body frame for tensor */

      for (m=1; m<3; m++){
           p = 0;
           for (l=0; l<3; l++){
               r = 0;
               for (k=0; k<3; k++){
                   r += matR[l][k]*matD[k][l];
               }
               matDT[m][l] += r * matR[m][l];
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

    T_x = (16*pi*property.nu*pow(a,3)*aspr)/(3*(beta0 + (aspr*aspr)*gamma0)) * ((1-(aspr*aspr)) * matSB[2][1] + (1+ (aspr*aspr)) * ( matWB[1][2] -  ox)) ;

    T_y = (16*pi*property.nu*pow(a,3)*aspr)/(3*((aspr*aspr)*gamma0 + alpha0))*(((aspr*aspr) - 1) * matSB[0][2] + ((aspr*aspr)+1) * ( matWB[0][2] -  oy)) ;

    T_z = (32*pi*property.nu*pow(a,3)*aspr)/(3*(alpha0 + gamma0))*( matWB[1][0] -  oz) ;
    

    /* compute angular velocity */
    //   obx = lbx/ Ixx;
    //  oby = lby/ Iyy;
    //   obz = lbz/ Izz;

    do_x = (T_x +  oy* oz*(I_yy - I_zz))/I_xx;
    do_y = (T_y +  oz*  ox *(I_zz - I_xx))/I_yy;

    /* I_xx = I_yy hence the equation is simplified */
    do_z = (T_z)/I_zz;


    /* taking average to determine the middle value */
    o12x = 0.5*( ox_old+ ox);
    o12y = 0.5*( oy_old+ oy);
    o12z = 0.5*( oz_old+ oz);


   /*assigning the old values*/
     ox_old= ox;
     oy_old= oy;
     oz_old= oz;

    /* estimating the value of omega at time t */
     ox = o12x + (do_x * (property.time_dt*0.5));
     oy = o12y + (do_y * (property.time_dt*0.5));
     oz = o12z + (do_z * (property.time_dt*0.5));   

    /* transformation: From angular momentum in the lab frame to body frame */
    obx = matR[0][0] * ox + matR[0][1] * oy + matR[0][2] * oz;
    oby = matR[1][0] * ox + matR[1][1] * oy + matR[1][2] * oz;
    obz = matR[2][0] * ox + matR[2][1] * oy + matR[2][2] * oz;

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
     ox = o12x + do_x *property.time_dt;
     oy = o12y + do_y *property.time_dt;
     oz = o12z + do_z *property.time_dt;


     /* rotation matrix of quaternions */
    matR[0][0] = 1-2*(q2*q2+q3*q3);  matR[1][0] = 2*(q1*q2-q0*q3);   matR[2][0] = 2*(q1*q3+q0*q2); 
    matR[0][1] = 2*(q1*q2+q0*q3);   matR[1][1] = 1-2*(q1*q1+q3*q3);  matR[2][1] = 2*(q2*q3-q0*q1);
    matR[0][2] = 2*(q1*q3-q0*q2);   matR[1][2] = 2*(q2*q3+q0*q1);   matR[2][2] = 1-2*(q1*q1+q2*q2);  


    /* multiplication matrix of quaternions to advance in time */
    B[0][0] =  q0;    B[0][1] = -q1;    B[0][2] = -q2;    B[0][3] = -q3;    
    B[1][0] =  q1;    B[1][1] =  q0;    B[1][2] = -q3;    B[1][3] =  q2;
    B[2][0] =  q2;    B[2][1] =  q3;    B[2][2] =  q0;    B[2][3] = -q1;
    B[3][0] =  q3;    B[3][1] = -q2;    B[3][2] =  q1;    B[3][3] =  q0;
    
    /* transformation: From angular momentum in the solid frame to angular momentum in the body frame */
    obx = matR[0][0] *  ox + matR[0][1] *  oy + matR[0][2] *  oz;
    oby = matR[1][0] *  ox + matR[1][1] *  oy + matR[1][2] *  oz;
    obz = matR[2][0] *  ox + matR[2][1] *  oy + matR[2][2] *  oz;



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

     ox = 0.5 * (o12x + 3* ox);
     oy = 0.5 * (o12y + 3* oy);
     oz = 0.5 * (o12z + 3* oz);


/* conversion from quaternions to axis angles */
theta = 2 * acos(q0) ;
 px = q1 / sqrt(1-q0*q0) ;
 py = q2 / sqrt(1-q0*q0) ;
 pz = q3 / sqrt(1-q0*q0) ;


   /* assigning the values back to variables */  


}


