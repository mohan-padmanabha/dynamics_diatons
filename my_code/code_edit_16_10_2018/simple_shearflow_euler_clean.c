#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



int main() {

      /* translation */
    double q0, q1, q2, q3;
    double matR[3][3], matK[3][3], matKT[3][3];
    double aspr, mp, a, rhop, mu, pi, nu, time_dt;
    double k_tranx, k_trany, k_tranz;
    double fx, fy, fz;
    double vx,vy,vz;
    double v2x, v2y, v2z;
    double x, y, z;
    double ux, uy, uz;

/* orinetation */
    double I_xx, I_yy, I_zz;
    double matD[3][3], matSB[3][3], matDT[3][3], matWB[3][3], B[4][4] ;
    double alpha0, beta0, gamma0, theta;
    double T_x, T_y, T_z;
    double do_x, do_y, do_z;
    double obx, oby, obz;
    double dq0dt, dq1dt, dq2dt, dq3dt;
    double qsq, invq, p, r;
    int l, m, k, i;
    double qt0, qt1, qt2, qt3;
    double ox,oy,oz;
    double o2x, o2y, o2z;
    double dlo_x, dlo_y, dlo_z;


   FILE * qdata = fopen("quat_data_red.txt","w");
   FILE * vdata = fopen("vec_data_red.txt","w");


//  LAGRANGE_ORIENTATION_TRANSLATION

   pi = 3.1415;
   theta = 0. ;
   aspr = 4;
   nu = 0.0000015 ;
   double bx = 0.001;
   rhop = 1.1;
   time_dt = 0.1;
   a = aspr * bx ;
   x = 1;
   y = 1;
   z = 0;

   ux = 0;
   uy = 0;
   uz = 0;
   // double du_dx = 1;
   double x_max = 5;
   double y_max = 5;
   int itr = 10000;
   double px = 0;
   double py = 0;
   double pz = 0;

   vx = 0;
   vy = 0;
   vz = 0;

   ox = 1.0e-5;
   oy = 1.0e-5;
   oz = 0.0;


    double dx_ux = 0 ; double dy_ux = 1e-3 ; double dz_ux = 0;
    double dx_uy = 0 ; double dy_uy = 0; double dz_uy = 0;
    double dx_uz = 0 ; double dy_uz = 0 ; double dz_uz = 0;

    
/* conversion from axis angles to quaternions */

 q1 =  px * sin(theta/2);
 q2 =  py * sin(theta/2);
 q3 =  pz * sin(theta/2);
 q0 = cos(theta/2);


 // aspect ratio is major axis / minor axis
//     aspect_ratio = b/a ;
// mp is the mass of the particle (zhang et al 2001)
    mp = 4/3*pi*pow(a,3)*aspr*rhop ;

    matK[0][0] = (16*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr - 3)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));

    matK[1][1] = matK[0][0];

    matK[2][2] = (8*pi*a*pow((aspr*aspr - 1),3/2))/((2*aspr - 1)*log(aspr+pow((aspr*aspr-1),1/2)) + aspr*pow((aspr*aspr-1),1/2));

     printf("matK=%e %e %e\n", matK[0][0], matK[1][1], matK[2][2]);
	
    
     
     /* moment of inertia (diagnal in body frame)*/
    I_xx = ((1+(aspr*aspr))*(a*a)*mp)/(5) ;
    I_yy = I_xx ;
    I_zz = (2*a*a*mp)/(5) ;



    /* computation of constants alpha0 beta0 and gamma0 to compute torque */

    alpha0 = (2*(aspr*aspr)*pow(((aspr*aspr) - 1),1/2) + aspr*log((aspr-pow(((aspr*aspr)-1),1/2))/(aspr+pow(((aspr*aspr)-1),1/2)))) / (2* pow(((aspr*aspr)-1),3/2));
    beta0 = alpha0 ;
    gamma0 = (2*pow(((aspr*aspr) - 1),1/2) + aspr*log((aspr-pow(((aspr*aspr)-1),1/2))/(aspr+pow(((aspr*aspr)-1),1/2)))) / (pow(((aspr*aspr)-1),3/2));

for(i=0; i < itr ; i++){



    uy = dy_ux * y;
    ux = uy;
     /* rotation matrix of quaternions */
    matR[0][0] = (q0*q0)+(q1*q1)-(q2*q2)-(q3*q3);  matR[1][0] = 2*(q1*q2+q0*q3);   matR[2][0] = 2*(q1*q3-q0*q2); 
    matR[0][1] = 2*(q1*q2-q0*q3);   matR[1][1] = (q0*q0)-(q1*q1)+(q2*q2)-(q3*q3);  matR[2][1] = 2*(q2*q3+q0*q1);
    matR[0][2] = 2*(q1*q3+q0*q2);   matR[1][2] = 2*(q2*q3-q0*q1);   matR[2][2] = (q0*q0)-(q1*q1)-(q2*q2)+(q3*q3);  
   
     printf("matR %e %e %e %e \n", matR[0][0], matR[0][1], matR[0][2], matR[1][0]);
     


/* translation dyadic or resistance tensor */

          /* For transformation of co-ordinates from Lab frome to body frame */
          
               for (m=0; m<3; m++){
                     p = 0;
                 for (l=0; l<3; l++){
                     r = 0;
                   for (k=0; k<3; k++){
                     r += matR[m][k]*matK[k][l];                    
                    }
                     p += r * matR[m][l];
                     matKT[m][l] = p;
 
                      }  
                     
                  }
                       

/* hydrodynamic drag force (brenner 1964 && zhang et al 2001) */
 
    fx = (matKT[0][0]*(ux - vx) +  matKT[0][1]*(uy - vy) +  matKT[0][2]*(uz - vz))*nu; 
    fy = (matKT[1][0]*(ux - vx) +  matKT[1][1]*(uy - vy) +  matKT[1][2]*(uz - vz))*nu;
    fz = (matKT[2][0]*(ux - vx) +  matKT[2][1]*(uy - vy) +  matKT[2][2]*(uz - vz))*nu;
     printf("force =%e %e %e\n", fx, fy, fz);

 

    /* Here compute v particle at time n+1/2 */

     v2x = vx + (fx * (0.5*time_dt)) / mp;
     v2y = vy + (fy * (0.5*time_dt)) / mp;
     v2z = vz + (fz * (0.5*time_dt)) / mp;
     
    /* Here compute v particle at time n */

      /* computation of the position */
     x += ( v2x * time_dt);
     y += ( v2y * time_dt);
     z += ( v2z * time_dt);
 

/* computing the velocity at time n+1 using euler */
     vx += (fx * (time_dt)) / mp;
     vy += (fy * (time_dt)) / mp;
     vz += (fz * (time_dt)) / mp;
     
     
    /* to find the value of v  by extrapolation at time n +1 for next iteration */ 
//     vx = 0.5 * (v12x + 3* vx);
//     vy = 0.5 * (v12y + 3* vy);
//     vz = 0.5 * (v12z + 3* vz);   




     /* -------------------------------------------------------------------------*/
     /* orientation computation */

    /* multiplication matrix of quaternions to advance in time */
    B[0][0] =  q0;    B[0][1] = -q1;    B[0][2] = -q2;    B[0][3] = -q3;    
    B[1][0] =  q1;    B[1][1] =  q0;    B[1][2] = -q3;    B[1][3] =  q2;
    B[2][0] =  q2;    B[2][1] =  q3;    B[2][2] =  q0;    B[2][3] = -q1;
    B[3][0] =  q3;    B[3][1] = -q2;    B[3][2] =  q1;    B[3][3] =  q0;


    /* velocity gradient matrix */
    matD[0][0]= dx_ux ; matD[0][1]= dy_ux; matD[0][2]= dz_ux;
    matD[1][0]= dx_uy ; matD[1][1]= dy_uy; matD[1][2]= dz_uy;
    matD[2][0]= dx_uz ; matD[2][1]= dy_uz; matD[2][2]= dz_uz;	 
    
    
    
    /* transformation: From angular momentum in the lab frame to body frame */
    obx = matR[0][0] * ox + matR[0][1] * oy + matR[0][2] * oz;
    oby = matR[1][0] * ox + matR[1][1] * oy + matR[1][2] * oz;
    obz = matR[2][0] * ox + matR[2][1] * oy + matR[2][2] * oz;

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
                matWB[k][l] = 0.5*(matDT[k][l]-matDT[l][k]);
               }



    /* computation of torque using the values computed */

    T_x = (16*pi*nu*pow(a,3)*aspr)/(3*(beta0 + (aspr*aspr)*gamma0)) * ((1-(aspr*aspr)) * matSB[2][1] + (1+ (aspr*aspr)) * ( matWB[2][1] -  obx)) ;

    T_y = (16*pi*nu*pow(a,3)*aspr)/(3*((aspr*aspr)*gamma0 + alpha0))*(((aspr*aspr) - 1) * matSB[0][2] + ((aspr*aspr)+1) * ( matWB[0][2] -  oby)) ;

    T_z = (32*pi*nu*pow(a,3)*aspr)/(3*(alpha0 + beta0))*( matWB[0][1] -  obz) ;
    

    /* computing the deravative of angular momentum for computing the next time step */    

    do_x = (T_x +  oby* obz*(I_yy - I_zz))/I_xx;
    do_y = (T_y +  obz*  obx *(I_zz - I_xx))/I_yy;

    /* I_xx = I_yy hence the equation is simplified */
    do_z = (T_z)/I_zz;

    /* transformation: From angular momentum in the body frame to lab frame using do_x * R(transpose)  */
    dlo_x = matR[0][0] * do_x + matR[1][0] * do_y + matR[2][0] * do_z;
    dlo_y = matR[0][1] * do_x + matR[1][1] * do_y + matR[2][1] * do_z;
    dlo_z = matR[0][2] * do_x + matR[1][2] * do_y + matR[2][2] * do_z;

    



   /* Now compute dq/dt = 0.5* B*(0,obx,oby,obz) , this is just a guess of dq/dt(t) */
    dq0dt = 0.5 *  (B[0][1] * obx + B[0][2] * oby + B[0][3] * obz );
    dq1dt = 0.5 *  (B[1][1] * obx + B[1][2] * oby + B[1][3] * obz );
    dq2dt = 0.5 *  (B[2][1] * obx + B[2][2] * oby + B[2][3] * obz );
    dq3dt = 0.5 *  (B[3][1] * obx + B[3][2] * oby + B[3][3] * obz );

    
    /* advancement in time  q(t+0.5*time_dt) = q(t) + 0.5*time_dt*dqdt(t) */
    q0 += 0.5*time_dt*dq0dt;
    q1 += 0.5*time_dt*dq1dt;
    q2 += 0.5*time_dt*dq2dt;
    q3 += 0.5*time_dt*dq3dt;

    
    /* quaternion normalization */
    
     qsq =q0*q0 + q1*q1 + q2*q2 + q3*q3;
    
    invq = 1.0/ sqrt(qsq);

    q0 *= invq;
    q1 *= invq;
    q2 *= invq;
    q3 *= invq;

    
    /* assigning values for computing later */
    
    o2x = ox;
    o2y = oy;
    o2z = oz;
    
    /* estimating the value of omega at time t+dt/2  */
     ox += dlo_x *(time_dt*0.5);
     oy += dlo_y *(time_dt*0.5);
     oz += dlo_z *(time_dt*0.5);


     /* rotation matrix of quaternions updated time step t + dt/2 */
    matR[0][0] = (q0*q0)+(q1*q1)-(q2*q2)-(q3*q3);  matR[1][0] = 2*(q1*q2+q0*q3);   matR[2][0] = 2*(q1*q3-q0*q2); 
    matR[0][1] = 2*(q1*q2-q0*q3);   matR[1][1] = (q0*q0)-(q1*q1)+(q2*q2)-(q3*q3);  matR[2][1] = 2*(q2*q3+q0*q1);
    matR[0][2] = 2*(q1*q3+q0*q2);   matR[1][2] = 2*(q2*q3-q0*q1);   matR[2][2] = (q0*q0)-(q1*q1)-(q2*q2)+(q3*q3); 


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

    
    /* advancement in time  q(t+time_dt) = q(t+0.5*time_dt) + 0.5*time_dt*dqdt(t+0.5*time_dt) */
    q0 += 0.5*time_dt*dq0dt;
    q1 += 0.5*time_dt*dq1dt;
    q2 += 0.5*time_dt*dq2dt;
    q3 += 0.5*time_dt*dq3dt;


    /* quaternion normalization */
    
    qsq =q0*q0 + q1*q1 + q2*q2 + q3*q3;
    
    invq = 1.0/ sqrt(qsq);

    q0 *= invq;
    q1 *= invq;
    q2 *= invq;
    q3 *= invq;


    
    
     /* estimating the value of omega at time t+dt using euler values at time t  */
     ox = o2x + dlo_x *(time_dt);
     oy = o2y + dlo_y *(time_dt);
     oz = o2z + dlo_z *(time_dt);   
    
    
    /* extrapolation to predict the value of ox at time t+dt */

//     ox = 0.5 * (o12x + 3* ox);
//     oy = 0.5 * (o12y + 3* oy);
//     oz = 0.5 * (o12z + 3* oz);


/* conversion from quaternions to axis angles */
 theta = 2 * acos(q0) ;
 px = q1 / sqrt(1-q0*q0) ;
 py = q2 / sqrt(1-q0*q0) ;
 pz = q3 / sqrt(1-q0*q0) ;




 printf("q0=%e q1=%e q2=%e q3=%e \n", q0, q1, q2, q3);
 printf("mp = %e\t KT = %e\t ux = %e\t fx = %e\t vx = %e\t x = %e\t ixx = %e\t tx = %e\t ox = %e\n", mp, matKT[0][0], ux, fx, vx, x, I_xx, T_x, ox); 
 printf("mp = %e\t KT = %e\t uy = %e\t fy = %e\t vy = %e\t y = %e\t iyy = %e\t ty = %e\t oy = %e\n", mp, matKT[1][1], uy, fy, vy, y, I_yy, T_y, oy); 
 printf("mp = %e\t KT = %e\t uz = %e\t fz = %e\t vz = %e\t z = %e\t izz = %e\t tz = %e\t oz = %e\n", mp, matKT[2][2], uz, fz, vz, z, I_zz, T_z, oz); 
 fprintf(qdata, "%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n", q0, q1, q2, q3, px, py, pz, x, y, z);
 fprintf(vdata, "%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t  %e\t %e\t %e\n", ox, oy, oz, vx, vy, vz, x, y, z, fx, fy, fz, T_x, T_y, T_z);

 if(x > 5) {x = 0;}
 if(y > 5) {y = 0;}
 if(z > 5) {z = 0;}

}

   /* assigning the values back to variables */  
return 0;

}



