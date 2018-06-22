

/*! 
  \file md.c 
  \brief This is the integrator for finite size objects
  \defgroup md Molecular Dynamics
*/

#define LEAPFROG_ROTATION
//#define AB2
//#define FLORIAN_TORQUE

#include "md.h"
static char rev[] = "$Id$"; 

int border_nodes=0;
int *border_idx=NULL;
int border_idx_size=0;

void pbc_int(int *);

/*! Euler angle to quaternion conversion 
 \ingroup md */
quaternion euler2quaternion(euler e)
{
  quaternion q;
  q.q[0] = cos(0.5*e.theta)*cos(0.5*(e.phi+e.psi));
  q.q[1] = sin(0.5*e.theta)*cos(0.5*(e.phi-e.psi));
  q.q[2] = sin(0.5*e.theta)*sin(0.5*(e.phi-e.psi));
  q.q[3] = cos(0.5*e.theta)*sin(0.5*(e.phi+e.psi));
  return q;
}

void quaternion2euler()
{  
}


/*! 
  Initialization of MD particles configuration  
  \ingroup md 
*/
void md_init()
{
  int i;
  p = (particle_type *) malloc(sizeof(particle_type)*NMD);

  for (i=0;i<NMD;i++) {
    p[i].x = md_x+(double)i*SX/2.0;
    p[i].y = md_y;
    p[i].z = md_z;
    p[i].vx  = p[i].vy  = p[i].vz  = 0.0;
    p[i].vx0 = p[i].vy0 = p[i].vz0 = 0.0;

    p[i].ax  = p[i].ay  = p[i].az  = 0.0;
    p[i].ax0 = p[i].ay0 = p[i].az0 = 0.0;

    p[i].theta = M_PI/1.;
    p[i].psi   = 0.0;
    p[i].phi   = M_PI/2.;
    p[i].q[0] = cos(0.5*p[i].theta)*cos(0.5*(p[i].phi+p[i].psi));
    p[i].q[1] = sin(0.5*p[i].theta)*cos(0.5*(p[i].phi-p[i].psi));
    p[i].q[2] = sin(0.5*p[i].theta)*sin(0.5*(p[i].phi-p[i].psi));
    p[i].q[3] = cos(0.5*p[i].theta)*sin(0.5*(p[i].phi+p[i].psi));

    p[i].lx = p[i].ly = p[i].lz = 0.0;
    p[i].tx = p[i].ty = p[i].tz = 0.0;

    p[i].obx = p[i].oby = p[i].obz = 0.0;
    p[i].ox = p[i].oy = p[i].oz = 0.0;

#ifdef AB2
    p[i].qf[0] = 0.0;         
    p[i].qf[1] = 0.0;       
    p[i].qf[2] = 0.0;     
    p[i].qf[3] = 0.0;
    p[i].ofx = p[i].ofy = p[i].ofz = 0.0;   
#endif

#define ELLIPSOID
#define r_x (10.)
#define r_y (5.)
#define r_z (10.)
    
    p[i].m = 0.6*(4./3.*M_PI)*r_x*r_y*r_z;


#ifdef ELLIPSOID
    p[i].Ixx = p[i].m * (r_y*r_y+r_z*r_z) / 5.0;
    p[i].Iyy = p[i].m * (r_x*r_x+r_z*r_z) / 5.0;
    p[i].Izz = p[i].m * (r_x*r_x+r_y*r_y) / 5.0;
    fprintf(stderr,"Ixx = %g Iyy = %g Izz = %g\n",p[i].Ixx,p[i].Iyy,p[i].Izz);
#endif
  }

  MPI_Bcast(p,NMD,MPI_particle_type,0,MPI_COMM_WORLD);
}

/*! This subroutine will advance the particles in space and rotate them

  It implements the leap-frog integrator 

  From the LBM I get the velocity torque and acceleration at time 
  \f[ f_n = f(t_n)\;  \tau_n = \tau(t_n)\;  \f]

  \ingroup md
*/


void md_advance()
{
  int pp;
  double dt = 1.0;
  double A[3][3];
  double B[4][4];
  double q[4];
  double x, y, z;
  double ux, uy, uz;
  double vx, vy, vz;
  double vx0, vy0, vz0;
  double fx, fy, fz;
  double k_xx, k_yy, k_zz;
  double k_tranx, k_trany, k_tranz;
  double i_xx, i_yy, i_zz;
  double s_zy, s_xz;
  double w_zy, w_xz, w_yx;
  double dy_uz, dz_uy, dz_ux, dx_uz, dx_uz, dy_ux;
  double alfa, beta, gama;
  double t_x, t_y, t_z;
 
  double q0, q1, q2, q3;
  double dq0dt, dq1dt, dq2dt, dq3dt;
//  double lx, ly, lz;
 // double lbx,lby,lbz;
  double obx,oby,obz;
  double obx12p,oby12p,obz12p;
  double obx12,oby12,obz12;
  double do_x, do_y, do_z;
  double qsq, invq;
  double mu, m, a, aspr, rho;
  double v12x, v12y, v12z;
  double v12xp, v12yp, v12zp;

#ifdef AB2
  double ofx,ofy,ofz;                                         
  double qf0,qf1,qf2,qf3;   
#endif

#ifdef TIMING
  profile_on(7);
#endif
  
  for (pp=0;pp<NMD;pp++) {
    fprintf(stderr,"%d %g %g %g\n",pp,p[pp].x,p[pp].y,p[pp].z);
    /*! Euler for center of mass 
     \f[
     v_{n+1/2}=v_{n-1/2} + a_i dt
     \f]

     /* my part of the program */ 
     /* needs a lot of changes*/
     /* has to be optimized */
     /* here we compute the force required for translation of the element */
    


    /* assigning values of quaternions*/

    p[pp].q[0] = q0;
    p[pp].q[1] = q1;
    p[pp].q[2] = q2;
    p[pp].q[3] = q3;
    
    A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    A[0][1] = 2*(q1*q2+q0*q3);
    A[0][2] = 2*(q1*q3-q0*q2);
    
    A[1][0] = 2*(q1*q2-q0*q3);
    A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A[1][2] = 2*(q2*q3+q0*q1);
    
    A[2][0] = 2*(q1*q3+q0*q2);
    A[2][1] = 2*(q2*q3-q0*q1);
    A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
     


    aspr = b/a ;
    mp = 4/3*pi*pow(a,3)*aspr*rho_p ;

    K_xx = (16(pow(aspr,2) - 1))/((2*pow(aspr,2)-3) ln(aspr + sqrt(pow(aspr,2)-1)/(sqrt(pow(aspr,2) - 1))+ aspr);

    K_yy = K_xx;

    K_zz = (8(pow(aspr,2) - 1))/((2*pow(aspr,2)-1) ln(aspr + sqrt(pow(aspr,2)-1)/(sqrt(pow(aspr,2) - 1))- aspr);



/* translation dyadic or resistance tensor */
    K_tranx = (A[0][0] * K_xx)* A[0][0]+ (A[0][1] * K_xx)* A[0][1]+(A[0][2] * K_xx)* A[0][2];

    K_trany = (A[1][0] * K_yy )* A[1][0]+ ( A[1][1] * K_yy)* A[1][1]+(A[1][2] * K_yy)* A[1][2];

    K_tranz = (A[2][0] * K_zz)* A[2][0]+ ( A[2][1] * K_zz)* A[2][1]+(A[2][2] * K_zz)* A[2][2];

/* hydrodynamic drag force (brenner 1964) */

    fx = mu K_tranx * (ux - vx);
    fy = mu K_trany * (uy - vy);
    fz = mu K_tranz * (uz - vz);   

    /* Here compute v particle at time n-1/2 */
    v12x =  0.5 * (p[pp].vx + p[pp].vx0);
    v12y =  0.5 * (p[pp].vy + p[pp].vy0);
    v12z =  0.5 * (p[pp].vz + p[pp].vz0);


    /* here we assign the present time value to previous variable */
    p[pp].vx0 = p[pp].vx;
    p[pp].vy0 = p[pp].vy;
    p[pp].vz0 = p[pp].vz;

    /* Here compute v particle at time n+1/2 */

    vx = v12x + (p[pp].fx * dt) / p[pp].m;
    vy = v12y + (p[pp].fy * dt) / p[pp].m;
    vz = v12z + (p[pp].fz * dt) / p[pp].m;

    /* Here compute v particle at time n */

    /* extrapolated velocity value to 

    /*!
      Euler for center of mass 
      \f[
      x_{n+1}=x_{n} + v_{n+1/2} dt
      \f]
    */
    /*!
      Leap-frog for position 
     */
    p[pp].x += (vx * dt);
    p[pp].y += (vy * dt);
    p[pp].z += (vz * dt);
 
    /* to find the value of v at time n +1 for next iteration */ 
    p[pp].vx = 0.5 * (v12x + 3*vx);
    p[pp].vy = 0.5 * (v12y + 3*vy);
    p[pp].vz = 0.5 * (v12z + 3*vz);   

    #ifdef LEAPFROG_ROTATION
    /*Hydrodynamic torque (jeffery 1922 *for linear shear flow ) */

    I_xx = ((1+pow(aspr,2)*pow(a,2))/5*mp ;
    I_yy = I_xx ;
    I_zz = (2*pow(a,2))/5*mp ;


    S_zy = 1/2*(dy_uz + dz_uy) ;
    S_xz = 1/2*(dz_ux + dx_uz) ;


    W_zy = 1/2*(dy_uz - dz_uy) ;
    W_xz = 1/2*(dz_ux - dx_uz) ;
    W_yx = 1/2*(dx_uy - dy_ux) ;



    alfa = (pow(aspr,2)/(pow(aspr,2) - 1)) + aspr/(2*pow(pow(aspr,2)-1),3/2))* ln ((aspr - sqrt(pow(aspr,2) - 1))/(aspr + sqrt(pow(aspr,2)-1)));
    beta = alfa_0 ;
    gama = - (2/(pow(aspr,2) - 1)) - aspr/(pow(pow(aspr,2)-1),3/2))* ln ((aspr - sqrt(pow(aspr,2) - 1))/(aspr + sqrt(pow(aspr,2)-1))) ;


    T_x = (16*pi*mu*pow(a,3)*aspr)/(3*(beta + pow(aspr,2)*gama))*((1-pow(aspr,2))*S_zy + (1+ pow(aspr,2))*(W_yz - W_x)) ;

    T_y = (16*pi*mu*pow(a,3)*aspr)/(3*(pow(aspr,2)*gama + alfa))*((pow(aspr,2) - 1)*S_xz + (pow(aspr,2)+1)*(W_xz - W_y)) ;

    T_z = (32*pi*mu*pow(a,3)*aspr)/(3*(alfa + gama))*(W_yx-W_z) ;
    

    /* compute angular velocity */
 //   obx = lbx/p[pp].Ixx;
 //  oby = lby/p[pp].Iyy;
 //   obz = lbz/p[pp].Izz;

    do_x = (T_x + oy*oz*(I_yy - I_zz))/I_xx;
    do_y = (T_y + oy*oz*(I_zz - I_xx))/I_yy;
    
    /* I_xx = I_yy hence the equation is simplified */
    do_z = (T_z)/I_zz;


    /* taking average to determine the middle value */
    o12x = 0.5(p[pp].ox0+p[pp].ox);
    o12y = 0.5(p[pp].oy0+p[pp].oy);
    o12z = 0.5(p[pp].oz0+p[pp].oz);


   /*assigning the old values*/
    p[pp].ox0=p[pp].ox;
    p[pp].oy0=p[pp].oy;
    p[pp].oz0=p[pp].oz;

    /* estimating the value of omega at time t */
    p[pp].ox = o12x + (do_x * (dt*0.5));
    p[pp].oy = o12y + (do_y * (dt*0.5));
    p[pp].oz = o12z + (do_z * (dt*0.5));   

    /* assigning the value to local variable for multiple use */
    p[pp].q[0] = q0;
    p[pp].q[1] = q1;
    p[pp].q[2] = q2;
    p[pp].q[3] = q3;

    /* to transform the angular velocity to body frame */

    A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    A[0][1] = 2*(q1*q2+q0*q3);
    A[0][2] = 2*(q1*q3-q0*q2);
    
    A[1][0] = 2*(q1*q2-q0*q3);
    A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A[1][2] = 2*(q2*q3+q0*q1);
    
    A[2][0] = 2*(q1*q3+q0*q2);
    A[2][1] = 2*(q2*q3-q0*q1);
    A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    /* transformation: From angular momentum in the solid frame to angular momentum in the body frame */
    obx = A[0][0] *p[pp].ox + A[0][1] *p[pp].oy + A[0][2] *p[pp].oz;
    oby = A[1][0] *p[pp].ox + A[1][1] *p[pp].oy + A[1][2] *p[pp].oz;
    obz = A[2][0] *p[pp].ox + A[2][1] *p[pp].oy + A[2][2] *p[pp].oz;

    B[0][0] =  q0;
    B[0][1] = -q1;
    B[0][2] = -q2;
    B[0][3] = -q3;
    B[1][0] =  q1;
    B[1][1] =  q0;
    B[1][2] = -q3;
    B[1][3] =  q2;
    B[2][0] =  q2;
    B[2][1] =  q3;
    B[2][2] =  q0;
    B[2][3] = -q1;
    B[3][0] =  q3;
    B[3][1] = -q2;
    B[3][2] =  q1;
    B[3][3] =  q0;

   /* Now compute dq/dt = 0.5* B*(0,obx,oby,obz) , this is just a guess of dq/dt(t) */
    dq0dt = 0.5 *  (B[0][1] * obx + B[0][2] * oby + B[0][3] * obz );
    dq1dt = 0.5 *  (B[1][1] * obx + B[1][2] * oby + B[1][3] * obz );
    dq2dt = 0.5 *  (B[2][1] * obx + B[2][2] * oby + B[2][3] * obz );
    dq3dt = 0.5 *  (B[3][1] * obx + B[3][2] * oby + B[3][3] * obz );

    
    /* q(t+0.5*dt) = q(t) + 0.5*dt*dqdt(t) */
    p[pp].q[0] += 0.5*dt*dq0dt;
    p[pp].q[1] += 0.5*dt*dq1dt;
    p[pp].q[2] += 0.5*dt*dq2dt;
    p[pp].q[3] += 0.5*dt*dq3dt;

    
    /* quaternion normalization */
    
    qsq = 
      p[pp].q[0]*p[pp].q[0]  +
      p[pp].q[1]*p[pp].q[1]  +  
      p[pp].q[2]*p[pp].q[2]  +  
      p[pp].q[3]*p[pp].q[3];
    
    invq = 1.0/ sqrt(qsq);

    p[pp].q[0] *= invq;
    p[pp].q[1] *= invq;
    p[pp].q[2] *= invq;
    p[pp].q[3] *= invq;

    /* estimating the value of omega at time t+dt/2  */
    p[pp].ox = o12x + do_x *dt;
    p[pp].oy = o12y + do_y *dt;
    p[pp].oz = o12z + do_z *dt;
    
   /* to transform the angular velocity to body frame */
    
    q0 = p[pp].q[0];
    q1 = p[pp].q[1];
    q2 = p[pp].q[2];
    q3 = p[pp].q[3];

    A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    A[0][1] = 2*(q1*q2+q0*q3);
    A[0][2] = 2*(q1*q3-q0*q2);
    
    A[1][0] = 2*(q1*q2-q0*q3);
    A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A[1][2] = 2*(q2*q3+q0*q1);
    
    A[2][0] = 2*(q1*q3+q0*q2);
    A[2][1] = 2*(q2*q3-q0*q1);
    A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
    
    /* transformation: From angular momentum in the solid frame to angular momentum in the body frame */
    obx = A[0][0] * p[pp].ox + A[0][1] * p[pp].oy + A[0][2] * p[pp].oz;
    oby = A[1][0] * p[pp].ox + A[1][1] * p[pp].oy + A[1][2] * p[pp].oz;
    obz = A[2][0] * p[pp].ox + A[2][1] * p[pp].oy + A[2][2] * p[pp].oz;

    B[0][0] =  q0;
    B[0][1] = -q1;
    B[0][2] = -q2;
    B[0][3] = -q3;
    B[1][0] =  q1;
    B[1][1] =  q0;
    B[1][2] = -q3;
    B[1][3] =  q2;
    B[2][0] =  q2;
    B[2][1] =  q3;
    B[2][2] =  q0;
    B[2][3] = -q1;
    B[3][0] =  q3;
    B[3][1] = -q2;
    B[3][2] =  q1;
    B[3][3] =  q0;

   /* Now compute dq/dt = 0.5* B*(0,obx,oby,obz) , this is just a guess of dq/dt(t+0.5*dt) */
    dq0dt = 0.5 *  (B[0][1] * obx + B[0][2] * oby + B[0][3] * obz );
    dq1dt = 0.5 *  (B[1][1] * obx + B[1][2] * oby + B[1][3] * obz );
    dq2dt = 0.5 *  (B[2][1] * obx + B[2][2] * oby + B[2][3] * obz );
    dq3dt = 0.5 *  (B[3][1] * obx + B[3][2] * oby + B[3][3] * obz );

    p[pp].q[0] += 0.5*dt*dq0dt;
    p[pp].q[1] += 0.5*dt*dq1dt;
    p[pp].q[2] += 0.5*dt*dq2dt;
    p[pp].q[3] += 0.5*dt*dq3dt;

    /* quaternion normalization */
    
    qsq = 
      p[pp].q[0]*p[pp].q[0]  +
      p[pp].q[1]*p[pp].q[1]  +  
      p[pp].q[2]*p[pp].q[2]  +  
      p[pp].q[3]*p[pp].q[3];
    
    invq = 1.0/ sqrt(qsq);

    p[pp].q[0] *= invq;
    p[pp].q[1] *= invq;
    p[pp].q[2] *= invq;
    p[pp].q[3] *= invq;

    p[pp].ox = 0.5 * (o12x + 3*p[pp].ox);
    p[pp].oy = 0.5 * (o12y + 3*p[pp].oy);
    p[pp].oz = 0.5 * (o12z + 3*p[pp].oz);
#endif
  }
/*end o ffor loop */

void quaternion2euler()
{
phi = atan2(2*(q0*q1+q2*q3),(1-2*(pow(q1,2)+pow(q2,2))));
theta = asin(2*(q0*q2-q3*q1));
psi = atan2(2*(q0*q3+q1*q2),(1-2*(pow(q2,2)+pow(q3,2))));
}

#ifdef TIMING
  profile_off(7);
#endif
}
