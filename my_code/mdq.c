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
  double q0, q1, q2, q3;
  double dq0dt, dq1dt, dq2dt, dq3dt;
  double lx, ly, lz;
  double lbx,lby,lbz;
  double obx,oby,obz;
  double obx12p,oby12p,obz12p;

  double qsq, invq;

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

     \code
     p[pp].vx += (p[pp].ax * dt) / p[pp].m;
     p[pp].vy += (p[pp].ay * dt) / p[pp].m;
     p[pp].vz += (p[pp].az * dt) / p[pp].m;
     \endcode

    */

    /* Here compute v particle at time n-1/2 */
    v12x =  0.5 * (p[pp].vx + p[pp].vx0);
    v12y =  0.5 * (p[pp].vy + p[pp].vy0);
    v12z =  0.5 * (p[pp].vz + p[pp].vz0);

    p[pp].vx0 = p[pp].vx;
    p[pp].vy0 = p[pp].vy;
    p[pp].vz0 = p[pp].vz;

    /* Here compute v particle at time n+1/2 */

    v12xp = v12x + (p[pp].ax * dt) / p[pp].m;
    v12yp = v12y + (p[pp].ay * dt) / p[pp].m;
    v12zp = v12z + (p[pp].az * dt) / p[pp].m;

    /* Here compute v particle at time n */
    p[pp].vx = 0.5 * (v12x + v12xp);
    p[pp].vy = 0.5 * (v12y + v12yp);
    p[pp].vz = 0.5 * (v12z + v12zp);


    /*!
      Euler for center of mass 
      \f[
      x_{n+1}=x_{n} + v_{n+1/2} dt
      \f]
    */
    /*!
      Leap-frog for position 
     */
    p[pp].x += (v12xp * dt);
    p[pp].y += (v12yp * dt);
    p[pp].z += (v12zp * dt);
    
    /*
      p[pp].q[0] = cos(0.5*p[pp].theta)*cos(0.5*(p[pp].phi+p[pp].psi));
      p[pp].q[1] = sin(0.5*p[pp].theta)*cos(0.5*(p[pp].phi-p[pp].psi));
      p[pp].q[2] = sin(0.5*p[pp].theta)*sin(0.5*(p[pp].phi-p[pp].psi));
      p[pp].q[3] = cos(0.5*p[pp].theta)*sin(0.5*(p[pp].phi+p[pp].psi));
    */
#ifdef LEAPFROG_ROTATION
    /*! \From here integration as from the book:
         "Computer Simulation of Liquids, M.P.Allen & D.J.Tildesley (OUP, USA, 1988)"
    */

    /*
    fprintf(stderr,"Leapfrog for rotation\n");
    */


    /* ls(t) = ls(t-0.5*dt) + 0.5*dt*ts(t) */
    lx = p[pp].lx + 0.5*dt*p[pp].tx;
    ly = p[pp].ly + 0.5*dt*p[pp].ty;
    lz = p[pp].lz + 0.5*dt*p[pp].tz;

    /* Bulding rotation martix */
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
    lbx = A[0][0] *lx + A[0][1] *ly + A[0][2] *lz;
    lby = A[1][0] *lx + A[1][1] *ly + A[1][2] *lz;
    lbz = A[2][0] *lx + A[2][1] *ly + A[2][2] *lz;

    /* compute angular velocity */
    obx = lbx/p[pp].Ixx;
    oby = lby/p[pp].Iyy;
    obz = lbz/p[pp].Izz;

    /* Build B matrix */
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

    /* q(t+0.5*dt) = q(t) + 0.5*dt*dqdt(t) */
    q0 = p[pp].q[0] + 0.5*dt*dq0dt;
    q1 = p[pp].q[1] + 0.5*dt*dq1dt;
    q2 = p[pp].q[2] + 0.5*dt*dq2dt;
    q3 = p[pp].q[3] + 0.5*dt*dq3dt;

    /* end of the first half step , we go now with the second half */

    /* ls(t+0.5*dt) = ls(t-0.5*dt) + dt*ts(t) , full time step */
    p[pp].lx = p[pp].lx + dt*p[pp].tx;
    p[pp].ly = p[pp].ly + dt*p[pp].ty;
    p[pp].lz = p[pp].lz + dt*p[pp].tz;


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
    p[pp].lbx = A[0][0] *p[pp].lx + A[0][1] *p[pp].ly + A[0][2] *p[pp].lz;
    p[pp].lby = A[1][0] *p[pp].lx + A[1][1] *p[pp].ly + A[1][2] *p[pp].lz;
    p[pp].lbz = A[2][0] *p[pp].lx + A[2][1] *p[pp].ly + A[2][2] *p[pp].lz;

    /* compute angular velocity */
    /*
    p[pp].obx = p[pp].lbx/p[pp].Ixx;
    p[pp].oby = p[pp].lby/p[pp].Iyy;
    p[pp].obz = p[pp].lbz/p[pp].Izz;
    */
    obx12p = p[pp].lbx/p[pp].Ixx;
    oby12p = p[pp].lby/p[pp].Iyy;
    obz12p = p[pp].lbz/p[pp].Izz;

    /* average the angular velocity */
    p[pp].obx = 0.5*(p[pp].obx+obx12p);
    p[pp].oby = 0.5*(p[pp].oby+oby12p);
    p[pp].obz = 0.5*(p[pp].obz+obz12p);

    /* transformation: from omega in the body to omega in the solid frame */
    p[pp].ox = A[0][0]*p[pp].obx + A[1][0]*p[pp].oby + A[2][0]*p[pp].obz; 
    p[pp].oy = A[0][1]*p[pp].obx + A[1][1]*p[pp].oby + A[2][1]*p[pp].obz; 
    p[pp].oz = A[0][2]*p[pp].obx + A[1][2]*p[pp].oby + A[2][2]*p[pp].obz; 




    /* Build B matrix */
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

   /* Now compute dq/dt = 0.5* B*(0,obx,oby,obz) */
    dq0dt = 0.5 *  (B[0][1] * p[pp].obx + B[0][2] * p[pp].oby + B[0][3] * p[pp].obz );
    dq1dt = 0.5 *  (B[1][1] * p[pp].obx + B[1][2] * p[pp].oby + B[1][3] * p[pp].obz );
    dq2dt = 0.5 *  (B[2][1] * p[pp].obx + B[2][2] * p[pp].oby + B[2][3] * p[pp].obz );
    dq3dt = 0.5 *  (B[3][1] * p[pp].obx + B[3][2] * p[pp].oby + B[3][3] * p[pp].obz );

    /* q(t+dt) = q(t) + dt*dqdt(t) */
    p[pp].q[0] = p[pp].q[0] + dt*dq0dt;
    p[pp].q[1] = p[pp].q[1] + dt*dq1dt;
    p[pp].q[2] = p[pp].q[2] + dt*dq2dt;
    p[pp].q[3] = p[pp].q[3] + dt*dq3dt;

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
    
#else

    q0 = p[pp].q[0];
    q1 = p[pp].q[1];
    q2 = p[pp].q[2];
    q3 = p[pp].q[3];

    fprintf(stderr,"q_in %g %g %g %g\n",p[pp].q[0],p[pp].q[1],p[pp].q[2],p[pp].q[3]);

    A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    A[0][1] = 2*(q1*q2+q0*q3);
    A[0][2] = 2*(q1*q3-q0*q2);
    
    A[1][0] = 2*(q1*q2-q0*q3);
    A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A[1][2] = 2*(q2*q3+q0*q1);
    
    A[2][0] = 2*(q1*q3+q0*q2);
    A[2][1] = 2*(q2*q3-q0*q1);
    A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    /* transformation: From torque in the solid frame to torque in the body frame */

    p[pp].tbx = A[0][0] *p[pp].tx + A[0][1] *p[pp].ty + A[0][2] *p[pp].tz;
    p[pp].tby = A[1][0] *p[pp].tx + A[1][1] *p[pp].ty + A[1][2] *p[pp].tz;
    p[pp].tbz = A[2][0] *p[pp].tx + A[2][1] *p[pp].ty + A[2][2] *p[pp].tz;
    
    fprintf(stderr,"t  %g %g %g |%g|\n",p[pp].tx,p[pp].ty,p[pp].tz,
	    p[pp].tx*p[pp].tx+p[pp].ty*p[pp].ty+p[pp].tz*p[pp].tz);
    fprintf(stderr,"tb %g %g %g |%g|\n",p[pp].tbx,p[pp].tby,p[pp].tbz,
	    p[pp].tbx*p[pp].tbx+p[pp].tby*p[pp].tby+p[pp].tbz*p[pp].tbz);

    /* Update omega in body frame  */
#ifdef AB2_NO
    ofx = p[pp].tbx / p[pp].Ixx + (p[pp].Iyy - p[pp].Izz) / p[pp].Ixx * p[pp].oby * p[pp].obz;
    ofy = p[pp].tby / p[pp].Iyy + (p[pp].Izz - p[pp].Ixx) / p[pp].Iyy * p[pp].obz * p[pp].obx;
    ofz = p[pp].tbz / p[pp].Izz + (p[pp].Ixx - p[pp].Iyy) / p[pp].Izz * p[pp].obx * p[pp].oby;

    p[pp].obx += dt*(1.5*ofx - 0.5*p[pp].ofx);
    p[pp].oby += dt*(1.5*ofy - 0.5*p[pp].ofy);
    p[pp].obz += dt*(1.5*ofz - 0.5*p[pp].ofz);
    
    p[pp].ofx = ofx;
    p[pp].ofy = ofy;
    p[pp].ofz = ofz;
    
#else
    p[pp].obx += dt * ( p[pp].tbx / p[pp].Ixx +
                        (p[pp].Iyy - p[pp].Izz) / p[pp].Ixx * p[pp].oby * p[pp].obz);
    p[pp].oby += dt * ( p[pp].tby / p[pp].Iyy +
                        (p[pp].Izz - p[pp].Ixx) / p[pp].Iyy * p[pp].obz * p[pp].obx);
    p[pp].obz += dt * ( p[pp].tbz / p[pp].Izz +
                        (p[pp].Ixx - p[pp].Iyy) / p[pp].Izz * p[pp].obx * p[pp].oby);
#endif

    fprintf(stderr,"o  %g %g %g\n",p[pp].ox,p[pp].oy,p[pp].oz);
    fprintf(stderr,"ob %g %g %g\n",p[pp].obx,p[pp].oby,p[pp].obz);
    
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


#ifdef DEBUG_IMPOSE_ROTATION
    p[pp].obx = 0.0;
    p[pp].oby = 5.e-4;
    p[pp].obz = 0.0;
#endif

    /* transformation: from omega in the body to omega in the solid frame */

    p[pp].ox = A[0][0]*p[pp].obx + A[1][0]*p[pp].oby + A[2][0]*p[pp].obz; 
    p[pp].oy = A[0][1]*p[pp].obx + A[1][1]*p[pp].oby + A[2][1]*p[pp].obz; 
    p[pp].oz = A[0][2]*p[pp].obx + A[1][2]*p[pp].oby + A[2][2]*p[pp].obz; 
    
    /* update quaternion */
#ifdef AB2

    qf0 = 0.5*(B[0][1] * p[pp].obx + B[0][2] * p[pp].oby + B[0][3] * p[pp].obz);
    qf1 = 0.5*(B[1][1] * p[pp].obx + B[1][2] * p[pp].oby + B[1][3] * p[pp].obz);
    qf2 = 0.5*(B[2][1] * p[pp].obx + B[2][2] * p[pp].oby + B[2][3] * p[pp].obz);
    qf3 = 0.5*(B[3][1] * p[pp].obx + B[3][2] * p[pp].oby + B[3][3] * p[pp].obz);

    p[pp].q[0] += dt*(1.5*qf0 - 0.5*p[pp].qf[0]);
    p[pp].q[1] += dt*(1.5*qf1 - 0.5*p[pp].qf[1]);
    p[pp].q[2] += dt*(1.5*qf2 - 0.5*p[pp].qf[2]);
    p[pp].q[3] += dt*(1.5*qf3 - 0.5*p[pp].qf[3]);

    p[pp].qf[0] = qf0;
    p[pp].qf[1] = qf1;
    p[pp].qf[2] = qf2;
    p[pp].qf[3] = qf3;
#else
    p[pp].q[0] += 0.5 * dt * (B[0][1] * p[pp].obx + B[0][2] * p[pp].oby + B[0][3] * p[pp].obz );
    p[pp].q[1] += 0.5 * dt * (B[1][1] * p[pp].obx + B[1][2] * p[pp].oby + B[1][3] * p[pp].obz );
    p[pp].q[2] += 0.5 * dt * (B[2][1] * p[pp].obx + B[2][2] * p[pp].oby + B[2][3] * p[pp].obz );
    p[pp].q[3] += 0.5 * dt * (B[3][1] * p[pp].obx + B[3][2] * p[pp].oby + B[3][3] * p[pp].obz );
#endif

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

#endif
  }
/*end o ffor loop */

#ifdef TIMING
  profile_off(7);
#endif
}
/* end of the function void md advance */











//#define DEBUG_OBJECT

/*! Describe MD object n */
int object(double xl, double yl, double zl, int n)
{
#ifdef DEBUG_OBJECT
  FILE *fout;
#endif
  double d2;
  double drx=(10.*10.);
  double dry=(5.*5.);
  double drz=(10.*10.);

  double x0, y0, z0;
  double dx, dy, dz;
  double d2x1, d2x2;
  double d2y1, d2y2;
  double d2z1, d2z2;

  double theta, phi, psi;
  double q0, q1, q2, q3;
  double A[3][3];
  double x,y,z;
  double x0n, y0n, z0n;

#ifdef TIMING
  profile_on(8);
#endif

  //  fprintf(stderr,"%d %g %g %g\n",n,p[n].x,p[n].y,p[n].z);

  /* remap center of mass on the lattice  */
  x0 = wrap(p[n].x-0.5,SX)+0.5;
  y0 = wrap(p[n].y-0.5,SY)+0.5;
  z0 = wrap(p[n].z-0.5,SZ)+0.5;

  /* each processor treat its own domain  */
  xl += NX*mex;
  yl += NY*mey;
  zl += NZ*mez;

  /*
    q0 = cos(0.5*p[n].theta)*cos(0.5*(p[n].phi+p[n].psi));
    q1 = sin(0.5*p[n].theta)*cos(0.5*(p[n].phi-p[n].psi));
    q2 = sin(0.5*p[n].theta)*sin(0.5*(p[n].phi-p[n].psi));
    q3 = cos(0.5*p[n].theta)*sin(0.5*(p[n].phi+p[n].psi));
  */

  q0 = p[n].q[0];
  q1 = p[n].q[1];
  q2 = p[n].q[2];
  q3 = p[n].q[3];


  // fprintf(stderr,"check q %g %g %g %g - %g\n",q0,q1,q2,q3,q0*q0+q1*q1+q2*q2+q3*q3);
  A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  A[0][1] = 2*(q1*q2+q0*q3);
  A[0][2] = 2*(q1*q3-q0*q2);

  A[1][0] = 2*(q1*q2-q0*q3);
  A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
  A[1][2] = 2*(q2*q3+q0*q1);

  A[2][0] = 2*(q1*q3+q0*q2);
  A[2][1] = 2*(q2*q3-q0*q1);
  A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;

  if (fabs(xl-x0)<SX/2) xl -= x0; else xl = -SX*(xl-x0)/fabs(xl-x0) + xl-x0;
  if (fabs(yl-y0)<SY/2) yl -= y0; else yl = -SY*(yl-y0)/fabs(yl-y0) + yl-y0;
  if (fabs(zl-z0)<SZ/2) zl -= z0; else zl = -SZ*(zl-z0)/fabs(zl-z0) + zl-z0;

  /*
    x = A[0][0] * xl + A[1][0] * yl + A[2][0] * zl;
    y = A[0][1] * xl + A[1][1] * yl + A[2][1] * zl;
    z = A[0][2] * xl + A[1][2] * yl + A[2][2] * zl;
  */

  /*
    Before */
  x = A[0][0] * xl + A[0][1] * yl + A[0][2] * zl;
  y = A[1][0] * xl + A[1][1] * yl + A[1][2] * zl;
  z = A[2][0] * xl + A[2][1] * yl + A[2][2] * zl;
  
  // fprintf(stderr,"%g %g %g -> %g %g %g\n",d2x1,d2y1,d2z1,x,y,z);
  dx = x*x;
  dy = y*y;
  dz = z*z;
  d2 = dx/drx + dy/dry + dz/drz;
 









 
#ifdef DEBUG_OBJECT
  fout = fopen("diagno.dat","a");
  if ((x==32) && (y==32) && (z==32))
    fprintf(fout,"x0 %g %g %g %g %g %g %g\n",x0,d2x1,d2x2,dx,dy,dz,d2);
  fclose(fout);
#endif
#ifdef TIMING
  profile_off(8);
#endif
  
    if (d2<1) return 1; else return 0;
}

/*!  Before this subroutine flag and populations should be in sync
  over processors to make things simple, for the moment this is
  enforced internally
*/
void tag_border_nodes(int *flag, pop_type *f1, pop_type *f2, 
		      poptype *rho1, poptype *rho2, int istep)
{
  int i, j, k, pp, ii, idx[NPOP], count=0;
  double m_x[10], m_y[10], m_z[10];
  double x,y,z, f_x, f_y, f_z;
  double vx, vy, vz, ox, oy, oz;
  double f1_0, f2_0, f1_p, f2_p;
  double xl, yl, zl, x0, y0, z0;
  int particle;
  double *ax=NULL;
  double *ay=NULL;
  double *az=NULL;
  double *tx=NULL;
  double *ty=NULL;
  double *tz=NULL;
  FILE *fout;

#ifdef FLORIAN_TORQUE
  int itx;
  float tmp_tx, tmp_ty,tmp_tz;
#endif


#ifdef TIMING
  profile_on(9);
#endif


  /*! Here I allocate space for NMD particles, if not done already */
  if (ax==NULL) {
    ax = malloc(sizeof(double)*NMD);
    ay = malloc(sizeof(double)*NMD);
    az = malloc(sizeof(double)*NMD);
    tx = malloc(sizeof(double)*NMD);
    ty = malloc(sizeof(double)*NMD);
    tz = malloc(sizeof(double)*NMD);
  }
  
  /*! Put to zero arrays containing MD force and torque on particles */
  for (pp=0;pp<NMD;pp++) {
    ax[pp] = ay[pp] = az[pp] = 0.0;
    tx[pp] = ty[pp] = tz[pp] = 0.0;
    p[pp].ax = p[pp].ay = p[pp].az = 0.0; 
    p[pp].tx = p[pp].ty = p[pp].tz = 0.0;
  }

  
  //  fout = fopen("flag.dat","w");    
  /*! 
    Second we tag inside / outside of particles
    
    Convention is as follow:
    - 0 means fluid node
    - 1 means particles interior
    - (ii>1) means shell of particles (ii-2)
  */
  for (i=1;i<=NX;i++) {
    for (j=1;j<=NY;j++) {
      for (k=1;k<=NZ;k++) {
	ii = IDX(i,j,k);
	flag[ii]=0;
	for (pp=0;pp<NMD;pp++) 
	  if (object(i,j,k,pp)) flag[ii] |= (2+pp); 
	//if (flag[ii]!=0) fprintf(fout,"%d %d %d %d\n",i,j,k,flag[ii]);
      }
    }
  }
  //  fclose(fout);

  /*! Syncronize flag frames */
  pbc_int(flag);

  for (i=1;i<=NX;i++) {
    for (j=1;j<=NY;j++) {
      for (k=1;k<=NZ;k++) {
	idx[0] = IDX(i,j,k);
	/*
	  for (pp=1;pp<NPOP;pp++) {
	  idx[pp] = idx[0] + delta[pp];
	  if  ( (flag[idx[0]] > 1 ) && (flag[idx[pp]] == 0 ) ) {
	    if (border_nodes >= border_idx_size) {
	      border_idx_size += 100;
	      border_idx = 
		(int *) realloc(border_idx,sizeof(int)*border_idx_size);
	    }
	    border_idx[border_nodes] = idx[0];
	    border_nodes++;
	  }
	}
	*/

	/*!  Here we detect if we are in the interior of a particle.
	  I am not sure this will be necessary in the future...
	  Indeed probably all can be done with the logic of
	  border_nodes
	*/
	if ( (flag[idx[0]] > 0 ) &&
	     (flag[idx[1]] > 0 ) &&
	     (flag[idx[2]] > 0 ) &&
	     (flag[idx[3]] > 0 ) &&
	     (flag[idx[4]] > 0 ) &&
	     (flag[idx[5]] > 0 ) &&
	     (flag[idx[6]] > 0 ) &&
	     (flag[idx[7]] > 0 ) &&
	     (flag[idx[8]] > 0 ) &&
	     (flag[idx[9]] > 0 ) &&
	     (flag[idx[10]] > 0 ) &&
	     (flag[idx[11]] > 0 ) &&
	     (flag[idx[12]] > 0 ) &&
	     (flag[idx[13]] > 0 ) &&
	     (flag[idx[14]] > 0 ) &&
	     (flag[idx[15]] > 0 ) &&
	     (flag[idx[16]] > 0 ) &&
	     (flag[idx[17]] > 0 ) &&
	     (flag[idx[18]] > 0 ) ) {
	  flag[idx[0]] = 1; 
	  /*! 
	     Putting to zero the populations inside the object is a
	     way to test the correctness of boundary conditions but it
	     would skrew up the momentum balance, beware !
	     
	     \code 
	     for (pp=0;pp<NPOP;pp++) {
	     f1[idx[0]].p[pp] = 0.0;
	     f2[idx[0]].p[pp] = 0.0;
	     }
	     \endcode
	  */
	}
      }
    }
  }

  fout = fopen("border.dat","w");
  for (i=0;i<border_nodes;i++)
    fprintf(fout,"%d %d\n",i,border_idx[i]);
  fclose(fout);
  
  /*! Syncronize flag frames */
  pbc_int(flag);

  /*! \bug Not sure, this should be set in main, before entering here */
  /*
    pbc(f1);
    pbc(f2);
  */
  fout = fopen("shape.dat","w");

  for (i=1;i<=NX;i++) {
    for (j=1;j<=NY;j++) {
      for (k=1;k<=NZ;k++) {
	idx[0] = IDX(i,j,k);
	for (pp=1;pp<NPOP;pp++) {
	  idx[pp] = idx[0] + delta[pp];
	}
	
	/*! put forces to zero: f_x = f_y = f_z = 0.0; */
	for (pp=1;pp<NPOP;pp++) {
	 
	  /* Notice these are adjusted prestreaming */
	  if  ( (flag[idx[0]] > 1 ) && (flag[idx[pp]] == 0 ) ) {

	    
	    particle = flag[idx[0]]-2;

	    //fprintf(stderr,"Part = %d\n",particle);
	    /*! This was the wrong way to do it */
	    x = (physx(i) + 0.5 * cx[pp] - p[particle].x);
	    y = (physy(j) + 0.5 * cy[pp] - p[particle].y);
	    z = (physz(k) + 0.5 * cz[pp] - p[particle].z);
	    
	    xl = (physx(i) + 0.5 * cx[pp]);
	    yl = (physy(j) + 0.5 * cy[pp]);
	    zl = (physz(k) + 0.5 * cz[pp]);
	    x0 = wrap(p[particle].x-0.5,SX)+0.5;
	    y0 = wrap(p[particle].y-0.5,SY)+0.5;
	    z0 = wrap(p[particle].z-0.5,SZ)+0.5;
	      
	    if (fabs(xl-x0)<SX/2) xl -= x0; 
	    else xl = -SX*(xl-x0)/fabs(xl-x0) + xl-x0;
	    if (fabs(yl-y0)<SY/2) yl -= y0; 
	    else yl = -SY*(yl-y0)/fabs(yl-y0) + yl-y0;
	    if (fabs(zl-z0)<SZ/2) zl -= z0; 
	    else zl = -SZ*(zl-z0)/fabs(zl-z0) + zl-z0;
	    
	    x = xl;
	    y = yl;
	    z = zl;
	    
	    /*! Velocity components at the interface */
	    vx = p[particle].vx - (y*p[particle].oz - z*p[particle].oy);
	    vy = p[particle].vy - (z*p[particle].ox - x*p[particle].oz);
	    vz = p[particle].vz - (x*p[particle].oy - y*p[particle].ox);

	    /*! Here we take into account the exchanges of momentum */
	    f1_0 = f1[idx[0]].p[pp];
	    f2_0 = f2[idx[0]].p[pp];
	    
	    f1_p = f1[idx[pp]].p[inv[pp]];
	    f2_p = f2[idx[pp]].p[inv[pp]];

	    ax[particle] += -2.*cx[pp]*(f1[idx[pp]].p[inv[pp]] + 
				   invcs2*rho1[idx[pp]]*ww[pp] * 
				   (vx*cx[pp]+vy*cy[pp]+vz*cz[pp]));
	    ay[particle] += -2.*cy[pp]*(f1[idx[pp]].p[inv[pp]] + 
				   invcs2*rho1[idx[pp]]*ww[pp] * 
				   (vx*cx[pp]+vy*cy[pp]+vz*cz[pp]));
	    az[particle] += -2.*cz[pp]*(f1[idx[pp]].p[inv[pp]] + 
				   invcs2*rho1[idx[pp]]*ww[pp] * 
				   (vx*cx[pp]+vy*cy[pp]+vz*cz[pp]));
	    
	    ax[particle] += -2.*cx[pp]*(f2[idx[pp]].p[inv[pp]] + 
				   invcs2*rho2[idx[pp]]*ww[pp] * 
				   (vx*cx[pp]+vy*cy[pp]+vz*cz[pp]));
	    ay[particle] += -2.*cy[pp]*(f2[idx[pp]].p[inv[pp]] + 
				   invcs2*rho2[idx[pp]]*ww[pp] * 
				   (vx*cx[pp]+vy*cy[pp]+vz*cz[pp]));
	    az[particle] += -2.*cz[pp]*(f2[idx[pp]].p[inv[pp]] + 
				   invcs2*rho2[idx[pp]]*ww[pp] * 
				   (vx*cx[pp]+vy*cy[pp]+vz*cz[pp]));

	    /*! Adjust the population entering into the solid particle */
	    f1[idx[pp]].p[inv[pp]] = f1_0 - 2.0*invcs2*rho1[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp]);
	    f2[idx[pp]].p[inv[pp]] = f2_0 - 2.0*invcs2*rho2[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp]);

	    /*! Here adjust the momentum going out of the solid particles */
	    f1[idx[0]].p[pp] = f1_p + 2.0*invcs2*rho1[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp]);
	    f2[idx[0]].p[pp] = f2_p + 2.0*invcs2*rho2[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp]);

	    
	    f_x = 2.*cx[pp]*( f1_0 - f1_p - 
			      2.0*invcs2*rho1[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp])
			      +
			      f2_0 - f2_p - 
			      2.0*invcs2*rho2[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp])
			      );
	    f_y = 2.*cy[pp]*( f1_0 - f1_p - 
			      2.0*invcs2*rho1[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp])
			      +
			      f2_0 - f2_p - 
			      2.0*invcs2*rho2[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp])
			      );
	    f_z = 2.*cz[pp]*( f1_0 - f1_p - 
			      2.0*invcs2*rho1[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp])
			      +
			      f2_0 - f2_p - 
			      2.0*invcs2*rho2[idx[pp]]*ww[pp]*(vx*cx[pp]+vy*cy[pp]+vz*cz[pp])
			      );

	    /*! Here cumulate the torque */
	    tx[particle] += (y*f_z - z*f_y);
	    ty[particle] += (z*f_x - x*f_z);
	    tz[particle] += (x*f_y - y*f_x);


	    // fprintf(stderr,"%g %g %g\n",f_x,f_y,f_z);
	    fprintf(fout,"%g %g %g %g %g %g %d\n",x,y,z,vx,vy,vz,particle);
	    
	    count++;
	  } /* if */
	}
      }
    }
  }
  
  fclose(fout);


  /*! The call to invpbc ensure that modified populations on the frame
      are synced back on the appropriate lattice position */

  invpbc(f1);
  invpbc(f2);

  
  /*!  MPI_Allreduce ensure that contributions to acceleration and
      torque from different processors are all summed up and visible
      to all processors */
  for (pp=0;pp<NMD;pp++) {
    MPI_Allreduce(tx+pp,&(p[pp].tx),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(ty+pp,&(p[pp].ty),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(tz+pp,&(p[pp].tz),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(ax+pp,&(p[pp].ax),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(ay+pp,&(p[pp].ay),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(az+pp,&(p[pp].az),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }


#ifdef FLORIAN_TORQUE
  tmp_tx = tmp_ty = tmp_tz = 0.0;
  pp=0; /* for the first particle? */
  fout = fopen("torque.in","r");
  for(itx=0; itx<=istep; itx++){
    fscanf(fout,"%e %e %e",&tmp_tx,&tmp_ty,&tmp_tz);
  }
  fprintf(stderr,"istep %d Florian torque %e %e %e\n",istep, tmp_tx, tmp_ty, tmp_tz);
  fclose(fout);
  p[pp].tx = (double)tmp_tx;
  p[pp].ty = (double)tmp_ty;
  p[pp].tz = (double)tmp_tz;

  tmp_tx = tmp_ty = tmp_tz = 0.0;
  fout = fopen("force.in","r");
  for(itx=0; itx<=istep; itx++){
    fscanf(fout,"%e %e %e",&tmp_tx,&tmp_ty,&tmp_tz);
  }
  fprintf(stderr,"istep %d Florian force %e %e %e\n",istep, tmp_tx, tmp_ty, tmp_tz);
  fclose(fout);
  p[pp].ax = (double)tmp_tx;
  p[pp].ay = (double)tmp_ty;
  p[pp].az = (double)tmp_tz;
  /*
  for (pp=0;pp<NMD;pp++) {
    p[pp].tx = 0.0;
    p[pp].ty = 0.0;
    p[pp].tz = 0.0;
    p[pp].ax = 0.0;
    p[pp].ay = 0.0;
    p[pp].az = 0.0;
  }
  */
#endif




  /* if (istep==2) exit(99); */

  if (AMIROOT) {
    for (pp=0;pp<NMD;pp++) {
      
#ifdef DEBUG
      fprintf(stderr,"Solid-fluid Links: %d momentum %g %g %g  %g %g %g a %g %g %g t %g %g %g o %g %g %g links %d\n",count,
	      p[pp].x,p[pp].y,p[pp].z,
	      p[pp].vx,p[pp].vy,p[pp].vz,
	      p[pp].ax,p[pp].ay,p[pp].az,
	      p[pp].tx,p[pp].ty,p[pp].tz,
	      p[pp].ox,p[pp].oy,p[pp].oz,
	      count);
#endif
      
      fout = fopen("p.dat","a");
      fprintf(fout,"%d x %g y %g z %g vx %g vy %g vz %g ax %g ay %g az %g torque  %g %g %g o %g %g %g ob %g %g %g links %d\n",istep,
	      p[pp].x,p[pp].y,p[pp].z,
	      p[pp].vx,p[pp].vy,p[pp].vz,
	      p[pp].ax,p[pp].ay,p[pp].az,
	      p[pp].tx,p[pp].ty,p[pp].tz,
	      p[pp].ox,p[pp].oy,p[pp].oz,
	      p[pp].obx,p[pp].oby,p[pp].obz,
	      count);
      fclose(fout);
    }
  }

#ifdef TIMING
  profile_off(9);
#endif

}


