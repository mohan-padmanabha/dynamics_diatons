/*initial parameters */
a =  ;
b =  ;
rho_p = ;
pi = 3.1415;
beta = a/b ;
mu = ;
mp = 4/3*pi*pow(a,3)*beta*rho_p ;

/*inertial forces*/

I_xx = ((1+pow(beta,2)*pow(a,2))/5*mp ;
I_yy = I_xx ;
I_zz = (2*pow(a,2))/5*mp ;

/* hydrodynamic forceses*/
K_xx = (16(pow(beta,2) - 1))/((2*pow(beta,2)-3) ln(beta + sqrt(pow(beta,2)-1)/(sqrt(pow(beta,2) - 1))+ beta);

K_yy = K_xx;

K_zz = (8(pow(beta,2) - 1))/((2*pow(beta,2)-1) ln(beta + sqrt(pow(beta,2)-1)/(sqrt(pow(beta,2) - 1))- beta);

/*rotation matrix*/

quaternion euler2quaternion(euler e)
(
  quaternion q;
  q.q[0] = cos(0.5*e.theta)*cos(0.5*(e.phi+e.psi));
  q.q[1] = sin(0.5*e.theta)*cos(0.5*(e.phi-e.psi));
  q.q[2] = sin(0.5*e.theta)*sin(0.5*(e.phi-e.psi));
  q.q[3] = cos(0.5*e.theta)*sin(0.5*(e.phi+e.psi));
  return q;
)
    A[0][0] = 1-2(pow(q2,2) - pow(q3,2));
    A[0][1] = 2*(q1*q2+q0*q3);
    A[0][2] = 2*(q1*q3-q0*q2);
    
    A[1][0] = 2*(q1*q2-q0*q3);
    A[1][1] = 1-2(pow(q1,2) + pow(q3,2));
    A[1][2] = 2*(q2*q3+q0*q1);
    
    A[2][0] = 2*(q1*q3+q0*q2);
    A[2][1] = 2*(q2*q3-q0*q1);
    A[2][2] = 1-2(pow(q1,2) + pow(q2,2));


alfa_0 = (pow(beta,2)/(pow(beta,2) - 1)) + beta/(2*pow(pow(beta,2)-1),3/2))* ln ((beta - sqrt(pow(beta,2) - 1))/(beta + sqrt(pow(beta,2)-1)));
beta_0 = alfa_0 ;
gama_0 = - (2/(pow(beta,2) - 1)) - beta/(pow(pow(beta,2)-1),3/2))* ln ((beta - sqrt(pow(beta,2) - 1))/(beta + sqrt(pow(beta,2)-1))) ;


S_zy = 1/2*(dy_uz + dz_uy) ;
S_xz = 1/2*(dz_ux + dx_uz) ;

W_zy = 1/2*(dy_uz - dz_uy) ;
W_xz = 1/2*(dz_ux - dx_uz) ;
W_yx = 1/2*(dx_uy + dy_ux) ;



/*Hydrodynamic torque (jeffery 1922 *for linear shear flow ) */

T_x = (16*pi*mu*pow(a,3)*beta)/(3*(beta_0 + pow(beta,2)*gama_0))*((1-pow(beta,2))*S_zy + (1+ pow(beta,2))*(W_yz - W_x)) ;

T_y = (16*pi*mu*pow(a,3)*beta)/(3*(pow(beta,2)*gama_0 + alfa_0))*((pow(beta,2) - 1)*S_xz + (pow(beta,2)+1)*(W_xz - W_y)) ;

T_z = (32*pi*mu*pow(a,3)*beta)/(3*(alfa_0 + gama_0))*(W_yx-W_z) ;

F_x = mu A[][] * K_xx * A_t[][] * (ux_f - vx_p)
F_y = mu A[][] * K_yy * A_t[][] * (uy_f - vy_p)
F_z = mu A[][] * K_zz * A_t[][] * (uz_f - vz_p)


dt_upx = F_x / m
dt_upy = F_y / m
dt_upz = F_z / m

dt_Wx = (T_x + W_y*W_z*(I_y-I_z))/I_x
dt_Wy = (T_y + W_z*W_x*(I_z-I_x))/I_y
dt_Wz = (T_z + W_x*W_y*(I_x-I_y))/I_z





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
     /* my part of the program */ 
     /* needs a lot of changes*/
     /* has to be optimized */
     /* here we compute the force required for translation of the element */
    A[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    A[0][1] = 2*(q1*q2+q0*q3);
    A[0][2] = 2*(q1*q3-q0*q2);
    
    A[1][0] = 2*(q1*q2-q0*q3);
    A[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    A[1][2] = 2*(q2*q3+q0*q1);
    
    A[2][0] = 2*(q1*q3+q0*q2);
    A[2][1] = 2*(q2*q3-q0*q1);
    A[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
     


    beta = b/a ;
    mp = 4/3*pi*pow(a,3)*beta*rho_p ;

    K_xx = (16(pow(beta,2) - 1))/((2*pow(beta,2)-3) ln(beta + sqrt(pow(beta,2)-1)/(sqrt(pow(beta,2) - 1))+ beta);

    K_yy = K_xx;

    K_zz = (8(pow(beta,2) - 1))/((2*pow(beta,2)-1) ln(beta + sqrt(pow(beta,2)-1)/(sqrt(pow(beta,2) - 1))- beta);



/* translation dyadic or resistance tensor */
    K_tranx = (A[0][0] * K_xx)* A[0][0]+ (A[0][1] * K_xx)* A[0][1]+(A[0][2] * K_xx)* A[0][2]

    K_trany = (A[1][0] * K_yy )* A[1][0]+ ( A[1][1] * K_yy)* A[1][1]+(A[1][2] * K_yy)* A[1][2]

    K_tranz = (A[2][0] * K_zz)* A[2][0]+ ( A[2][1] * K_zz)* A[2][1]+(A[2][2] * K_zz)* A[2][2]

/* hydrodynamic drag force (brenner 1964) */

    fx = mu K_tranx * (ux_f - vx_p);
    fy = mu K_trany * (uy_f - vy_p);
    fz = mu K_tranz * (uz_f - vz_p);   



    /* Here compute v particle at time n-1/2 */
    v12x =  0.5 * (p[pp].vx + p[pp].vx0);
    v12y =  0.5 * (p[pp].vy + p[pp].vy0);
    v12z =  0.5 * (p[pp].vz + p[pp].vz0);

    p[pp].vx0 = p[pp].vx;
    p[pp].vy0 = p[pp].vy;
    p[pp].vz0 = p[pp].vz;

    /* Here compute v particle at time n+1/2 */

    v12xp = v12x + (p[pp].fx * dt) / p[pp].m;
    v12yp = v12y + (p[pp].fy * dt) / p[pp].m;
    v12zp = v12z + (p[pp].fz * dt) / p[pp].m;

    /* Here compute v particle at time n */
    /* for now we are not using this as the value has to 
      be extrapolated to find the value of force 
    p[pp].vx = 0.5 * (v12x + v12xp);
    p[pp].vy = 0.5 * (v12y + v12yp);
    p[pp].vz = 0.5 * (v12z + v12zp);
    */

    /* extrapolated velocity value to 
    to find the value of force at time n +1 */ 
    p[pp].vxp = 0.5 * (v12x + 3*v12xp);
    p[pp].vyp = 0.5 * (v12y + 3*v12yp);
    p[pp].vzp = 0.5 * (v12z + 3*v12zp);



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


#ifdef LEAPFROG_ROTATION
    /*Hydrodynamic torque (jeffery 1922 *for linear shear flow ) */

    I_xx = ((1+pow(beta,2)*pow(a,2))/5*mp ;
    I_yy = I_xx ;
    I_zz = (2*pow(a,2))/5*mp ;


    S_zy = 1/2*(dy_uz + dz_uy) ;
    S_xz = 1/2*(dz_ux + dx_uz) ;



    W_zy = 1/2*(dy_uz - dz_uy) ;
    W_xz = 1/2*(dz_ux - dx_uz) ;
    W_yx = 1/2*(dx_uy + dy_ux) ;



    alfa_0 = (pow(beta,2)/(pow(beta,2) - 1)) + beta/(2*pow(pow(beta,2)-1),3/2))* ln ((beta - sqrt(pow(beta,2) - 1))/(beta + sqrt(pow(beta,2)-1)));
    beta_0 = alfa_0 ;
    gama_0 = - (2/(pow(beta,2) - 1)) - beta/(pow(pow(beta,2)-1),3/2))* ln ((beta - sqrt(pow(beta,2) - 1))/(beta + sqrt(pow(beta,2)-1))) ;



    T_x = (16*pi*mu*pow(a,3)*beta)/(3*(beta_0 + pow(beta,2)*gama_0))*((1-pow(beta,2))*S_zy + (1+ pow(beta,2))*(W_yz - W_x)) ;

    T_y = (16*pi*mu*pow(a,3)*beta)/(3*(pow(beta,2)*gama_0 + alfa_0))*((pow(beta,2) - 1)*S_xz + (pow(beta,2)+1)*(W_xz - W_y)) ;

    T_z = (32*pi*mu*pow(a,3)*beta)/(3*(alfa_0 + gama_0))*(W_yx-W_z) ;




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

    ox = o12x + (doxdt * (dt*0.5));
    oy = o12y + (doydt * (dt*0.5));
    oz = o12z + (dozdt * (dt*0.5));   


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
    obx = A[0][0] *ox + A[0][1] *oy + A[0][2] *oz;
    oby = A[1][0] *ox + A[1][1] *oy + A[1][2] *oz;
    obz = A[2][0] *ox + A[2][1] *oy + A[2][2] *oz;

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



    ox12p = o12x + doxdt *dt;
    oy12p = o12y + doydt *dt;
    oz12p = o12z + dozdt *dt;

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
    obx12p = A[0][0] *ox12p + A[0][1] *oy12p + A[0][2] *oz12p;
    oby12p = A[1][0] *ox12p + A[1][1] *oy12p + A[1][2] *oz12p;
    obz12p = A[2][0] *ox12p + A[2][1] *oy12p + A[2][2] *oz12p;

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
    dq0dt = 0.5 *  (B[0][1] * obx12p + B[0][2] * oby12p + B[0][3] * obz12p );
    dq1dt = 0.5 *  (B[1][1] * obx12p + B[1][2] * oby12p + B[1][3] * obz12p );
    dq2dt = 0.5 *  (B[2][1] * obx12p + B[2][2] * oby12p + B[2][3] * obz12p );
    dq3dt = 0.5 *  (B[3][1] * obx12p + B[3][2] * oby12p + B[3][3] * obz12p );

    q0 = p[pp].q[0] + 0.5*dt*dq0dt;
    q1 = p[pp].q[1] + 0.5*dt*dq1dt;
    q2 = p[pp].q[2] + 0.5*dt*dq2dt;
    q3 = p[pp].q[3] + 0.5*dt*dq3dt;

    opx = 0.5 * (o12x + 3*o12xp);
    opy = 0.5 * (o12y + 3*o12yp);
    opz = 0.5 * (o12z + 3*o12zp);
    




