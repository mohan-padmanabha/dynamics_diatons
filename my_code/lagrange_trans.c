

    my_double pi;
    my_double q0, q1, q2, q3;
    my_double matQ[2][2];
    my_double aspr, mp, a, rhop, mu;
    my_double k_xx, k_yy, k_zz;
    my_double k_tranx, k_trany, k_tranz;
    my_double fx, fy, fz;
    my_double v12x, v12y, v12z;



    pi = 3.1415;
    
    (tracer+ipart)->q[0] = q0;
    (tracer+ipart)->q[1] = q1;
    (tracer+ipart)->q[2] = q2;
    (tracer+ipart)->q[3] = q3;
    



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





