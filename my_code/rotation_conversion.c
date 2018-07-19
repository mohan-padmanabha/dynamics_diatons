

/* conversion from quaternions to axis angles */
theta = 2 * acos(q0) ;
x = q1 / sqrt(1-q0*q0) ;
y = q2 / sqrt(1-q0*q0) ;
z = q3 / sqrt(1-q0*q0) ;

/* conversion from axis angles to quaternions */
q1 = x * sin(theta/2);
q2 = y * sin(theta/2);
q3 = z * sin(theta/2);
q0 = cos(theta/2);


/* conversion from quaternions to euler angles */
heading = atan2(2*q2*q0-2*q1*q3 , 1 - 2*q2*q2 - 2*q3*q3);
attitude = asin(2*q1*q2 + 2*q3*q0) ;
bank = atan2(2*q1*q0-2*q2*q3 , 1 - 2*q1*q1 - 2*q3*q3);


/*
except when qx*qy + qz*qw = 0.5 (north pole)
which gives:
heading = 2 * atan2(x,w)
bank = 0
and when qx*qy + qz*qw = -0.5 (south pole)
which gives:
heading = -2 * atan2(x,w)
bank = 0

*/

/* conversion from euler angles to quaternions */

c1 = cos(theta / 2) ;
c2 = cos(phi / 2) ;
c3 = cos(psi / 2) ;
s1 = sin(theta / 2) ;
s2 = sin(phi / 2) ;
s3 = sin(psi / 2) ;

q0 = c1 c2 c3 - s1 s2 s3 ;
q1 = s1 s2 c3 + c1 c2 s3 ;
q2 = s1 c2 c3 + c1 s2 s3 ;
q3 = c1 s2 c3 - s1 c2 s3 ;





