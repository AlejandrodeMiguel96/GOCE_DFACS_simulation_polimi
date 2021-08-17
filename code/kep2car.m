function [ rr_GE, vv_GE ] = kep2car(kep,mu)
% Transformation from Keplerian elements to cartesian coordinates
%
% [ rr, vv ] = kep2car( a, e, i, OM, om, th, mu )
%
% -------------------------------------------------------------------------
% Input arguments:
% a [1x1] semi-major axis [km]
% e [1x1] eccentricity [-]
% i [1x1] inclination [rad]
% OM [1x1] RAAN (Right Ascension of the Ascending Node) [rad]
% om [1x1] argument of periapsis [rad]
% th [1x1] true anomaly [rad]
% mu [1x1] gravitational parameter [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% rr [3x1] position vector [km]
% vv [3x1] velocity vector [km/s]

a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
th = kep(6);

R3_OM=[cos(OM) sin(OM) 0
    -sin(OM) cos(OM) 0
    0 0 1];
R1_i=[1 0 0
    0 cos(i) sin(i)
    0 -sin(i) cos(i)];
R3_om=[cos(om) sin(om) 0
    -sin(om) cos(om) 0
    0 0 1];

p=a*(1-e^2); %semi-latus rectum
r=p/(1+e*cos(th)); 

rr_PF=r*[cos(th); sin(th); 0];
vv_PF=sqrt(mu/p)*[-sin(th); e+cos(th);0];

T=(R3_OM)'*(R1_i)'*(R3_om)';

rr_GE=T*rr_PF;
vv_GE=T*vv_PF;