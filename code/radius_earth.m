function r_earth = radius_earth(a,c,rr_sc)
%Computes the distance from the center to the surface of the Earth considering 
%it an ellipsoid (more precisely a spheroid) defined by a and b, based on
%the position of a satellite in an orbit around it.
%
% PROTOTYPE r_earth = radius_earth(a,c,rr_sc)
% 
% INPUT
% a[1]         Semi-major axis length
% c[1]         Semi-minor axis length
% 
% rr_sc[3x1]   Position of the spacecraft in Cartesian coordinates (x,y,z)  
% 
% OUTPUT
% r_earth[3x1] Position of a point in the Earth's surface in Cartesian
%              coordinates (x,y,z)


%% Compute the polar angle (phi) and azimuth angle (lambda) of the ellipsoid
% from the satellite position rr_sc

r_sc = norm(rr_sc);
phi = asin(rr_sc(3)/r_sc);
lambda = atan(rr_sc(2)/rr_sc(1));

%% Equations of the ellipsoid's surface points
% The equations are:
% x = r_earth*cos(phi)*cos(lambda);
% y = r_earth*cos(phi)*sin(lambda);
% z = r_earth*sin(phi);
% x^2/a^2+y^2/b^2+z^2/c^2=1

b = a; %for a spheroid (two sides equal)

termx = cos(phi)^2*cos(lambda)^2/a^2;
termy = cos(phi)^2*sin(lambda)^2/b^2;
termz = sin(phi)^2/c^2;
r_earth = sqrt(1/(termx+termy+termz));

end

