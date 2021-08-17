function cond0 = GOCE_initcond(data)

% Calculates initial values for the integration starting from an
% equilibrium point where T=D.
%
% INPUT:
% data [struct]             Structure with all the data necessary.
%
% OUTPUT:
% cond0 [vector 1x12]       Vector of initial variables.
%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Authors: Alejandro de Miguel, Rafael Felix, Carmen Salas
% Last modification: 09/12/2019
% Politecnico di Milano, Modeling and Simulation of Aerospace Systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rr0,vv0] = kep2car(data.kep0, data.muE);                                   % [km, km/s]
vv0_rel   = vv0 - cross(data.ww_earth, rr0);
v0_rel    = norm(vv0_rel);
rho       = atm_model(data.h0);
A_sc      = data.A*1e-6;                                                    % Frontal area in [km^2]
D0        = 1e3*0.5*A_sc*data.cd*rho*v0_rel^2*vv0_rel/v0_rel;               % Drag force modulus [N]
D0_sc     = dot(D0,vv0)/norm(vv0);
T0        = D0_sc;                                                          % Thrust force modulus [N]
mdot0     = T0/sqrt(2*data.e*data.DV/data.m_i);   
A_valv0   = mdot0/(data.cdis*sqrt(2/data.R/data.T2*data.k/(data.k - 1))*... % Valve area [m^2]
            data.p2*sqrt((2/(data.k + 1))^(2/(data.k - 1)) - (2/(data.k ...
            + 1))^((data.k + 1)/(data.k - 1))));
options   = optimset('Display', 'off', 'TolFun', 1e-9);
alpha0    = fzero(@(alpha) alpha - sin(alpha) - 2*pi*A_valv0/data.A0, ...   % Angle defining geometry
          3*pi/2, options);
xvalv0    = 10*data.A0*cos(alpha0/4)^2;
I0        = data.Kfcv*(10*data.A0 - xvalv0)/data.Ki;

% Equilibrium assumption

xa0    = 0;
va0    = 0;
vvalv0 = 0; 
Vout0  = 0;

cond0 = [data.kep0, xa0, va0, Vout0, xvalv0, vvalv0, I0]';

end