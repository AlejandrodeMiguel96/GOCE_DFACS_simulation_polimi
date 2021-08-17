function [value, isterminal ,direction] = GOCE_event(~ , xx)

% Function that actuates inside the ode in case several limit situations 
% occur.
%
% INPUT:
% xx    [vector 1x12]       Vector of variables.
%
% OUTPUT:
% value [vector 1x4]        Vector of boolean values indicating whether or
%                           not an event has occured.
% isterminal [vector 1x4]   Vector indicating whether integration must stop 
%                           should an event occur.
% direction  [vector 1x4]   Vector indicating when to take an event into 
%                           consideration.
%
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Authors: Alejandro de Miguel, Rafael Felix, Carmen Salas
% Last modification: 09/12/2019
% Politecnico di Milano, Modeling and Simulation of Aerospace Systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

value = [1; 1; 1; 1];

g  = 5e-4;                                                                  % Electrodes-Mass Gap [m]
A0 = 1e-5;                                                                  % Valve Orifice Area [m^2]

if     xx(7)  >=  g && xx(8) > 0                                            % Accelerometer's upper limit
    value(1) = 0;
elseif xx(7)  <= -g && xx(8) < 0                                            % Accelerometer's lower limit
    value(2) = 0;
end

if     xx(10) >= 10*A0 && xx(11) > 0                                        % Valve's upper limit
    value(3) = 0;
elseif xx(10) <= 0     && xx(11) < 0                                        % Valve's lower limit
    value(4) = 0;
end

isterminal = [1; 1; 1; 1];                                                  % Integration is stopped should any condition be reached
direction  = [0; 0; 0; 0];

end