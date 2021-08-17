function [xx, tt] = GOCE_int(t_ode, data, cond0, options, fail_T, fail_valve, odeselect)

% Sets the integration for the given values, iterates in case of actuation 
% of any event.
%
% INPUT:
% t_ode   [vector]          Vector of integration time.
% data    [struct]          Structure with all the data necessary.
% cond0   [vector]          Vector of initial variables.
% options [struct]          Set of options for the ODE.
% fail_T  [bool]            Boolean with value 0 if there is no engine
%                           failure, 1 otherwise.
% fail_valv [bool]          Boolean with value 0 if there is no valve
%                           failure, 1 otherwise.
% odeselect [constant]      Constant with value 1, 2 or 3 if the desired           
%                           integrator is ode23s, ode15s or ode23tb.
%
% OUTPUT:
% xx [matrix]               Matrix of the simulated variables.
% tt [vector]               Vector of the simulated time.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Authors: Alejandro de Miguel, Rafael Felix, Carmen Salas
% Last modification: 15/01/2020
% Politecnico di Milano, Modeling and Simulation of Aerospace Systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 6
    odeselect = 2;
end

t0     = t_ode(1);                                                          % t0 and tf defined such that t_ode can be defined as
tf     = t_ode(end);                                                        % initial and final times and as a time vector
t_int  = t_ode;                                                             % Time integrated
cont_t = 1;

while t0 < tf

    if odeselect == 1
    [t, aux ,~ ,xe , ie] = ode23s(@(tt, xx) GOCE_ode(tt, xx, data, ...
                           fail_T, fail_valve), t_int, cond0, options);
    end 
    if odeselect == 2
    [t, aux ,~ ,xe , ie] = ode15s(@(tt, xx) GOCE_ode(tt, xx, data, ...
                           fail_T, fail_valve), t_int, cond0, options);
    end
    if odeselect == 3
    [t, aux ,~ ,xe , ie] = ode23tb(@(tt, xx) GOCE_ode(tt, xx, data, ...
                           fail_T, fail_valve), t_int, cond0, options);
    end    
    
    
    xx(cont_t:cont_t+length(t)-1,:) = aux;                                  % Solutions are saved
    tt(cont_t:cont_t+length(t)-1)   = t;                                    % Time is saved
    
    i1 = find(ie == 1 | ie == 2);                                           % Index of ie if the accelerometer reaches the limits
    i2 = find(ie == 3 | ie == 4);                                           % Index of ie if the valve reaches the limits
    
    if     ie(i1) == 1
        xe(7) = data.g;                                                     % Accelerometer has reached a limit, velocity is reset
        xe(8) = 0;                                                          % to zero, position to the limit 
    elseif ie(i1) == 2                                                                 
        xe(7) = -data.g;
        xe(8) = 0;
    end
    
    if     ie(i2) == 3
        xe(10) = 10*data.A0;                                                % Valve has reached a limit, velocity is reset to
        xe(11) = 0;                                                         % zero, position to the limit
    elseif ie(i2) == 4                                                          
        xe(10) = 0;
        xe(11) = 0;
    end
    
    cond0  = xe;                                                            % Corrected solutions at the time of the event are
    cont_t = cont_t + length(t) - 1;                                        % used as next initial conditions
    t0     = t(end);
    
    if length(t_ode) > 2                                                    % Time of integration is redefined
        t_int = t_ode(cont_t:end);
    else
        t_int = [t0 tf];
    end
    
end
end