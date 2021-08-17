function GOCE_plots(t, df, data)

% Plots the evolution of keplerian elements, accelerometer and valve 
% parameters, thrust and drag for the given simulation.
%
% INPUT:
% t    [vector]             Vector of simulated time.
% xx   [matrix]             Matrix of simulated variables.
% data [struct]             Structure with all the data necessary.
%
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Authors: Alejandro de Miguel, Rafael Felix, Carmen Salas
% Last modification: 15/01/2020
% Politecnico di Milano, Modeling and Simulation of Aerospace Systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rr = zeros(length(t),3);
for i = 1:length(t)
    rr(i,:) = kep2car(df(i,:),data.muE);
end

figure(1)
plot(t, df(:,1))
xlabel('Time $[s]$'    , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Semimajor axis', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(2)
plot(t, df(:,2))
xlabel('Time $[s]$'  , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Eccentricity', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(3)
plot(t, df(:,3))
xlabel('Time $[s]$' , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Inclination', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(4)
plot(t, df(:,4))
xlabel('Time $[s]$'    , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Ascending node', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(5)
plot(t, df(:,5))
xlabel('Time $[s]$'            , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Argument of periapsis' , 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(6)
plot(t, df(:,6))
xlabel('Time $[s]$'  , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('True anomaly', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(7) %orbit plot 3D with perturbations
plot3(rr(:,1)  ,rr(:,2)  ,rr(:,3));
hold on
plot3(rr(1,1)  ,rr(1,2)  ,rr(1,3)  ,'g*');
plot3(rr(end,1),rr(end,2),rr(end,3),'r*')
grid on

figure(9)
plot(t, df(:,7))
xlabel('Time $[s]$', 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Accelerometer''s position $[m]$', 'FontSize', 12, 'Interpreter',...
    'Latex')
grid on
    
figure(10)
plot(t, df(:,8))
xlabel('Time $[s]$', 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Accelerometer''s velocity $[m/s]$', 'FontSize', 12, ...
    'Interpreter', 'Latex')
grid on

figure(11)
plot(t, df(:,9))
xlabel('Time $[s]$'          , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Output voltage $[V]$', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(12)
plot(t, df(:,10))
xlabel('Time $[s]$'             , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Valve''s position $[m]$', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(13)
plot(t, df(:,11))
xlabel('Time $[s]$'               , 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Valve''s velocity $[m/s]$', 'FontSize', 12, 'Interpreter', 'Latex')
grid on

figure(14)
[T, D] = GOCE_TD(df, t, data);
plot(t, T, t, D)
xlabel('Time $[s]$', 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Spacecraft''s thrust, drag $[N]$', 'FontSize', 12, ...
    'Interpreter', 'Latex')
legend('Thrust', 'Drag', 'FontSize', 12, 'Interpreter', 'Latex', ...
    'Location', 'Best')
grid on

figure(15)
alpha  = 2*pi - 2*acos(1-2*df(:,10)/10/data.A0);
A_valv = data.A0*(alpha-sin(alpha))/2/pi;                                    
plot(t, A_valv)
xlabel('Time $[s]$', 'FontSize', 12, 'Interpreter', 'Latex')
ylabel('Valve''s aperture area $[m^2]$', 'FontSize', 12, 'Interpreter', ...
    'Latex')
grid on