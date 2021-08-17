%% ORBITAL DATA

data.muE      = 3.98600433*1e5;                                             % Planetary constant of Earth [km^3/s^2]                                                   
data.Re       = 6378.16;                                                    % Earth's radius at Equator [km]
w_earth       = 15.04*pi/180/3600;                                          % Earth's angular velocity [rad/s]
data.ww_earth = [0; 0; w_earth];
data.ra       = data.Re;                                                    % Earth'S (considered ellipsoid) semi-major axis [km]
data.rc       = 6356.778;                                                   % Earth'S (considered ellipsoid) semi-minor axis [km]
data.h0       = 254.9;                                                      % Initial orbit altitude [km]
data.rp       = data.h0+data.ra; %check what Earth radius to add            % Initial orbit perigee radius [km]                             
e0            = 0.0045;                                                     % Initial orbit eccentricity [-]
a0            = data.rp/(1 - e0);                                           % Initial semi-major axis[km]
i0            = 90*pi/180;                                                  % Initial orbit inclination [rad]
O0            = 0;                                                          % Initial orbit RAAN[rad]
o0            = 0;                                                          % Initial orbit argument of perigee[rad]
nu0           = 0;                                                          % Initial true anomaly [rad]
data.kep0     = [a0, e0, i0, O0, o0, nu0];                                  % Keplerian initial elements vector
data.T0       = 2*pi*sqrt(a0^3/data.muE);                                   % Initial orbital period [s]

%% SPACECRAFT DATA

data.M    = 300;                                                            % GOCE Mass [kg]
data.A    = 1.1;                                                            % GOCE Frontal Section [m^2]
data.beta = 300;                                                            % GOCE Ballistic Coefficient
data.cd   = data.M/data.A/data.beta;                                        % Drag coefficient

%% ACCELEROMETER DATA

data.ep    = 8.85e-12;                                                      % Accelerometer Permittivity [F/m]
data.Aa    = 1.6e-3;                                                        % Accelerometer Seismic Mass Section [m^2]
data.m     = 0.32;                                                          % Accelerometer Seismic Mass [kg]
data.g     = 5e-4;                                                          % Electrodes-Mass Gap [m]
data.Vbias = 10;                                                            % Accelerometer Bias Voltage [V]
data.Cf    = 2e-12;                                                         % Capacitance [F]
data.Kpa   = 1e6;                                                           % Controller Proportional Gain (Accelerometer)
data.Kda   = 5e4;                                                           % Controller Derivative Gain (Accelerometer)

%% VALVE DATA

data.A0    = 1e-5;                                                          % Valve Orifice Area [m^2]
data.m_fcv = 2e-1;                                                          % Valve Spool Mass [kg]
data.Ki    = 0.2;                                                           % Proportionality Coefficient Current-Spool
data.Kpv   = 0.1;                                                           % Controller Proportional Gain (Valve)
data.Kiv   = 3;                                                             % Controller Integral Gain (Valve)
data.Kfcv  = 7e3;                                                           % Valve Spring Coefficient 
data.c     = 30;                                                            % Valve Friction Coefficient
data.cdis  = 0.61;                                                          % Discharge coefficient of the Valve
Ru         = 8.3145;                                                        % Gas constant [J/K/mol]
m_xen      = 131.293*1e-3;                                                  % Molar mass of Xenon [kg/mol]
data.R     = Ru/m_xen;

%% THRUSTER DATA

data.T2    = 240;                                                           % Xenon Working Temperature
data.p2    = 2e5;                                                           % Xenon Working Pressure
data.k     = 1.66;                                                          % Xenon Specific Heat Ratio
data.m_i   = 2.188e-25;                                                     % Xenon Ion Mass
data.e     = 1.6e-19;                                                       % Electron Charge
data.DV    = 2e3;                                                           % Acceleration Grid Voltage