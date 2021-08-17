%% Final project
close all
clearvars
clc
run('GOCE_Data.m');

%% INTEGRATION  
t        = [0, 10*data.T0];                                                 % 10 periods are simulated
cond0    = GOCE_initcond(data);                                             % Initial conditions obtained
options  = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', 10, ...
           'Events', @GOCE_event);
[xx, tt] = GOCE_int(t, data, cond0, options, 0, 0);

GOCE_plots(tt, xx, data)

%% INTEGRATORS COMPARISON 

t        = [0, 10*data.T0];                                                 % 10 periods are simulated
cond0    = GOCE_initcond(data);                                             % Initial conditions obtained

% Eignvalues
eign0 = GOCE_eigenvalues(data.Kpa, data.Kda, data.Kpv, data.Kiv, ...
        data.m_fcv, data.Ki); 
eign0 = sort(eign0);

% Ode23s
fprintf ('for ode23s \n');

options  = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 10, ...        % Set tolerance at 10^-6
           'Events', @GOCE_event);
tic;
sol1 = GOCE_int(t, data, cond0, options, 0, 0, 1);                          % Integration for ode23s 
cputime1 = toc;                                                             % (last input of the integrator function selectes the model)
fprintf ('for tol 10e-6 \n');
fprintf('No. points = %d, \t CPU time = %f \n',size(sol1,1), cputime1);

options  = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', 10, ...        % Set tolerance at 10^-9
           'Events', @GOCE_event);                                          % No more restricted tolerances are tested as the computational time 
tic;                                                                        % for 10^-9 is already high
sol2 = GOCE_int(t, data, cond0, options, 0, 0, 1);                          % Integration for ode23s 
cputime2 = toc; 
fprintf ('for tol 10e-9 \n');
fprintf('No. points = %d, \t CPU time = %f \n',size(sol2,1), cputime2);

% Ode15s
fprintf ('for ode15s \n');

options  = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 10, ...        % Set tolerance at 10^-6
           'Events', @GOCE_event);
tic;
sol1 = GOCE_int(t, data, cond0, options, 0, 0, 2);                          % Integration for ode15s
cputime1 = toc; 
fprintf ('for tol 10e-6 \n');
fprintf('No. points = %d, \t CPU time = %f  \n',size(sol1,1), cputime1);

options  = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', 10, ...        % Set tolerance at 10^-9
           'Events', @GOCE_event);
tic;
sol2 = GOCE_int(t, data, cond0, options, 0, 0, 2);                          % Integration for ode15s
cputime2 = toc; 
fprintf ('for tol 10e-9 \n');
fprintf('No. points = %d, \t CPU time = %f  \n',size(sol2,1), cputime2);

options  = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'MaxStep', 10, ...      % Set tolerance at 10^-12
           'Events', @GOCE_event);
tic;
sol3     = GOCE_int(t, data, cond0, options, 0, 0, 2);                      % Integration for ode15s
cputime3 = toc; 
fprintf ('for tol 10e-12 \n');
fprintf('No. points = %d, \t CPU time = %f  \n',size(sol3,1), cputime3);

% Ode23tb
fprintf ('for ode23tb \n');

options  = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 10, ...        % Set tolerance at 10^-6
           'Events', @GOCE_event);
tic;
sol1     = GOCE_int(t, data, cond0, options, 0, 0, 3);                      % Integration for ode23tb
cputime1 = toc; 
fprintf('for tol 10e-6 \n');
fprintf('No. points = %d, \t CPU time = %f  \n',size(sol1,1), cputime1);

options  = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', 10, ...        % Set tolerance at 10^-9
           'Events', @GOCE_event);                                          % No more restricted tolerances are tested as the computational time 
tic;                                                                        % for 10^-9 is already high
sol2     = GOCE_int(t, data, cond0, options, 0, 0, 3);                      % Integration for ode23tb
cputime2 = toc; 
fprintf ('for tol 10e-9 \n');
fprintf('No. points = %d, \t CPU time = %f\n',size(sol2,1), cputime2);

%% SENSITIVITY ANALYSIS

clearvars -except data
close all
clc

change = (-20:4:20)/100;                                                    % Change (percentage) of the chosen variable
err    = zeros(1, length(change));                                          % err vector is preallocated
lambda = zeros(5,6,length(change));

Kpa0   = data.Kpa;                                                          % Initial variables are saved
Kda0   = data.Kda;
Kpv0   = data.Kpv;
Kiv0   = data.Kiv;
m_fcv0 = data.m_fcv;
Ki0    = data.Ki;
eign_init = GOCE_eigenvalues(data.Kpa, data.Kda, data.Kpv, data.Kiv, data.m_fcv, data.Ki);
eign_init = sort(eign_init);
counter = 0;

for var = 1:6                                                               % Loop for every variable
    for cont = 1:length(change)
        switch var
            case 1
                data.Kpa   = Kpa0  *(1 + change(cont));                     % At each iteration, the variable is incremented in
            case 2                                                          % the percentage indicated by the value of change
                data.Kpa   = Kpa0;                                          % Variables already executed are restored
                data.Kda   = Kda0  *(1 + change(cont));
            case 3
                data.Kda   = Kda0;
                data.Kpv   = Kpv0  *(1 + change(cont));
            case 4
                data.Kpv   = Kpv0;
                data.Kiv   = Kiv0  *(1 + change(cont));
            case 5
                data.Kiv   = Kiv0;
                data.m_fcv = m_fcv0*(1 + change(cont));
            case 6
                data.m_fcv = m_fcv0;
                data.Ki    = Ki0   *(1 + change(cont));
        end
        
        t        = linspace(0, 10*data.T0, 1e4);
        cond0    = GOCE_initcond(data);
        options  = odeset('RelTol',1e-9,'AbsTol',1e-9, 'MaxStep', 10, ...
                   'Events', @GOCE_event);
        [df, tt] = GOCE_int(t, data, cond0, options, 0, 0);
        
        [T, D] = GOCE_TD(df, tt, data);                                     % Thrust and drag for this iteration
        
        err(cont) = norm(T - D);                                            % Error is calculated as the norm of the difference
    end
    
    figure(1)
    plot(change, err)
    hold on
end

figure(1)
grid on
legend('$K_{pa}$', '$K_{da}$', '$K_{pv}$', '$K_{iv}$', '$m_{fcv}$', ...
    '$K_i$', 'Interpreter', 'Latex', 'FontSize', 12, 'Location', 'North')
xlabel('Variation $[\%]$', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('$||T-D||$', 'Interpreter', 'Latex', 'FontSize', 12)

% saveas(figure(1), 'Sensitivity.fig')
% saveas(figure(1), 'Sensitivity.eps', 'epsc')

data.Ki = Ki0;

%% OPTIMIZATION

clearvars -except data
close all
clc

% optimize = 'PD';                                                            % Optimization of Kpa, Kda
% optimize = 'PI';                                                            % Optimization of Kpv, Kiv
% optimize = 'solenoid';                                                      % Optimization of Ki, m_fcv
% optimize = 'all';                                                           % Optimization of Kpa, Kda, Kpv, Kiv, Ki, m_fcv
optimize = 'sensitive';                                                     % Optimization of Kpa, Kiv, Ki

Kpa0   = data.Kpa;                                                          % Initial variables are saved
Kda0   = data.Kda;
Kpv0   = data.Kpv;
Kiv0   = data.Kiv;
m_fcv0 = data.m_fcv;
Ki0    = data.Ki;

% Integration 
t        = linspace(0, 5*data.T0, 1e3);                                     % Five periods
init0    = GOCE_initcond(data);
options  = odeset('RelTol',1e-9,'AbsTol',1e-9, 'events', @GOCE_event);
[df, tt] = GOCE_int(t, data, init0, options, 0, 0);

[T, D] = GOCE_TD(df, tt, data);

norm0 = norm(T - D);                                                        % Value for the initial variables

switch optimize
    case 'PD'
        cond0 = [Kpa0, Kda0];                                               % cond0 is modified for each optimization
    case 'PI'
        cond0 = [Kpv0, Kiv0];
    case 'solenoid'
        cond0 = [m_fcv0, Ki0];
    case 'all'
        cond0 = [Kpa0, Kda0, m_fcv0, Ki0, Kpv0, Kiv0];
    case 'sensitive'
        cond0 = [Kpa0, Ki0, Kiv0];
end

options = optimoptions('fmincon', 'Display', 'iter', 'TolCon', 1e-9, ...
          'TolFun', 1e-9, 'TolX', 1e-9, 'Algorithm', 'active-set');
% options = optimoptions('fmincon', 'Display', 'iter', 'TolCon', 1e-9, ...
%           'TolFun', 1e-9, 'TolX', 1e-9, 'Algorithm', 'interior-point');
lb      = zeros(1, length(cond0));

[sol_opt, fval, exitflag] = fmincon(@(xx) opt(xx, data, optimize, norm0,...
                            t), cond0, [], [], [], [], lb, [], [],options);

switch optimize
    case 'PD'
        data.Kpa = sol_opt(1);                                              % Original variables are modified to the optimized
        data.Kda = sol_opt(2);
    case 'PI'
        data.Kpv = sol_opt(1);
        data.Kiv = sol_opt(2);
    case 'solenoid'
        data.m_fcv = sol_opt(1);
        data.Ki    = sol_opt(2);
    case 'all'
        data.Kpa   = sol_opt(1);
        data.Kda   = sol_opt(2);
        data.m_fcv = sol_opt(3);
        data.Ki    = sol_opt(4);
        data.Kpv   = sol_opt(5);
        data.Kiv   = sol_opt(6);
    case 'sensitive'
        data.Kpa   = sol_opt(1);
        data.Ki    = sol_opt(2);
        data.Kiv   = sol_opt(3);
end

% Integration 
init_opt = GOCE_initcond(data);
options  = odeset('RelTol',1e-9,'AbsTol',1e-9, 'events', @GOCE_event);      % Simulation with optimized variables
[df, tt] = GOCE_int(t, data, init_opt, options, 0, 0);

[T, D] = GOCE_TD(df, tt, data);

norm_opt = norm(T - D);                                                     % norm for the optimized variables

figure(1)
plot(tt, D, tt, T)
grid on
legend('$D$', '$T$', 'Interpreter', 'Latex', 'FontSize', 12, 'Location',...
    'Best')
xlabel('Time $[s]$', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('Spacecraft''s optimized thrust, drag $[N]$', 'Interpreter', ...
    'Latex', 'FontSize', 12)

GOCE_plots(tt,df,data)

Initial   = [cond0   norm0   ]';
Optimized = [sol_opt norm_opt]';
switch optimize
    case 'PD'
        Variables = {'Kpa'; 'Kda'; '||T-D||'};
    case 'PI'
        Variables = {'Kpv'; 'Kiv'; '||T-D||'};
    case 'solenoid'
        Variables = {'m_fcv'; 'Ki'; '||T-D||'};
    case 'all'
        Variables = {'Kpa'; 'Kda'; 'm_fcv'; 'Ki'; 'Kpv'; 'Kiv'; '||T-D||'};
    case 'sensitive'
        Variables = {'Kpa'; 'Ki'; 'Kiv'; '||T-D||'};
end
Comparison = table(Initial, Optimized, 'RowNames', Variables);              % Table comparing initial and optimized variables
disp(Comparison)

data.Kpa   = Kpa0;                                                          % Variables are restored to their initial values
data.Kda   = Kda0;
data.Kpv   = Kpv0;
data.Kiv   = Kiv0;
data.m_fcv = m_fcv0;
data.Ki    = Ki0;

%% FAILURE

clearvars -except data
close all
clc

fail_T    = 0;                                                              % Failure: thruster not working (T = 0)
fail_valv = 1;                                                              % Failure: valve is blocked

t_initial = 0;
t_fail    = 1.8* data.T0;                                                   % Time of failure
t_recover = 1.9* data.T0;                                                   % Time of recovery
t_final   =   4* data.T0;

tv_fail1 = [t_initial       , t_fail   ];                                   % Time before the failure
tv_fail2 = [t_fail    + 1e-9, t_recover];                                   % Duration of the failure
tv_fail3 = [t_recover + 1e-9, t_final  ];                                   % Time after  the failure

cond1      = GOCE_initcond(data);
options    = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', ...
             @GOCE_event, 'MaxStep', 10);
[xx1, tt1] = GOCE_int(tv_fail1, data, cond1, options, 0, 0);
cond2      = xx1(end,:);

if fail_valv == 1
    cond2(11)  = 0;                                                         % Should the valve be blocked, zero initial 
end                                                                         % velocity is imposed

[xx2, tt2] = GOCE_int(tv_fail2, data, cond2, options, fail_T, fail_valv);
cond3      = xx2(end,:);
[xx3, tt3] = GOCE_int(tv_fail3, data, cond3, options, 0, 0);

GOCE_plots([tt1 tt2 tt3]', [xx1; xx2; xx3], data)

[T, D] = GOCE_TD([xx1; xx2; xx3], [tt1 tt2 tt3]', data);

if fail_T == 1
    close 14
    figure(14)
    T(length(xx1)+1:length(xx1) + length(xx2)) = 0;
    plot([tt1 tt2 tt3]', T, [tt1 tt2 tt3]', D)
    xlabel('Time $[s]$', 'FontSize', 12, 'Interpreter', 'Latex')
    ylabel('Spacecraft''s thrust, drag $[N]$', 'FontSize', 12, ...
        'Interpreter', 'Latex')
    legend('Thrust', 'Drag', 'FontSize', 12, 'Interpreter', 'Latex', ...
        'Location', 'Best')
    grid on
end

%% FUNCTIONS

function f = opt(xx, data, optimize, norm0, t)

switch optimize
    case 'PD'
        data.Kpa   = xx(1);
        data.Kda   = xx(2);
    case 'PI'
        data.Kpv   = xx(1);
        data.Kiv   = xx(2);
    case 'solenoid'
        data.m_fcv = xx(1);
        data.Ki    = xx(2);
    case 'all'
        data.Kpa   = xx(1);
        data.Kda   = xx(2);
        data.m_fcv = xx(3);
        data.Ki    = xx(4);
        data.Kpv   = xx(5);
        data.Kiv   = xx(6);
    case 'sensitive'
        data.Kpa   = xx(1);
        data.Ki    = xx(2);
        data.Kiv   = xx(3);
end

cond0    = GOCE_initcond(data);
options  = odeset('RelTol',1e-9,'AbsTol',1e-9, 'events', @GOCE_event);
[df, tt] = GOCE_int(t, data, cond0, options, 0, 0);

[T, D] = GOCE_TD(df, tt, data);

f = norm(T - D)/norm0;                                                      % The norm of the difference divided by the 
                                                                            % initial norm is the function to optimize
end