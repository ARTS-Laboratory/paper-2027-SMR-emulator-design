%====================================================================
% steadySolve2.m
%  Steady‐state driver for the SG‐reactor model.
%  Uses fsolve to find 13 unknowns X = [T_HL; T_CL; T_f; L1; L2; P_s;
%               h_out; T_m3; T_m2; T_m1; T_p5; T_p3; T_p1]
%  for a given fixed m_dot_s.  Then prints all unknowns plus 
%  the derived quantities: m_s, T_out, T_p2, T_p4, T_sh, T_sc, L3.
%====================================================================

clear
clc

%% 1) Build all constant parameters
params = buildParams4();

%% 2) Specify the steady‐state SG mass flow (kg/s)
%m_dot_s_fixed = 67.07;  % e.g. from plant design or operator setpoint

%% 3) Provide an initial guess X0 for the 13 unknowns
%   Order must match reactorSteadyResidual:
%     X(1)=T_HL, X(2)=T_CL, X(3)=T_f, 
%     X(4)=L1,   X(5)=L2,   X(6)=P_s, 
%     X(7)=h_out, X(8)=T_m3, X(9)=T_m2, 
%     X(10)=T_m1, X(11)=T_p5, X(12)=T_p3, 
%     X(13)=T_p1.
X0 = [ ...
    310;      % T_HL guess [°C]
    258;      % T_CL guess [°C]
    520;      % T_f  guess [°C]
    2.5;        % L1   guess [m]
    18.5;        % L2   guess [m]
    34.48;    % P_s  guess [bar]
    2.998e6;  % h_out guess [J/kg]
    256;      % T_m3 guess [°C]
    282;      % T_m2 guess [°C]
    304;      % T_m1 guess [°C]
    263;      % T_p5 guess [°C]
    287;      % T_p3 guess [°C]
    308       % T_p1 guess [°C]
];

%% 4) Call fsolve to find the steady‐state X_ss
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-9, ...
    'StepTolerance',1e-9);

[X_ss, fval, exitflag, output] = fsolve( ...
    @(X) reactorSteadyResidual2(X, params), ...
    X0, opts);

if exitflag <= 0
    error('Steady‐state fsolve did not converge.');
end

%% 5) Unpack the steady‐state solution X_ss
T_HL_ss  = X_ss(1);    % [°C]
T_c1_ss  = X_ss(2);    % [°C]
T_f_ss   = X_ss(3);    % [°C]
L1_ss    = X_ss(4);    % [m]
L2_ss    = X_ss(5);    % [m]
P_s_ss   = X_ss(6);    % [bar]
h_out_ss = X_ss(7);    % [J/kg]
T_m3_ss  = X_ss(8);    % [°C]
T_m2_ss  = X_ss(9);    % [°C]
T_m1_ss  = X_ss(10);   % [°C]
T_p5_ss  = X_ss(11);   % [°C]
T_p3_ss  = X_ss(12);   % [°C]
T_p1_ss  = X_ss(13);   % [°C]

%% 6) Compute the extra “on‐the‐fly” / dependent variables

% 6.1 SG secondary mass flow (fixed)
m_s_ss = 67.07;  % [kg/s]

% 6.2 L3 = L_T - (L1 + L2)
L3_ss = params.L_T - (L1_ss + L2_ss);  % [m]

% 6.3 Primary‐side intermediate temperatures
T_p2_ss = 2*T_p1_ss - T_HL_ss;  % [°C]
T_p4_ss = 2*T_p3_ss - T_p2_ss;  % [°C]
T_CL_ss = 2*T_p5_ss - T_p4_ss;

% 6.4 Secondary‐side properties from XSteam
T_sat_ss = XSteam('Tsat_p', P_s_ss);          % saturation T at P_s [°C]
h_in_ss  = XSteam('h_pT', P_s_ss, 148.0045) * 1000;  % feed enthalpy [J/kg]
h_L_ss   = XSteam('hL_T', T_sat_ss) * 1000;  % [J/kg]
h_V_ss   = XSteam('hV_T', T_sat_ss) * 1000;  % [J/kg]

% 6.5 Two‐phase mixture enthalpies
h_sh_ss = 0.5*(h_out_ss + h_V_ss);  % [J/kg]
h_sc_ss = 0.5*(h_in_ss  + h_L_ss);  % [J/kg]

% 6.6 Two‐phase mixture temperatures
T_sh_ss = XSteam('T_ph', P_s_ss, (h_sh_ss/1000));  % [°C]
T_sc_ss = XSteam('T_ph', P_s_ss, (h_sc_ss/1000));  % [°C]

% 6.7 “Outlet” fluid temperature
T_out_ss = XSteam('T_ph', P_s_ss, (h_out_ss/1000)); % [°C]

% --- values that were “fixed” in the steady solver
phi_ss  = 1.0000;           % dimensionless
C_ss    = 3500;             % J/(kg·K)

%% 7) Evaluate all 13 residuals at the final solution
F_final = reactorSteadyResidual2(X_ss, params);  % a 13×1 vector

%% 8) Print all results

fprintf('\n=== Steady‐State Solution (fsolve) ===\n\n');

% 8.1 Print the “unknowns vector”
fprintf(' Unknowns vector (X_ss):\n');
fprintf('   T_HL  = %.4f °C\n', T_HL_ss);
fprintf('   T_CL  = %.4f °C\n', T_CL_ss);
fprintf('   T_c1  = %.4f °C\n', T_c1_ss);
fprintf('   T_f   = %.4f °C\n', T_f_ss);
fprintf('   L1    = %.4f m\n',    L1_ss);
fprintf('   L2    = %.4f m\n',    L2_ss);
fprintf('   P_s   = %.4f bar\n',  P_s_ss);
fprintf('   h_out = %.4f J/kg\n',h_out_ss);
fprintf('   T_m3  = %.4f °C\n',  T_m3_ss);
fprintf('   T_m2  = %.4f °C\n',  T_m2_ss);
fprintf('   T_m1  = %.4f °C\n',  T_m1_ss);
fprintf('   T_p5  = %.4f °C\n',  T_p5_ss);
fprintf('   T_p3  = %.4f °C\n',  T_p3_ss);
fprintf('   T_p1  = %.4f °C\n',  T_p1_ss);
fprintf('\n');

% 8.2 Print the “on‐the‐fly” / dependent variables

fprintf(' Derived variables:\n');
fprintf('   m_s (secondary mass flow) = %.4f kg/s\n', m_s_ss);
fprintf('   L3 (third tube segment)   = %.4f m\n',    L3_ss);
fprintf('   T_p2 (node‐2 fluid T)     = %.4f °C\n',  T_p2_ss);
fprintf('   T_p4 (node‐4 fluid T)     = %.4f °C\n',  T_p4_ss);
fprintf('   T_sh  (two‐phase mix T)   = %.4f °C\n',  T_sh_ss);
fprintf('   T_sc  (two‐phase lqd T)   = %.4f °C\n',  T_sc_ss);
fprintf('   T_out (SG outlet fluid T) = %.4f °C\n',  T_out_ss);

fprintf('\n');
fprintf(' Secondary‐side enthalpies & temperatures at P_s:\n');
fprintf('   T_sat = %.4f °C   (sat. temp at P_s)\n', T_sat_ss);
fprintf('   h_in  = %.4f J/kg (feed enthalpy at T_in)\n', h_in_ss);
fprintf('   h_L   = %.4f J/kg (liquid enthalpy at sat)\n', h_L_ss);
fprintf('   h_V   = %.4f J/kg (vapor enthalpy at sat)\n', h_V_ss);
fprintf('   h_sh  = %.4f J/kg (two‐phase ave)\n',       h_sh_ss);
fprintf('   h_sc  = %.4f J/kg (two‐phase ave)\n',       h_sc_ss);

fprintf('\nSteady‐state solve completed.\n');
fprintf('fsolve exitflag = %d\n', exitflag);

% 8.3 Print the 13 residuals at the solution (should all be ≈0)
fprintf(' Residual vector F_final (should be ≈0):\n');
for k = 1:13
    fprintf('   F(%2d) = %+12.6e\n', k, F_final(k));
end
fprintf('\n');

fprintf(' fsolve exitflag = %d\n', exitflag);
fprintf(' fsolve output:\n');
disp(output);
fprintf('\n');

% Define the same index map your dynamic model uses
  idx.phi   = 1;  
  idx.C     = 2;
  idx.T_f   = 3;
  idx.T_c1  = 4;
  idx.T_HL  = 5;
  idx.T_CL  = 6;
  idx.L1    = 7;
  idx.L2    = 8;
  idx.P_s   = 9;
  idx.h_out = 10;
  idx.T_m3  = 11;
  idx.T_m2  = 12; 
  idx.T_m1  = 13;
  idx.T_p5  = 14; 
  idx.T_p3  = 15;
  idx.T_p1  = 16;  
  

% --- allocate and fill
X_dyn0           = zeros(16,1);
X_dyn0(idx.phi)  = phi_ss;
X_dyn0(idx.C)    = C_ss;
X_dyn0(idx.T_f)  = T_f_ss;
X_dyn0(idx.T_c1) = T_c1_ss;
X_dyn0(idx.T_HL) = T_HL_ss;
X_dyn0(idx.T_CL) = T_CL_ss;
X_dyn0(idx.P_s)  = P_s_ss;
X_dyn0(idx.h_out)= h_out_ss;
X_dyn0(idx.T_p1) = T_p1_ss;
X_dyn0(idx.T_p3) = T_p3_ss;
X_dyn0(idx.T_p5) = T_p5_ss;
X_dyn0(idx.T_m1) = T_m1_ss;
X_dyn0(idx.T_m2) = T_m2_ss;
X_dyn0(idx.T_m3) = T_m3_ss;
X_dyn0(idx.L1)   = L1_ss;
X_dyn0(idx.L2)   = L2_ss;

% (optional) quick sanity check
assert(all(~isnan(X_dyn0)),'Some entries are NaN')

format longG         % or  format longE  (12-digit scientific)
disp(X_dyn0)         % now shows full precision
%% 8) (Optional) Save results to a .mat file for later dynamic use
save('steadyStateResults.mat', ...
     'X_ss', 'X_dyn0', ...          % <-- new line!
     'm_s_ss', 'T_out_ss', 'T_p2_ss', 'T_p4_ss', ...
     'T_sh_ss', 'T_sc_ss', 'L3_ss', 'params');

%====================================================================
% end of steadySolve.m
%====================================================================
