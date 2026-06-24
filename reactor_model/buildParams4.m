function params = buildParams4()
% buildParams Returns a struct of all constant parameters for the reactor DAE model.
% Units noted in comments.

% — Neutron kinetics & reactivity —
params.alpha_f       = -2.16e-5;       % Fuel temperature reactivity coeff. [1/°C]
params.alpha_c       = -1.8e-4;        % Moderator temperature reactivity coeff. [1/°C]
params.Lambda        = 2e-5;           % Prompt neutron generation time [s]
params.beta          = 0.007;          % Delayed neutron fraction [–]
params.lambda_val    = 0.1;            % Precursor decay constant [1/s]
params.T_ref_f       = 518.6012;            % Reference fuel temperature [°C]
params.T_ref_c       = 284.05556;            % Reference coolant temperature [°C]

% — Reactor power & efficiency —
params.P_0           = 160e6;          % Reactor thermal power [W]
params.tau           = 0.97;           % Fraction of power deposited in fuel [–]
params.eta_T         = 0.69243719;     % Turbine efficiency [–]
params.eta_SG        = 0.9946;         % Steam generator efficiency [–]

% — Mass flow & inlet —
params.m_dot_p_rated = 587;            % Rated primary mass flow [kg/s]

% — Tube & SG geometry —
params.n             = 1380;           % Number of SG tubes [–]
params.pi_val        = pi;             % Pi constant
params.D_o           = 0.0159;         % Tube outer diameter [m]
params.D_i           = 0.01336;        % Tube inner diameter [m]
params.L_T           = 24.2;           % Total tube length [m]
params.A_w           = (params.pi_val/4)*(params.D_o.^2 - params.D_i.^2); % Tube wall area [m^2]
params.alpha_o       = 19093.9714;     % Outside heat transfer coeff. [W/m^2·K]
params.alpha_i       = 2697.1930;      % Inside heat transfer coeff. [W/m^2·K]

% — Heat transfer surfaces —
params.h_fc          = 1135;           % Fuel-clad heat transfer coeff. [W/m^2·K]
params.A_fc          = 583;            % Fuel-clad heat transfer area [m^2]

% — Flow areas —
params.A_p           = 0.7266;         % Primary flow area [m^2]
params.A_s           = 0.1935;         % Secondary flow area [m^2]

% — System volumes —
params.V_core        = 2.5202;         % Core coolant volume [m^3]
params.V_HL          = 17.9812;        % Hot leg volume [m^3]
params.V_CL          = 16.3671;        % Cold leg volume [m^3]
params.V_SG_primary  = 17.5848;        % SG primary side volume [m^3]

% — Primary loop pressure —
params.P_primary     = 127.6;          % Primary loop pressure [bar]

% — Fuel properties —
params.m_f           = 11252;          % Fuel mass [kg]
params.c_pf          = 0.467e3;        % Fuel specific heat [J/kg·K]

% — Tube-metal properties —
params.rho_w         = 8192;           % Tube wall density [kg/m^3]
params.c_w           = 463;            % Tube wall specific heat [J/kg·K]

% — Saturated mixture void fraction —
params.gamma_bar     = 0.251913;       % Assumed constant void fraction [–]

%params.C0 = (beta / (Lambda * lambda)) * 1;

% — Mechanical power output of Turbine Reference —
%params.P_mech_ref = 45e6;  % whatever your design point is
end
