function F = reactorSteadyResidual2(X, params)
% reactorSteadyResidual  Steady‐state residuals for SG‐reactor model
%
%   F = reactorSteadyResidual(X, params) returns a 13×1 vector of steady‐state
%   equations, given the 13 unknowns in X and the constant parameters in params.
%
%   Unknown vector X (all units in SI, as described below):
%     X(1)  = T_HL    = Hot‐leg coolant temperature [°C]
%     X(2)  = T_CL    = Cold‐leg coolant temperature [°C]
%     X(3)  = T_f     = Fuel temperature [°C]
%     X(4)  = L1      = Tube segment‐1 length [m]
%     X(5)  = L2      = Tube segment‐2 length [m]
%     X(6)  = P_s     = SG secondary pressure [bar]
%     X(7)  = h_out   = Outlet enthalpy from SG [J/kg]
%     X(8)  = T_m3    = Metal‐wall (tube wall) temperature at segment 3 [°C]
%     X(9)  = T_m2    = Metal‐wall temperature at segment 2 [°C]
%    X(10)  = T_m1    = Metal‐wall temperature at segment 1 [°C]
%    X(11)  = T_p5    = Primary‐side fluid temperature at node 5 [°C]
%    X(12)  = T_p3    = Primary‐side fluid temperature at node 3 [°C]
%    X(13)  = T_p1    = Primary‐side fluid temperature at node 1 [°C]
%
%   The constant‐parameter struct "params" must include fields:
%     params.n, params.pi_val, params.D_o, params.D_i
%     params.alpha_o, params.alpha_i, params.m_dot_p_rated
%     params.P_primary, params.h_fc, params.A_fc
%     params.tau, params.P_0, params.eta_SG, params.L_T
%     params.gamma_bar (not used here, but may be needed elsewhere)
%
%   We also fix, during steady‐state, these three variables:
%     phi   = 1            (dimensionless)
%     C     = 3500         (J/(kg·K), unused in these equations but needed by dynamic model)
%     T_c1  = 284.05556    (°C)
%
%   The remaining intermediate variables computed inside:
%     T_sat   = saturation temperature at P_s  [°C]
%     L3      = params.L_T − (L1 + L2)         [m]
%     T_p2    = 2·T_p1 − T_HL                   [°C]
%     T_p4    = 2·T_p3 − T_p2                   [°C]
%     T_out   = XSteam('T_ph', P_s, h_out/1000) [°C]
%     h_in    = XSteam('h_pT', P_s, T_in)·1000  [J/kg]
%     h_L     = XSteam('hL_T', T_sat)·1000      [J/kg]
%     h_V     = XSteam('hV_T', T_sat)·1000      [J/kg]
%     h_sh    = (h_out + h_V)/2                 [J/kg]
%     h_sc    = (h_in  + h_L)/2                 [J/kg]
%     T_sh    = XSteam('T_ph', P_s, h_sh/1000)  [°C]
%     T_sc    = XSteam('T_ph', P_s, h_sc/1000)  [°C]
%     m_dot_s = 1.9452 · P_s                    [kg/s]
%
%   The 13 residuals F(1) … F(13) correspond to:
%     F(1–3)   = Primary‐loop heat‐transfer balances in SG segments 1, 2, 3
%     F(4–6)   = Wall heat‐transfer balances in SG segments 1, 2, 3
%     F(7–9)   = Secondary‐side energy balances in SG
%     F(10–11) = Coolant temperature balances (fuel‐clad + primary loop)
%     F(12)    = Governor‐valve (mass‐flow) equation
%     F(13)    = Fuel energy balance
%
%   At steady state, all F(i) must equal zero.  fsolve will move the
%   13 unknowns in X until these equations are all satisfied.
%
%   Units summary (all SI, except pressures in bar as required by XSteam):
%     Temperatures [°C], lengths [m], enthalpies [J/kg], pressures [bar],
%     mass flows [kg/s], heat‐transfer coefficients [W/(m²·K)], etc.

  %---------------------------------------------------------------------
  % 1) Unpack unknown vector X
  %---------------------------------------------------------------------
  T_HL  = X(1);    % Hot‐leg coolant temp [°C]
  T_c1  = X(2);    % Cold‐leg coolant temp [°C]
  T_f   = X(3);    % Fuel temperature [°C]
  L1    = X(4);    % Tube length segment 1 [m]
  L2    = X(5);    % Tube length segment 2 [m]
  P_s   = X(6);    % SG secondary pressure [bar]
  h_out = X(7);    % Outlet enthalpy [J/kg]
  T_m3  = X(8);    % Metal‐wall temp seg 3 [°C]
  T_m2  = X(9);    % Metal‐wall temp seg 2 [°C]
  T_m1  = X(10);   % Metal‐wall temp seg 1 [°C]
  T_p5  = X(11);   % Primary fluid temp node‐5 [°C]
  T_p3  = X(12);   % Primary fluid temp node‐3 [°C]
  T_p1  = X(13);   % Primary fluid temp node‐1 [°C]

  %---------------------------------------------------------------------
  % 2) Unpack constant parameters from params
  %---------------------------------------------------------------------
  n              = params.n;              % number of tubes [–]
  pi_val         = params.pi_val;         % π [–]
  D_o            = params.D_o;            % tube outer diameter [m]
  D_i            = params.D_i;            % tube inner diameter [m]
  alpha_o        = params.alpha_o;        % outside HTC [W/(m²·K)]
  alpha_i        = params.alpha_i;        % inside HTC [W/(m²·K)]
  m_dot_p_rated  = params.m_dot_p_rated;  % primary mass flow [kg/s]
  P_primary      = params.P_primary;      % primary loop pressure [bar]
  h_fc           = params.h_fc;           % fuel‐clad HTC [W/(m²·K)]
  A_fc           = params.A_fc;           % fuel‐clad area [m²]
  tau            = params.tau;            % fraction of power into fuel [–]
  P_0            = params.P_0;            % reactor thermal power [W]
  eta_SG         = params.eta_SG;         % SG efficiency [–]
  L_T            = params.L_T;            % total tube length [m]
  % (params.gamma_bar is not used here)

  %---------------------------------------------------------------------
  % 3) Assign “fixed” variables in steady‐state
  %---------------------------------------------------------------------
  phi  = 1.0000;        % dimensionless (held fixed)
  C    = 3500;          % J/(kg·K) (held fixed, not used in these eqns)
  %T_c1 = 284.05556;     % coolant reference temp [°C] (held fixed)
  T_in = 148.0045;      % feed inlet temp [°C]  (held fixed)
  m_dot_s = 67.07;      % mass flow rate [kg/s] (held fixed)


  %---------------------------------------------------------------------
  % 4) Compute intermediate (secondary) variables
  %---------------------------------------------------------------------
  % 4.1 Saturation temperature at P_s [°C]
  T_sat = XSteam('Tsat_p', P_s);

  % 4.2 Tube‐length segment #3 [m]
  L3 = L_T - (L1 + L2);

  % 4.3 Primary‐side interpolation nodes
  T_p2 = 2*T_p1 - T_HL;   % node‐2 temp [°C]
  T_p4 = 2*T_p3 - T_p2;   % node‐4 temp [°C]
  %T_p5 = (T_p4 + T_CL)/2;   % node‐5 temp [°C]
  T_CL = 2*T_p5 - T_p4;   % Cold Leg temp [°C]


  % 4.4 Secondary enthalpies [J/kg]
  T_out = XSteam('T_ph', P_s, h_out/1000);                % [°C]
  h_in  = XSteam('h_pT', P_s, T_in) * 1000;               % [J/kg]
  h_L   = XSteam('hL_T', T_sat) * 1000;                   % [J/kg]
  h_V   = XSteam('hV_T', T_sat) * 1000;                   % [J/kg]

  % 4.5 Mixed enthalpies & fluid temps
  h_sh = (h_out + h_V)/2;    % [J/kg]
  h_sc = (h_in  + h_L)/2;    % [J/kg]
  T_sh = XSteam('T_ph', P_s, h_sh/1000);  % [°C]
  T_sc = XSteam('T_ph', P_s, h_sc/1000);  % [°C]

  %---------------------------------------------------------------------
  % 5) Build the 13 steady‐state residual equations F(1) … F(13)
  %---------------------------------------------------------------------

  %-- 5.1 Primary‐loop heat‐transfer balances in SG segments
  %   F(1): Segment #1 (length L3, between T_m1 and T_p1)
  F(1) = n * pi_val * D_o * alpha_o * L3 * (T_m1 - T_p1) ...
         + m_dot_p_rated * (phi^(1/3)) * ( ...
           (XSteam('h_pT', P_primary, T_HL) * 1000) - ...
           (XSteam('h_pT', P_primary, T_p2) * 1000) );

  %   F(2): Segment #2 (length L2, between T_m2 and T_p3)
  F(2) = n * pi_val * D_o * alpha_o * L2 * (T_m2 - T_p3) ...
         + m_dot_p_rated * (phi^(1/3)) * ( ...
           (XSteam('h_pT', P_primary, T_p2) * 1000) - ...
           (XSteam('h_pT', P_primary, T_p4) * 1000) );

  %   F(3): Segment #3 (length L1, between T_m3 and T_p5)
  F(3) = n * pi_val * D_o * alpha_o * L1 * (T_m3 - T_p5) ...
         + m_dot_p_rated * (phi^(1/3)) * ( ...
           (XSteam('h_pT', P_primary, T_p4) * 1000) - ...
           (XSteam('h_pT', P_primary, T_CL)  * 1000) );

  %-- 5.2 Wall heat‐transfer balances (tube‐metal) in SG segments
  %   F(4): Segment #1 (length L3), wall between T_p1 and T_m1
  F(4) = eta_SG * n * pi_val * D_o * alpha_o * L3 * (T_p1 - T_m1) ...
         + n * pi_val * D_i * alpha_i * L3 * (T_sh - T_m1);

  %   F(5): Segment #2 (length L2), wall between T_p3 and T_m2
  F(5) = eta_SG * n * pi_val * D_o * alpha_o * L2 * (T_p3 - T_m2) ...
         + n * pi_val * D_i * alpha_i * L2 * (T_sat - T_m2);

  %   F(6): Segment #3 (length L1), wall between T_p5 and T_m3
  F(6) = eta_SG * n * pi_val * D_o * alpha_o * L1 * (T_p5 - T_m3) ...
         + n * pi_val * D_i * alpha_i * L1 * (T_sc - T_m3);

  %-- 5.3 Secondary‐side energy balances in SG
  %   F(7): Between two‐phase (h_V) and outlet (h_out), length L3
  F(7) = m_dot_s * (h_V - h_out) ...
         + n * pi_val * D_i * alpha_i * L3 * (T_m1 - T_sh);

  %   F(8): Between liquid (h_L) and vapor (h_V), length L2
  F(8) = m_dot_s * (h_L - h_V) ...
         + n * pi_val * D_i * alpha_i * L2 * (T_m2 - T_sat);

  %   F(9): Between feed (h_in) and liquid (h_L), length L1
  F(9) = m_dot_s * (h_in - h_L) ...
         + n * pi_val * D_i * alpha_i * L1 * (T_m3 - T_sc);

  %-- 5.4 Coolant‐temperature equations (primary loop + fuel‐clad)
  %   F(10): Energy into cold leg: (1−τ)P₀ + h_fc A_fc (T_f − T_c1)
  %           plus primary‐loop enthalpy change from CL→HL
  F(10) = ( (1 - tau)*P_0 * phi ) ...
          + h_fc * A_fc * (T_f - T_c1) ...
          + m_dot_p_rated * (phi^(1/3)) * ( ...
            (XSteam('h_pT', P_primary, T_CL) * 1000) - ...
            (XSteam('h_pT', P_primary, T_HL) * 1000) );

  %   F(11): Geometric constraint for coolant node T_c1:
  %          T_c1 = (T_HL + T_CL)/2  at steady state
  F(11) = T_c1 - ( (T_HL + T_CL)/2 );

  %-- 5.5 Governor valve (mass‐flow) equation
  F(12) = m_dot_s - ( 1.9452 * P_s );

  %-- 5.6 Fuel energy balance (fuel‐clad side)
  F(13) = tau * P_0 * phi + h_fc * A_fc * (T_c1 - T_f);

  %---------------------------------------------------------------------
  % 6) Return column vector
  %---------------------------------------------------------------------
  F = F(:);
end
