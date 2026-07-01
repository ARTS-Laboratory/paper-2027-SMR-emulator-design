function dy = reactor_ode4(t, y, params, m_dot_out1, m_dot_in1)

%% 1. Unpack parameters for readability
alpha_f       = params.alpha_f;
alpha_c       = params.alpha_c;
Lambda        = params.Lambda;
beta          = params.beta;
lambda_val    = params.lambda_val;
T_ref_f       = params.T_ref_f;
T_ref_c       = params.T_ref_c;
P_0           = params.P_0;
tau           = params.tau;
eta_T         = params.eta_T;
eta_SG        = params.eta_SG;
m_dot_p_rated = params.m_dot_p_rated;
n             = params.n;
pi_val        = params.pi_val;
D_o           = params.D_o;
D_i           = params.D_i;
L_T           = params.L_T;
A_w           = params.A_w;
alpha_o       = params.alpha_o;
alpha_i       = params.alpha_i;
h_fc          = params.h_fc;
A_fc          = params.A_fc;
A_p           = params.A_p;
A_s           = params.A_s;
V_core        = params.V_core;
V_HL          = params.V_HL;
V_CL          = params.V_CL;
V_SG_primary  = params.V_SG_primary;
P_primary     = params.P_primary;
m_f           = params.m_f;
c_pf          = params.c_pf;
rho_w         = params.rho_w;
c_w           = params.c_w;
gamma_bar     = params.gamma_bar;

% Input variables
T_in         = params.T_in_fun(t);
rho_ext      = params.rho_ext_fun(t);
Valve        = params.Valve_fun(t);


  % Index map ----
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
  
  
   

  % unpack states & params ---
  phi     = y(idx.phi);    C      = y(idx.C);
  T_f     = y(idx.T_f);    T_c1   = y(idx.T_c1);
  T_HL    = y(idx.T_HL);   T_CL   = y(idx.T_CL);
  P_s     = y(idx.P_s);    h_out  = y(idx.h_out);
  T_p1    = y(idx.T_p1);   T_p3   = y(idx.T_p3);
  T_p5    = y(idx.T_p5);   T_m1   = y(idx.T_m1);
  T_m2    = y(idx.T_m2);   T_m3   = y(idx.T_m3);
  L1      = y(idx.L1);     L2     = y(idx.L2);



  %% Derived parameters (must update each call) ---

  % moving‐boundary constraint
  L3      = L_T - (L1 + L2);

  % reactivity
  react   = rho_ext + alpha_f*(T_f-T_ref_f) + alpha_c*(T_c1-T_ref_c);

  % valve & mass‐flow
  %m_dot_out = Valve*P_s;       % governor law
  m_dot_out = m_dot_out1;
  %m_dot_in  = m_dot_out;       % feedwater ⇄ steam balance
  m_dot_in = m_dot_in1;


  % primary flow
  %m_p     = m_dot_p_rated * phi^(1/3);
  m_p     = m_dot_p_rated;


  % saturated‐steam temperature
  T_sat   = XSteam('Tsat_p', P_s);

  %Secondary enthalpies [J/kg]
  T_out = XSteam('T_ph', P_s, h_out/1000);                % [°C]
  h_in  = XSteam('h_pT', P_s, T_in) * 1000;               % [J/kg]
  h_L   = XSteam('hL_p',P_s)*1e3;                        % [J/kg]
  h_V   = XSteam('hV_p',P_s)*1e3;                        % [J/kg]
  h_sat = gamma_bar*h_V   + (1-gamma_bar)*h_L;

  %subcooled and superheated enthalpies & fluid temps
  h_sh = (h_out + h_V)/2;    % [J/kg]
  h_sc = (h_in  + h_L)/2;    % [J/kg]
  T_sh = XSteam('T_ph', P_s, h_sh/1000);  % [°C]
  T_sc = XSteam('T_ph', P_s, h_sc/1000);  % [°C]


  % averages / moving‐boundary temperatures
  T_c2    = 2*T_c1 - T_CL;
  T_p2    = 2*T_p1 - T_HL;
  T_p4    = 2*T_p3 - T_p2;
  T_p6    = 2*T_p5 - T_p4;

%% Precompute thermophysical properties via XSteam

%Densities in secondary side
rho_V  = XSteam('rhoV_p',P_s);      % kg m-3
rho_L  = XSteam('rhoL_p',P_s);
rho_sat     = gamma_bar*rho_V + (1-gamma_bar)*rho_L;
rho_sh = XSteam('rho_ph',P_s, h_sh/1e3);
rho_sc = XSteam('rho_ph',P_s, h_sc/1e3);

%Densities in primary side
rho_Core    = XSteam('rho_pT',   P_primary, T_c1);
rho_HL      = XSteam('rho_pT',   P_primary, T_HL);
rho_p1      = XSteam('rho_pT',   P_primary, T_p1);
rho_p3      = XSteam('rho_pT',   P_primary, T_p3);
rho_p5      = XSteam('rho_pT',   P_primary, T_p5);
rho_CL      = XSteam('rho_pT',   P_primary, T_CL);

%primary enthalpies [J/kg]
h_c2        = XSteam('h_pT',     P_primary, T_c2)*1000;
h_HL        = XSteam('h_pT',     P_primary, T_HL)*1000;
h_CL        = XSteam('h_pT',     P_primary, T_CL)*1000;
h_p2        = XSteam('h_pT',     P_primary, T_p2)*1000;
h_p4        = XSteam('h_pT',     P_primary, T_p4)*1000;
h_p6        = XSteam('h_pT',     P_primary, T_p6)*1000;

%primary Specific isobaric heat capacity [J/kg°C]
Cp_core     = XSteam('Cp_pT',    P_primary, T_c1)*1000;
Cp_p1       = XSteam('Cp_pT',    P_primary, T_p1)*1000;
Cp_p3       = XSteam('Cp_pT',    P_primary, T_p3)*1000;
Cp_p5       = XSteam('Cp_pT',    P_primary, T_p5)*1000;


%entropy and enthalpy of turbine for calculation of Pmech
S_Turbine   = XSteam('s_ph',     P_s, (h_out/1000));
h_isentropic = XSteam('h_ps', 0.08, S_Turbine)*1000;
   
% core coolant density & mass
m_c     = V_core * rho_Core;                  % [kg]

% hot‐leg & cold‐leg masses
m_HL    = V_HL  * rho_HL;                    % [kg]
m_CL    = V_CL  * rho_CL;                    % [kg]

% corresponding transit times
tau_HL  = m_HL  / m_p;                       % [s]
tau_CL  = m_CL  / m_p;                       % [s]

% turbine mechanical power
P_mech  = eta_T * (h_out - h_isentropic) * m_dot_out  ;


%====================================================================
% XSteam-based property slopes            (Pressure in bar, energy in J·kg⁻¹)
%
%  ρ  : density         [kg·m⁻³]
%  h  : specific enthalpy[J·kg⁻¹]
%  V  : saturated vapour
%  L  : saturated liquid
%====================================================================
% Required inputs already in workspace:
%   P_s        – current secondary pressure           [bar]
%   h_out      – outlet (superheated) enthalpy        [J·kg⁻¹]
%   h_in       – feed-water (sub-cooled) enthalpy     [J·kg⁻¹]
%   gamma_bar  – void fraction of the two-phase zone  [–]
%   T_sh, T_sc – mean super-heated / sub-cooled temps [°C]
%
fdiff = @(f,x,dx) (f(x+dx) - f(x-dx)) ./ (2*dx);   % central-difference helper

%--------------------------------------------------------------------
% 0) Adaptive finite-difference steps
%     (keep relative step ≈0.1 % but never smaller than an absolute floor)
%--------------------------------------------------------------------
dP = max(0.001*P_s  , 0.005);      % pressure step   [bar]
dh = max(0.001*h_out,    50);      % enthalpy step   [J·kg⁻¹]
dT = max(0.001*T_sh ,  0.01);      % temperature step[°C] (unused here)

%--------------------------------------------------------------------
% 1) Pure-phase saturation properties  — functions of P only
%--------------------------------------------------------------------
h_V   = XSteam('hV_p'  , P_s) * 1e3;   % sat-vapour enthalpy          [J·kg⁻¹]
h_L   = XSteam('hL_p'  , P_s) * 1e3;   % sat-liquid enthalpy          [J·kg⁻¹]
rho_V = XSteam('rhoV_p', P_s);         % sat-vapour density           [kg·m⁻³]
rho_L = XSteam('rhoL_p', P_s);         % sat-liquid density           [kg·m⁻³]

% ∂h/∂P at saturation (vapour & liquid)          [J·kg⁻¹·bar⁻¹]
dhV_dP    = fdiff(@(Pb) XSteam('hV_p',Pb)*1e3 , P_s, dP);
dhL_dP    = fdiff(@(Pb) XSteam('hL_p',Pb)*1e3 , P_s, dP);

% ∂ρ/∂P at saturation                           [kg·m⁻³·bar⁻¹]
d_rhoV_dP = fdiff(@(Pb) XSteam('rhoV_p',Pb)   , P_s, dP);
d_rhoL_dP = fdiff(@(Pb) XSteam('rhoL_p',Pb)   , P_s, dP);

% ∂(ρh)/∂P  = ρ·(∂h/∂P)+h·(∂ρ/∂P)             [J·m⁻³·bar⁻¹]
d_rhoh_V_dP = rho_V .* dhV_dP + h_V .* d_rhoV_dP;
d_rhoh_L_dP = rho_L .* dhL_dP + h_L .* d_rhoL_dP;

%--------------------------------------------------------------------
% 2) Mean values in super-heated (sh) & sub-cooled (sc) zones
%--------------------------------------------------------------------
h_sh = 0.5*(h_out + h_V);     % mean super-heated enthalpy  [J·kg⁻¹]
h_sc = 0.5*(h_in  + h_L);     % mean sub-cooled  enthalpy   [J·kg⁻¹]

% ρ = ρ(P,h)  ⇒  total derivative  dρ/dP = (∂ρ/∂P)|h + (∂ρ/∂h)|P · (∂h/∂P)
% First: hold h constant, vary P  → (∂ρ/∂P)|h           [kg·m⁻³·bar⁻¹]
d_rho_sh_dP = fdiff(@(Pb) XSteam('rho_ph',Pb, h_sh/1e3), P_s, dP);
d_rho_sc_dP = fdiff(@(Pb) XSteam('rho_ph',Pb, h_sc/1e3), P_s, dP);

% Second: hold P constant, vary h  → (∂ρ/∂h)|P          [kg·m⁻³·(J·kg⁻¹)⁻¹]
d_rho_sh_dh = fdiff(@(hh) XSteam('rho_ph',P_s, hh /1e3), h_sh, dh);
d_rho_sc_dh = fdiff(@(hh) XSteam('rho_ph',P_s, hh /1e3), h_sc, dh);

% Chain-rule contributions to dρ/dP via h(P) curve      [kg·m⁻³·bar⁻¹]
chain_sh = d_rho_sh_dh .* dhV_dP;   % for super-heated region
chain_sc = d_rho_sc_dh .* dhL_dP;   % for sub-cooled  region

%--------------------------------------------------------------------
% 3) Two-phase (vapour-liquid) mixture slopes at saturation
%--------------------------------------------------------------------
% Apparent density  ρ_sat = γ·ρ_V + (1-γ)·ρ_L
d_rho_sat_dP = gamma_bar .* d_rhoV_dP + (1-gamma_bar) .* d_rhoL_dP;   % [kg·m⁻³·bar⁻¹]

% Mixture enthalpy   h_sat = γ·h_V + (1-γ)·h_L
dh_sat_dP    = gamma_bar .* dhV_dP    + (1-gamma_bar) .* dhL_dP;      % [J·kg⁻¹·bar⁻¹]


%====================================================================
% MASS-BALANCE COEFFICIENTS — SECONDARY SIDE
%   Units check:
%     A_s          [m²]
%     ρ            [kg·m⁻³]
%     L_i          [m]
%     dρ/dP        [kg·m⁻³·bar⁻¹]
%     dρ/dh        [kg·m⁻³·(J·kg⁻¹)⁻¹]
%   ⇒  kg (for L-terms)   or   kg·bar⁻¹ / kg·(J·kg⁻¹)⁻¹ for P/h terms
%
%   The global balance is
%     MBL1Sec*dL1_dt + MBL2Sec*dL2_dt + MBPsSec*dP_s_dt +
%     MBhoSec*dh_out_dt  =  m_dot_in − m_dot_out
%====================================================================

%---------------------------
% 1) Sub-cooled region (SC)
%---------------------------
MBL1SC2  = A_s * (rho_sc - rho_L);                         % ∂m/∂L1
MBPsSC2  = A_s * L1 * ( 0.5*d_rho_sc_dh*dhL_dP + d_rho_sc_dP );
%              └─── ∂ρ/∂h ⋅ ∂h/∂P  +  ∂ρ/∂P  →  ∂ρ/∂P  (kg·m⁻³·bar⁻¹)

%---------------------------
% 2) Two-phase saturation zone
%---------------------------
MBL1SAT2 = A_s * (rho_L - rho_V);                          % ∂m/∂L1
MBL2SAT2 = A_s * (1-gamma_bar) * (rho_L - rho_V);          % ∂m/∂L2
MBPsSAT2 = A_s * L2 * d_rho_sat_dP;                        % ∂m/∂P_s

%---------------------------
% 3) Super-heated region (SH)
%---------------------------
MBL1SH2  = A_s * (rho_V - rho_sh);                         % ∂m/∂L1
MBL2SH2  = A_s * (rho_V - rho_sh);                         % ∂m/∂L2
MBPsSH2  = A_s * L3 * ( 0.5*d_rho_sh_dh*dhV_dP + d_rho_sh_dP );
MBhoSH2  = A_s * L3 * 0.5 * d_rho_sh_dh;                   % ∂m/∂h_out

%---------------------------
% 4) Aggregate coefficients used in the ODE row
%---------------------------
MBL1Sec = MBL1SC2  + MBL1SAT2 + MBL1SH2;   % multiplies  dL1/dt
MBL2Sec =               MBL2SAT2 + MBL2SH2;   % multiplies  dL2/dt
MBPsSec = MBPsSC2  + MBPsSAT2 + MBPsSH2;   % multiplies  dP_s/dt
MBhoSec = MBhoSH2;                         % multiplies  dh_out/dt

%---------------------------
% 5) Source/sink term  (kg·s⁻¹)
%---------------------------
MBSec = (m_dot_in - m_dot_out);

% FINAL ODE ROW (to be placed later):
%   MBL1Sec*dL1_dt + MBL2Sec*dL2_dt + MBPsSec*dP_s_dt + MBhoSec*dh_out_dt = MBSec

%===========================================================
% ENERGY BALANCE — SECONDARY SIDE : Sub-cooled region (SC)
%
% Governing ODE row
%   EBL1SC2 * dL1_dt   +   EBPsSC2 * dP_s_dt   =   EBSC2
%
% where
%   L1   : length of the SC zone                     [m]
%   P_s  : secondary pressure                        [bar]
%===========================================================

% ----------------------------------------------------------
% Coefficient on  dL1/dt   →  ∂(A_s ρ_sc h_sc L1) / ∂L1
% ----------------------------------------------------------
EBL1SC2 = A_s * rho_sc * (h_sc - h_L);                       % [J·m⁻¹]

% ----------------------------------------------------------
% Coefficient on  dP_s/dt  (chain rule for ρ_sc and h_sc)
%   Matches  A_s L1 [ (∂ρ/∂P + ½ ∂ρ/∂h ∂h_f/∂P)(h_sc-h_L)
%                      + ½ ρ_sc ∂h_f/∂P  − 1 ]
% ----------------------------------------------------------
EBPsSC2 = A_s * L1 * ( 0.5 * rho_sc * dhL_dP + (h_sc - h_L) * ( 0.5 * d_rho_sc_dh * dhL_dP + d_rho_sc_dP ) - 1 );                                                         % ½ ρ_sc ∂h_f/∂P % (h_sc-h_L)% dimensionless

% ----------------------------------------------------------
% Right-hand-side term (W)
%   • enthalpy flow  (ṁ_in)(h_in – h_L)
%   • tube-wall heat flux into SC zone
% ----------------------------------------------------------
EBSC2 = m_dot_in * (h_in - h_L) ...                          % [W]
      + n * pi_val * D_i * alpha_i * L1 * (T_m3 - T_sc);     % [W]

% Final ODE contribution:
%   EBL1SC2 * dL1_dt  +  EBPsSC2 * dP_s_dt  =  EBSC2
%==============================================================
% ENERGY BALANCE — Secondary side, Saturated zone
%==============================================================

% 1) Coefficients on dL1/dt and dL2/dt
EBL1SAT2 = A_s * ( rho_sc*h_L - rho_sh*h_V );                       % [J m⁻¹]
EBL2SAT2 = A_s * ( (1-gamma_bar)*(rho_L*h_L - rho_V*h_V) ...
                 + h_V*(rho_V - rho_sh) );                          % [J m⁻¹]

% 2) Coefficient on dP_s/dt  (chain-rule contributions)
%------------------------------------------------------------------
% Energy balance – saturated (two-phase) zone
%   PDF formula:  A_s [ … ]   ←→   EBPsSAT2
%------------------------------------------------------------------

EBPsSAT2 = A_s * ( ...
      h_L * L1 * ( d_rho_sc_dP + 0.5 * d_rho_sc_dh * dhL_dP ) ...            % h_f L1 ( ∂ρ₁/∂P + ½ ∂ρ₁/∂h ∂h_f/∂P )
    + L2 * ( gamma_bar * d_rhoh_V_dP + (1-gamma_bar) * d_rhoh_L_dP - 1 ) ... % L2 ( γ ∂(ρ_g h_g)/∂P + (1-γ) ∂(ρ_f h_f)/∂P – 1 )
    + h_V * L3 * ( d_rho_sh_dP + 0.5 * d_rho_sh_dh * dhV_dP ) );             % h_g L3 ( ∂ρ₃/∂P + ½ ∂ρ₃/∂h ∂h_g/∂P )


% 3) Coefficient on dh_out/dt
EBhoSAT2 = A_s * L3 * 0.5 * h_V * d_rho_sh_dh;                      % [J bar⁻¹]

% 4) Right-hand side (W)
EBSAT2 = m_dot_in*h_L - m_dot_out*h_V ...
       + n*pi_val*D_i*alpha_i*L2*(T_m2 - T_sat);
% Final ODE row:
%   EBL1SAT2*dL1_dt + EBL2SAT2*dL2_dt + EBPsSAT2*dP_s_dt + EBhoSAT2*dh_out_dt = EBSAT2

%==============================================================
% ENERGY BALANCE — Secondary side, Super-heated (SH) region
%==============================================================

% -- Coefficients on dL1/dt and dL2/dt ------------------------
EBL1SH2 = A_s * rho_sh * (h_V - h_sh);          % [J m⁻¹]
EBL2SH2 = A_s * rho_sh * (h_V - h_sh);          % [J m⁻¹]

% -- Coefficient on dP_s/dt  (chain-rule terms) ---------------
EBPsSH2 = A_s * L3 * ( ...
              0.5 * rho_sh * dhV_dP ...                                  % ½ ρ_sh ∂h_g/∂P_s
            + (h_sh - h_V) * ( 0.5 * d_rho_sh_dh * dhV_dP ...            % + (h_sh−h_g)(…)
                               + d_rho_sh_dP ) ...
            - 1 );                                                       % −1  (per PDF)

% -- Coefficient on dh_out/dt ---------------------------------
EBhoSH2 = 0.5 * A_s * L3 * ( (h_sh - h_V) * d_rho_sh_dh + rho_sh );      % [kg]

% -- Right-hand-side term (W) ---------------------------------
EBSH2 = m_dot_out * (h_V - h_out) ...                                    % NOTE:  m_dot_out
       + n * pi_val * D_i * alpha_i * L3 * (T_m1 - T_sh);

% Governing ODE row:
%   EBL1SH2*dL1_dt + EBL2SH2*dL2_dt + EBPsSH2*dP_s_dt + EBhoSH2*dh_out_dt = EBSH2

%--------------------------------------------------------------
% ENERGY BALANCE – Tube metal of SC segment (wall at T_m3)
%--------------------------------------------------------------

% Coefficient on dL1/dt
%   ∂/∂L1 [ ρ_w c_w (n A_w) (T_m3 − T_m2) ]      [J m⁻¹]
EBL1TMSC = n * rho_w * c_w * A_w * (T_m3 - T_m2);

% Coefficient on dT_m3/dt
%   ∂/∂T_m3 [ ρ_w c_w (n A_w) L1 T_m3 ]           [J K⁻¹]
EBTmTMSC = n * rho_w * c_w * A_w * L1;

% Right-hand-side: net heat flux into wall #3 (W)
EBTMSC = ...
      eta_SG * n * pi_val * D_o * alpha_o * L1 * (T_p5 - T_m3) ...   % from primary
    +          n * pi_val * D_i * alpha_i * L1 * (T_sc - T_m3);      % to secondary

% Governing ODE row
%   EBL1TMSC * dL1_dt   +   EBTmTMSC * dT_m3_dt   =   EBTMSC
%--------------------------------------------------------------
% ENERGY STORAGE COEFFICIENT on dT_m2/dt
%   C_m2 = n · ρ_w c_w A_w L2                 [J K⁻¹]
%   (n accounts for all tubes)
%--------------------------------------------------------------
EBTmTMSAT = n * rho_w * c_w * A_w * L2;

%--------------------------------------------------------------
% RIGHT-HAND-SIDE: net heat flux into wall #2  (W)
%   • from primary coolant at T_p3
%   • to saturated secondary at T_sat
%--------------------------------------------------------------
EBTMSat = ...
      eta_SG * n * pi_val * D_o * alpha_o * L2 * (T_p3 - T_m2) ...
    +          n * pi_val * D_i * alpha_i * L2 * (T_sat - T_m2);

% Governing ODE row (no moving boundary term here):
%   EBTmTMSAT * dT_m2_dt  =  EBTMSat

%------------------------------------------------------------------
% ENERGY STORAGE COEFFICIENTS
%   Wall #1 spans the SH region whose length changes with both L1 & L2
%------------------------------------------------------------------

% ∂/∂L1  [ n ρ_w c_w A_w (T_m2 − T_m1) ]   →  multiplies dL1/dt
EBL1TMSH = n * rho_w * c_w * A_w * (T_m2 - T_m1);   % [J m⁻¹]

% ∂/∂L2  [ n ρ_w c_w A_w (T_m2 − T_m1) ]   →  multiplies dL2/dt
EBL2TMSH = n * rho_w * c_w * A_w * (T_m2 - T_m1);   % [J m⁻¹]

% ∂/∂T_m1[ n ρ_w c_w A_w L3 T_m1 ]        →  multiplies dT_m1/dt
EBTmTMSH = n * rho_w * c_w * A_w * L3;             % [J K⁻¹]

%------------------------------------------------------------------
% RIGHT-HAND SIDE  (net heat flux into wall #1)      [W]
%------------------------------------------------------------------
EBTMSH = ...
      eta_SG * n * pi_val * D_o * alpha_o * L3 * (T_p1 - T_m1) ...  % from primary
    +          n * pi_val * D_i * alpha_i * L3 * (T_sh - T_m1);     % to secondary

% Governing ODE row:
%   EBL1TMSH*dL1_dt  +  EBL2TMSH*dL2_dt  +  EBTmTMSH*dT_m1_dt  =  EBTMSH
% ---------------------------------------------------------------
% ENERGY Balance / MOVING-BOUNDARY COEFFICIENTS
%   Primary coolant flows shell-side, Sub-ccoled region.
% ---------------------------------------------------------------

EBL1SC1  = rho_p5 * Cp_p5 * A_p * (T_p5 - T_p3);   % [J m⁻¹]
EBTp5SC1 = rho_p5 * Cp_p5 * A_p * L1;              % [J K⁻¹]

% ---------------------------------------------------------------
% RIGHT-HAND SIDE  (net heat into ∆x = L1 segment)  [W]
%   • from tube wall:   n π D_o α_o L1 (T_m3 − T_p5)
%   • by advection:     ṁ_p (h_p4 − h_p6)
% ---------------------------------------------------------------
EBSC1 = ...
    n * pi_val * D_o * alpha_o * L1 * (T_m3 - T_p5) ...
  + m_p * (h_p4 - h_p6);

% Governing ODE row:
%   EBL1SC1 * dL1_dt  +  EBTp5SC1 * dT_p5_dt  =  EBSC1

% --------------------------------------------------------------
% ENERGY BALANCE – Primary side, Saturated zone
%   Control volume between T_p2  (inlet)  and  T_p4  (outlet)
% --------------------------------------------------------------

EBL1SAT1  = 0;                                        % coefficient on dL1/dt
EBL2SAT1  = 0;                                        % coefficient on dL2/dt
EBTp3SAT1 = rho_p3 * Cp_p3 * A_p * L2;               % [J K⁻¹]  × dT_p3/dt

EBSAT1 = ...
      n * pi_val * D_o * alpha_o * L2 * (T_m2 - T_p3) ...   % heat from tube wall
    + m_p * (h_p2 - h_p4);                                  % enthalpy advection  (p2→p4)

% Governing ODE row
%   0*dL1_dt + 0*dL2_dt + EBTp3SAT1*dT_p3_dt  =  EBSAT1
%------------------------------------------------------------------
% ENERGY BALANCE – Primary side, Super-heated zone (L3)
%   Control volume between  T_HL  (inlet)  and  T_p2  (outlet)
%   Average state tracked by T_p1 = ½(T_HL + T_p2)
%------------------------------------------------------------------

% -- moving-boundary coefficients --------------------------------
EBL1SH1 = rho_p1 * Cp_p1 * A_p * (T_p3 - T_p1);   % ∂U/∂L1  [J m⁻¹]
EBL2SH1 = rho_p1 * Cp_p1 * A_p * (T_p3 - T_p1);   % ∂U/∂L2  [J m⁻¹]

% -- storage coefficient on dT_p1/dt -----------------------------
EBTp1SH1 = rho_p1 * Cp_p1 * A_p * L3;             % [J K⁻¹]

% -- right-hand side (heat in – heat out) ------------------------
EBSH1 = ...
      n * pi_val * D_o * alpha_o * L3 * (T_m1 - T_p1) ...   % heat from tube wall
    + m_p * (h_HL - h_p2);                                  % enthalpy advection

% Governing ODE:
%   EBL1SH1*dL1_dt  +  EBL2SH1*dL2_dt  +  EBTp1SH1*dT_p1_dt  =  EBSH1


  %--- 2) analytic derivatives ---
  dphi_dt  = ((react - beta)/Lambda)*phi + lambda_val*C;
  dC_dt    = (beta/Lambda)*phi - lambda_val*C;
  dT_f_dt  = (tau*P_0*phi + h_fc*A_fc*(T_c1-T_f))/(m_f*c_pf);
  dT_c1_dt = ((1-tau)*P_0*phi + h_fc*A_fc*(T_f-T_c1))/(m_c*Cp_core) + m_p*(h_CL-h_c2)/(m_c*Cp_core);
  dT_HL_dt = (T_c2 - T_HL)/tau_HL;
  dT_CL_dt = (T_p6 - T_CL)/tau_CL;

  core_dot = [ dphi_dt;
             dC_dt;
             dT_f_dt;
             dT_c1_dt;
             dT_HL_dt;
             dT_CL_dt ];

  %% --- Steam-generator coefficient matrix -----------------------
% Row / column order:
% [ L1̇  L2̇  Ṗs  ḣout  Ṫm3  Ṫm2  Ṫm1  Ṫp5  Ṫp3  Ṫp1 ]

M = [ ...
  EBL1SC2 ,     0     , EBPsSC2 ,     0     , 0 , 0 , 0 , 0 , 0 , 0 ;   % SC-energy
  EBL1SAT2, EBL2SAT2 , EBPsSAT2, EBhoSAT2 , 0 , 0 , 0 , 0 , 0 , 0 ;     % SAT-energy
  EBL1SH2 , EBL2SH2  , EBPsSH2 , EBhoSH2  , 0 , 0 , 0 , 0 , 0 , 0 ;     % SH-energy
  MBL1Sec , MBL2Sec  , MBPsSec , MBhoSec  , 0 , 0 , 0 , 0 , 0 , 0 ;     % mass balance
  EBL1TMSC,     0     ,    0    ,     0     , EBTmTMSC , 0 , 0 , 0 , 0 , 0 ;  % wall-3
       0   ,     0     ,    0    ,     0     ,    0    , EBTmTMSAT , 0 , 0 , 0 , 0 ;  % wall-2
  EBL1TMSH, EBL2TMSH ,    0    ,     0     ,    0    , 0 , EBTmTMSH , 0 , 0 , 0 ;     % wall-1
  EBL1SC1 ,     0     ,    0    ,     0     ,    0    , 0 , 0 , EBTp5SC1 , 0 , 0 ;     % primary SC
       0   ,     0     ,    0    ,     0     ,    0    , 0 , 0 , 0 , EBTp3SAT1 , 0 ;   % primary SAT
  EBL1SH1 , EBL2SH1  ,    0    ,     0     ,    0    , 0 , 0 , 0 , 0 , EBTp1SH1 ];    % primary SH

% Build the RHS vector R (same row order):
R = [EBSC2; EBSAT2; EBSH2; MBSec; EBTMSC; EBTMSat; EBTMSH; ...
     EBSC1; EBSAT1; EBSH1];

% check conditioning, then solve
%--------------------------------------------------------------
% MATLAB’s mldivide issues a warning (not an error) when M is singular,
% so we test the reciprocal condition number first.
rc = rcond(M);              % ≈ 1/cond(M).  Well-conditioned if rc ≳ 1e-12.

if rc < 1e-14       % threshold can be relaxed or tightened
    error('reactor_ode4:SingularSGMatrix', ...
          'SG coefficient matrix is singular/ill-conditioned (rcond = %.2e) at t = %.4f s', ...
          rc, t);
end

% Safe to solve
ySG_dot = M \ R;     % 10×1 column [L1̇; L2̇; Ṗs; ḣouṫ; Ṫm3; Ṫm2; Ṫm1; Ṫp5; Ṫp3; Ṫp1]

% Optional sanity-check: NaNs will still slip through if R contains NaN
if any(~isfinite(ySG_dot))
    error('reactor_ode4:NaNinSG', 'Non-finite SG derivative at t = %.4f s', t);
end

% Assemble the full 16-element derivative vector for the integrator:
dy          = zeros(16,1);
dy(1:6)     = core_dot;      % your six non-SG ODEs
dy(7:16)    = ySG_dot;       % SG derivatives in matching order
