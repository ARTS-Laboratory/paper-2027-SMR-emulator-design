function reactorBlock(block)
    setup(block);
end

function setup(block)
    % Ports
    block.NumInputPorts  = 5;    
    % CHANGE 1: Increased Output Ports from 20 to 21
    block.NumOutputPorts = 21;   
    block.SetPreCompPortInfoToDefaults;

    % Direct feedthrough flags:
    % We allow DF for 1..3 (used in dynamics); break DF for 4..5 to avoid algebraic loop
    block.InputPort(1).DirectFeedthrough = true;   % T_in
    block.InputPort(2).DirectFeedthrough = true;   % rho_ext (Needed for Reactivity Output)
    block.InputPort(3).DirectFeedthrough = true;   % Valve
    block.InputPort(4).DirectFeedthrough = false;  % m_dot_out1 (cached)
    block.InputPort(5).DirectFeedthrough = false;  % m_dot_in1  (cached)

    % States & sample time
    block.NumContStates = 16;
    block.SampleTimes   = [0 0];   % continuous

    % Register methods
    block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
    block.RegBlockMethod('InitializeConditions', @InitCond);
    block.RegBlockMethod('Outputs',              @Output);
    block.RegBlockMethod('Derivatives',          @Derivative);
    block.RegBlockMethod('Update',               @Update);
end

function DoPostPropSetup(block)
    % Create DWorks here (required by Simulink)
    block.NumDworks = 2;
    block.Dwork(1).Name            = 'mdot_out_last';
    block.Dwork(1).Dimensions      = 1;
    block.Dwork(1).DatatypeID      = 0;      
    block.Dwork(1).Complexity      = 'Real';
    block.Dwork(1).UsedAsDiscState = true;
    block.Dwork(2).Name            = 'mdot_in_last';
    block.Dwork(2).Dimensions      = 1;
    block.Dwork(2).DatatypeID      = 0;      
    block.Dwork(2).Complexity      = 'Real';
    block.Dwork(2).UsedAsDiscState = true;
end

function InitCond(block)
    % Initialize continuous states
    S = load('steadyStateResults.mat','X_dyn0');
    block.ContStates.Data = S.X_dyn0;

    % Initialize cached flows 
    block.Dwork(1).Data = 67.07;  % mdot_out_last 
    block.Dwork(2).Data = 67.07;  % mdot_in_last   
end

function Update(block)
    % Cache the latest inputs (breaks direct feedthrough on ports 4 & 5)
    block.Dwork(1).Data = block.InputPort(4).Data;  % m_dot_out1
    block.Dwork(2).Data = block.InputPort(5).Data;  % m_dot_in1
end

function Output(block)
    % 1..16: pass out states 
    y = block.ContStates.Data;
    for k = 1:16
        block.OutputPort(k).Data = y(k);
    end

    % Define indices for readability
    idx = struct('phi',1,'C',2,'T_f',3,'T_c1',4,'T_HL',5,'T_CL',6,...
                 'L1',7,'L2',8,'P_s',9,'h_out',10,'T_m3',11,'T_m2',12,...
                 'T_m1',13,'T_p5',14,'T_p3',15,'T_p1',16);

    % CHANGE 2: Read Input Port 2 (rho_ext)
    % Because DirectFeedthrough=true for Port 2, this is allowed here.
    rho_ext = block.InputPort(2).Data;

    % Load parameters
    % (For better performance, you could cache this in a DWork or userdata, 
    % but calling it here ensures it matches the ODE logic)
    params = buildParams4();

    % CHANGE 3: Calculate Reactivity (Algebraic Equation)
    % react = rho_ext + alpha_f*(T_f-T_ref_f) + alpha_c*(T_c1-T_ref_c);
    T_f  = y(idx.T_f);
    T_c1 = y(idx.T_c1);
    
    react = rho_ext + params.alpha_f * (T_f - params.T_ref_f) + ...
                      params.alpha_c * (T_c1 - params.T_ref_c);

    % Existing Derived Calculations
    P_s_val   = y(idx.P_s);
    h_out_val = y(idx.h_out);
    m_dot_out = block.Dwork(1).Data;

    %m_p    = params.m_dot_p_rated * y(idx.phi)^(1/3);
    m_p    = params.m_dot_p_rated;
    
    S_turb = XSteam('s_ph', P_s_val, h_out_val/1000);
    h_iso  = XSteam('h_ps', 0.08, S_turb)*1000;
    P_mech = params.eta_T * (h_out_val - h_iso) * m_dot_out;
    T_out  = XSteam('T_ph', P_s_val, h_out_val/1000);

    % Assign Derived Outputs
    block.OutputPort(17).Data = m_dot_out;
    block.OutputPort(18).Data = m_p;
    block.OutputPort(19).Data = P_mech;
    block.OutputPort(20).Data = T_out;
    
    % CHANGE 4: Assign Reactivity to Port 21
    block.OutputPort(21).Data = react;
end

function Derivative(block)
    t = block.CurrentTime;
    y = block.ContStates.Data;
    persistent params;
    if isempty(params), params = buildParams4(); end
    
    % Read inputs 
    T_in        = block.InputPort(1).Data;
    rho_ext     = block.InputPort(2).Data;
    Valve       = block.InputPort(3).Data;
    m_dot_out1  = block.InputPort(4).Data;
    m_dot_in1   = block.InputPort(5).Data;
    
    % Time functions
    params.T_in_fun    = @(~) T_in;
    params.rho_ext_fun = @(~) rho_ext;
    params.Valve_fun   = @(~) Valve;
    
    % Pass to ODE
    dy = reactor_ode4(t, y, params, m_dot_out1, m_dot_in1);
    block.Derivatives.Data = dy;
end
