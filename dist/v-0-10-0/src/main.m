%% INITIALIZATION
clear  
clc

addpath("models");
sim_start_time = string(datetime('now','Format','d-MMM-y_HH-mm-ss'));

fprintf("                                   \n");
fprintf("***********************************\n");
fprintf("* Numerical Validation for CONAC   *\n");
fprintf("***********************************\n");

SAVE_RESULT_FLAG    = 1;

%% SIMULATION SETTING
ctrl_dt = 1/250;            % control sampling time
dt = ctrl_dt/4;                % simulation time step
rpt_dt = 4;                 % report time step (for printing simulation progress)
T = 11+12*3;                       % with warmup time
% T = 23;
t = 0:dt:T;

%%
CTRL_LIST = {
    struct("CTRL_NUM", 1, "OPT_NUM", 3, ...
        "SAVE_NAME", "NAC", ...
        "Color", "#808080"), ...
    struct("CTRL_NUM", 3, "OPT_NUM", 1, ...
        "SAVE_NAME", "AUX", ...
            "Color", "magenta"), ...
    struct("CTRL_NUM", 1, "OPT_NUM", 2, ...
            "SAVE_NAME", "CONAC_LOW", ...
                "Color", "cyan") ...
    struct("CTRL_NUM", 1, "OPT_NUM", 1, ...
        "SAVE_NAME", "CONAC_HIGH", ...
        "Color", "blue"), ...
};

dataSet = cell(length(CTRL_LIST), 1);

for ctrl_idx = 1:1:length(CTRL_LIST)
    seed = 100; rng(seed);

    CTRL_INFO = CTRL_LIST{ctrl_idx};
    data = struct();
    data.CTRL_INFO = CTRL_INFO;

    fprintf("Running Simulation for Control Method:\n");
    fprintf(" %s, Option: %d\n", ...
        CTRL_INFO.SAVE_NAME, CTRL_INFO.OPT_NUM);

    %% SYSTEM DECLARE
    grad_x = model1_load();     % system dynamics (grad_x = f(x,u,t))
    r_func = ref4_load();       % ref. signal 

    x1 = [deg2rad(-90);0];      % initial state (link angle)
    x2 = [0;0];                 % initial state (link angular velocity)
    u = [0;0];                  % initial control input (torque)   

    num_x = length(x1); num_u = length(u); num_t = length(t);
    data.x1_hist = zeros(num_x, num_t); data.x1_hist(:, 1) = x1;
    data.x2_hist = zeros(num_x, num_t); data.x2_hist(:, 1) = x2;
    data.xd1_hist = zeros(num_x, num_t); data.xd2_hist = zeros(num_x, num_t);
    data.xd3_hist = zeros(num_x, num_t); data.xd3_hist = zeros(num_x, num_t);
    [xd1,xd2] = r_func(0);
    data.xd1_hist(:, 1) = xd1; data.xd2_hist(:, 1) = xd2;
    data.r_hist = zeros(num_x, num_t); data.r_hist(:, 1) = zeros(num_x, 1);
    data.u_hist = zeros(num_u, num_t); data.u_hist(:, 1) = u;
    data.uSat_hist = zeros(num_u, num_t); data.uSat_hist(:, 1) = u;

    %% CONTROLLER LOAD
    if CTRL_INFO.CTRL_NUM == 1    % CoNAC
        ctrl_path = "CoNAC";
    elseif CTRL_INFO.CTRL_NUM == 2 % Aux.
        ctrl_path = "CoNAC-AUX";
    elseif CTRL_INFO.CTRL_NUM == 3 % Complex Aux.
        ctrl_path = "CoNAC-AUX-Comp";
    else
        error("Invalid CONTROL_NUM. Must be 1, 2, or 3.");
    end

    addpath("controllers/"+ctrl_path);
    addpath(genpath('controllers/'+ctrl_path));

    % opt = loadOpts(ctrl_dt);
    addpath("controllers");
    opt = loadGlobalOpts(ctrl_dt, CTRL_INFO.CTRL_NUM, CTRL_INFO.OPT_NUM);
    initControl;

    %% RECORDER INITIALIZATION
    if CTRL_INFO.CTRL_NUM == 1
        % CoNAC
        data.lbd_hist = zeros(length(opt.beta), num_t); data.lbd_hist(:,1) = zeros(length(opt.beta),1);
        data.th_hist = zeros(opt.l_size-1, num_t); data.th_hist(:, 1) = nnWeightNorm(nn.th, opt);
    elseif CTRL_INFO.CTRL_NUM == 2 || CTRL_INFO.CTRL_NUM == 3
        % Aux.
        data.th_hist = zeros(opt.l_size-1, num_t); data.th_hist(:, 1) = nnWeightNorm(nn.th, opt);
        data.zeta_hist = zeros(num_u, num_t); data.zeta_hist(:, 1) = z;
    end

    %% MAIN SIMULATION
    % ********************************************************
    fprintf("Simulation Start\n");

    for t_idx = 2:1:num_t
        [xd1,xd2,xd3] = r_func(t(t_idx));
        
        e1 = x1 - xd1;
        e2 = x2 - xd2;

        r = e2 + opt.Lambda * e1;       % filtered error

        if t_idx==2 || rem(t(t_idx)/dt, ctrl_dt/dt) == 0
            preControl
            postControl
        end

        % control input saturation
        u_sat = sat(u, opt, CTRL_INFO.CTRL_NUM);

        % error check
        if isnan(norm(u))
            warning("NaN detected in control input at time %.2f. Simulation stopped.", t(t_idx));
            break;
        end
        
        % step forward
        % grad = grad_x([x1;x2], u_sat, t(t_idx));
        % x1 = x1 + grad(1:2) * dt;
        % x2 = x2 + grad(3:4) * dt;
        k1 = grad_x([x1;x2], u_sat, t(t_idx));
        k2 = grad_x([x1;x2] + k1*dt/2, u_sat, t(t_idx) + dt/2);
        k3 = grad_x([x1;x2] + k2*dt/2, u_sat, t(t_idx) + dt/2);
        k4 = grad_x([x1;x2] + k3*dt, u_sat, t(t_idx) + dt);
        del_x = (k1 + 2*k2 + 2*k3 + k4)/6*dt;
        x1 = x1 + del_x(1:2);
        x2 = x2 + del_x(3:4);

        % 
        data.x1_hist(:, t_idx) = x1;
        data.x2_hist(:, t_idx) = x2;
        data.xd1_hist(:, t_idx) = xd1;
        data.xd2_hist(:, t_idx) = xd2;
        data.xd3_hist(:, t_idx) = xd3;
        data.u_hist(:, t_idx) = u;
        data.uSat_hist(:, t_idx) = u_sat;
        data.th_hist(:, t_idx) = nnWeightNorm(nn.th, opt);
        data.r_hist(:, t_idx) = r;
        if CTRL_INFO.CTRL_NUM == 1
            data.lbd_hist(:, t_idx) = opt.lbd;
        elseif CTRL_INFO.CTRL_NUM == 2 || CTRL_INFO.CTRL_NUM == 3
            data.zeta_hist(:, t_idx) = z;
        end

% cpx_list(t_idx) = cpx_;
% bsc_hist(:,t_idx) =BSC_e2;
% bsc_hist(:,t_idx) = u_fix;

        % simulation report
        if rem(t(t_idx)/dt, rpt_dt/dt) == 0
            fprintf("Time Step %.2f/%.2fs (%.2f%%)\r",t(t_idx), T, t(t_idx)/T*100);
        end
        
    end
    fprintf("Simulation End\n");
    fprintf("\n");

    dataSet{ctrl_idx} = data;
    dataSet{ctrl_idx}.opt = opt;
    dataSet{ctrl_idx}.T = T;
    dataSet{ctrl_idx}.t = t;
    dataSet{ctrl_idx}.SIM = true;

end

if SAVE_RESULT_FLAG
    save("RESULT.mat", "dataSet", "-v7.3");
end

plotter