% ********************************************************
%    * Optimization based Neuro Adaptive Controller *
%
%   Crafted By - Ryu Myeongseok
%   Version 0.1
%   Date : 2024.08.-  
%
% ********************************************************
 
%% INITIALIZATION
% ********************************************************
clear
clc

addpath("dynamics_model/");
addpath("utils/");
addpath("utils_NN/");

%% SIMULATION SETTING
% ********************************************************
dt = 1e-4;
% T = 5;
T = 10;
% T = 15;
% T = 25;
t = 0:dt:T;
rpt_dt = 0.5;

seed = 1; rng(seed);

ctrl_name = "CTRL1"; % Backstepping
% ctrl_name = "CTRL2"; % Dixon
% ctrl_name = "CTRL3"; % Proposed no ball constraint
% ctrl_name = "CTRL4"; % Proposed

ANIMATION_FLAG = 0;
AINMATION_SAVE_FLAG = 0;
FIGURE_SAVE_FLAG = 0;
RESULT_SAVE_FLAG = 1;

fprintf("                                   \n");
fprintf("***********************************\n");
fprintf("*    Optimization based NN Ctrl   *\n");
fprintf("***********************************\n");

fprintf("                                 \n");
fprintf(" Controller:\n");
fprintf("     "+ctrl_name+"\n");
fprintf("\n");
fprintf(" FLAGs\n");
fprintf("     animation      : "+string(ANIMATION_FLAG)+"\n");
fprintf("     animation save : "+string(AINMATION_SAVE_FLAG)+"\n");
fprintf("     figure save    : "+string(FIGURE_SAVE_FLAG)+"\n");
fprintf("     result save    : "+string(RESULT_SAVE_FLAG)+"\n");
fprintf("\n");

%% SYSTEM DECLARE
% ********************************************************
% [grad_x, xd, IC] = system1();
% [grad_x, r_f, rd_f, rdd_f, IC] = system2();
% [grad_x, r_f, rd_f, rdd_f, IC] = system3();
[grad_x, r_f, rd_f, rdd_f, IC] = model1_load();

%% INITIAL CONDITION
% ********************************************************
x = IC.x;
u = IC.u;

%% NUERAL NETWORK DECLARE
% ********************************************************
nnOpt = nn_opt(ctrl_name, dt);
nn = nn_init(nnOpt);

%% RECORDER
% ********************************************************
num_x = length(x); num_u = length(u);
num_t = length(t);

x_hist = zeros(num_x, num_t); x_hist(:, 1) = x;
r_hist = zeros(num_x/2, num_t); r_hist(:, 1) = r_f(x, 0);
rd_hist = zeros(num_x/2, num_t); rd_hist(:, 1) = rd_f(x, 0);
rdd_hist = zeros(num_x/2, num_t); rdd_hist(:, 1) = rdd_f(x, 0);
u_hist = zeros(num_u, num_t); u_hist(:, 1) = u;
uSat_hist = zeros(num_u, num_t); uSat_hist(:, 1) = u;
uLPF_hist = zeros(num_u, num_t); uLPF_hist(:, 1) = u;
L_hist = zeros(8, num_t); L_hist(:,1) = zeros(8,1);
dot_L_hist = zeros(1, num_t); dot_L_hist(:,1) = 0;
V_hist = zeros(nnOpt.l_size-1, num_t); V_hist(:, 1) = nn_V_norm_cal(nn.V, nnOpt);

%% MAIN SIMULATION
% ********************************************************
fprintf("[INFO] Simulation Start\n");

% backstep check
k1 = nnOpt.ctrl_param.k1;
k2 = nnOpt.ctrl_param.k2;
M = nnOpt.ctrl_param.M;
C = nnOpt.ctrl_param.C;
G = nnOpt.ctrl_param.G;

u_LPF = zeros(2,1);

for t_idx = 2:1:num_t
    r1 = r_f(x, t(t_idx));
    r1d = rd_f(x, t(t_idx));
    r1dd = rdd_f(x, t(t_idx));
  
    % backstep check
    % e = x-[r;rd];
    e1 = x(1:2) - r1;
    r2 = r1d - k1*e1;
    e2 = x(3:4) - r2;
    e = [e1;e2];
    r2d = r1dd - k1*(x(3:4) - r1d);
   
    % control input
    if strcmp(ctrl_name, "CTRL1")
        % k1 = 1e2; k2 = 1e1;
        u = -M*k2*e2 -M*e1 + C*x(3:4) + G +M*r2d;
        % u = -[k1 0 k1 0 ; 0 k2 0 k2] * [e1;e2];
        dot_L = 0;
    else
        x_in = [r1];
        [nn, u_NN, info] = nn_forward(nn, nnOpt, x_in);
        [nn, nnOpt, dot_L, info] = nn_backward(nn, nnOpt, e, u_NN);

        % u_NN = zeros(2,1); 
        u = -u_NN;
    end

    LPF_tau = 1/(2*pi * 1000);
    LPF_Ts = 1e-4;
    % LPF_Ts = 1;

    u_LPF = (LPF_tau * u_LPF + LPF_Ts * u) / (LPF_tau + LPF_Ts);

    % if norm(u)^2 - nnOpt.cstr.u_ball^2 > 0
    if norm(u_LPF)^2 - nnOpt.cstr.u_ball^2 > 0
        u_sat = u_LPF/norm(u_LPF) * nnOpt.cstr.u_ball;
        % u_sat = u/norm(u) * nnOpt.cstr.u_ball;
    else
        % u_sat = u;
        u_sat = u_LPF;
    end

    % gradient calculation
    grad = grad_x(x, u_sat, t(t_idx));

    % error check
    assert(~isnan(norm(u)));
    assert(~isnan(norm(grad)));

    % step forward
    x = x + grad * dt;
    
    % record
    x_hist(:, t_idx) = x;
    r_hist(:, t_idx) = r1;
    rd_hist(:, t_idx) = r2;
%     rdd_hist(:, t_idx) = rdd;
    u_hist(:, t_idx) = u;
    uSat_hist(:, t_idx) = u_sat;
    uLPF_hist(:, t_idx) = u_LPF;
    if strcmp(nnOpt.alg,"Proposed")
        L_hist(:, t_idx) = nnOpt.Lambda;
    end
    V_hist(:, t_idx) = nn_V_norm_cal(nn.V, nnOpt);
    dot_L_hist(:, t_idx) = dot_L;

    % simulation report
    if rem(t(t_idx)/dt, rpt_dt/dt) == 0
        fprintf("[INFO] Time Step %.2f/%.2fs (%.2f%%)\r", ...
            t(t_idx), T, t(t_idx)/T*100);
    end
    
end
fprintf("[INFO] Simulation End\n");
fprintf("\n");

%% PLOT
% ********************************************************
fprintf("[INFO] Plotting...\n");
plot_result

%% ADDITAIONAL FUNCTIONS
% ********************************************************
if ANIMATION_FLAG
    fprintf("[INFO] Animation Generating...\n");

    animate
end

if FIGURE_SAVE_FLAG
    fprintf("[INFO] Figure Saving...\n");

    [~, ~] = mkdir("sim_result/"+ctrl_name);
    
    for idx = 1:1:fig_len
        f_name = "sim_result/"+ctrl_name+"/"+ctrl_name+"_fig"+string(idx);

        saveas(figure(idx), f_name+".png");
        exportgraphics(figure(idx), f_name+'.eps')
    end
end

if RESULT_SAVE_FLAG
    fprintf("[INFO] Result Saving...\n");

    save("sim_result/"+ctrl_name+"/"+ctrl_name+"_result.mat", ...
        "t", "x_hist", "r_hist", "rd_hist", "nnOpt", ...
        "u_hist", "uSat_hist", "L_hist", "V_hist", "dot_L_hist" ...
        );
end

% save_weights

%%
% ********************************************************
fprintf("\n");
fprintf("[INFO] Program Terminated\n");
beep()







