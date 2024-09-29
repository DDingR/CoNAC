
% ********************************************************
%    * Optimization based Neuro Adaptive Controller *
%
%   Crafted By - Ryu Myeongseok
%   Version 0.1
%   Date : 2024.08.-  
%
% ********************************************************

%% INITIALIZATION
clear  
clc

addpath("dynamics_model/");
addpath("utils/");
addpath("utils_NN/");
addpath("InpSatFunc/");

%% SIMULATION SETTING
dt = 1e-4;
% T = 5;
T = 10;
% T = 15;
% T = 200;
t = 0:dt:T;
rpt_dt = 1;


FIGURE_PLOT_FLAG    = 0;
ANIMATION_FLAG      = 0;
AINMATION_SAVE_FLAG = 0;
FIGURE_SAVE_FLAG    = 0;
RESULT_SAVE_FLAG    = 1;

%% CONTROLLER SETTING
opts = readtable("ctrlOpts.csv",'delimiter', ',', 'TextType', 'string','ReadVariableNames',1, 'ReadRowNames', 1);

num_ctrl = size(opts, 2);

% for c_idx = 1:1:num_ctrl
% for c_idx = 3
% for c_idx = 12
for c_idx = [3 12]
% for c_idx = [1 2 3 12]
    seed = 1; rng(seed);
    % ctrl_name = opts.Properties.VariableNames{c_idx};
    ctrl_name = "CTRL"+string(c_idx);
    [ctrlOpt, nnOpt] = opt_load(c_idx, opts, dt);

    simulate
end

%% TERMINATION
fprintf("\n");
fprintf("[INFO] Program Terminated\n");
beep()
