%% REPORTY SIMULATION
fprintf("                                   \n");
fprintf("***********************************\n");
fprintf("*    Optimization based NN Ctrl   *\n");
fprintf("***********************************\n");

fprintf("                                 \n");
fprintf(" Controller:\n");
fprintf("     "+string(c_idx)+"\n");
fprintf("\n");
fprintf(" FLAGs\n");
fprintf("     animation      : "+string(ANIMATION_FLAG)+"\n");
fprintf("     animation save : "+string(AINMATION_SAVE_FLAG)+"\n");
fprintf("     figure save    : "+string(FIGURE_SAVE_FLAG)+"\n");
fprintf("     result save    : "+string(RESULT_SAVE_FLAG)+"\n");
fprintf("\n");

%% SYSTEM DECLARE
[grad_x, IC] = model1_load();
[r_f, rd_f, rdd_f] = ref1_load(); % sin func
% [r_f, rd_f, rdd_f] = ref2_load(); % step func

%% INITIAL CONDITION
x = IC.x;
u = IC.u;
r1 = r_f(x, 0);
r2 = rd_f(x, 0);

%% NUERAL NETWORK DECLARE
nn = nn_init(nnOpt);


%% RECORDER
num_x = length(x); num_u = length(u);
num_t = length(t);

x_hist = zeros(num_x, num_t); x_hist(:, 1) = x;
r_hist = zeros(num_x/2, num_t); r_hist(:, 1) = r1;
rd_hist = zeros(num_x/2, num_t); rd_hist(:, 1) = r2;
u_hist = zeros(num_u, num_t); u_hist(:, 1) = u;
uSat_hist = zeros(num_u, num_t); uSat_hist(:, 1) = u;
L_hist = zeros(8, num_t); L_hist(:,1) = zeros(8,1);
dot_L_hist = zeros(1, num_t); dot_L_hist(:,1) = 0;
V_hist = zeros(nnOpt.l_size-1, num_t); V_hist(:, 1) = nn_V_norm_cal(nn.V, nnOpt);
aux_hist = zeros(num_x/2, num_t);


%% 
k1 = ctrlOpt.k1;
k2 = ctrlOpt.k2;
M = ctrlOpt.M;
if strcmp(ctrlOpt.type, "CTRL1")
    C = ctrlOpt.C;
    G = ctrlOpt.G;
    
elseif strcmp(ctrlOpt.type, "CTRL2")
    alg = nnOpt.alg;
    C = ctrlOpt.C;
    G = ctrlOpt.G;

elseif strcmp(ctrlOpt.type, "CTRL3")
    C = ctrlOpt.C;
    G = ctrlOpt.G;
    Azeta = ctrlOpt.Azeta;
    Bzeta = ctrlOpt.Bzeta;
    inv_M = ctrlOpt.inv_M;

elseif strcmp(ctrlOpt.type, "CTRL4")

else
    error()
end

%% MAIN SIMULATION
% ********************************************************
fprintf("[INFO] Simulation Start\n");

zeta = [0;0];
pre_zeta = zeta;
aux_hist(:,1) = zeta;   
z1 = 0; z2 = 0;

for t_idx = 2:1:num_t
    r1 = r_f(x, t(t_idx));
    r1d = rd_f(x, t(t_idx));
    r1dd = rdd_f(x, t(t_idx));
  
    % backstep check
    e1 = x(1:2) - r1;
    r2 = r1d - k1*e1;
    e2 = x(3:4) - r2;
    e = [e1;e2];
    r2d = r1dd - k1*(x(3:4) - r1d);
   
    % control input
    if strcmp(ctrlOpt.type, "CTRL1")
        u = -M*k2*e2 -M*e1 + C*x(3:4) + G +M*r2d;
        dot_L = 0;

    elseif strcmp(ctrlOpt.type, "CTRL2")
        x_in = r1;
        [nn, u_NN, info] = nn_forward(nn, nnOpt, x_in);
        [nn, nnOpt, dot_L, info] = nn_backward(nn, nnOpt, ctrlOpt, e, u_NN);

        % u_NN = zeros(size(u_NN));
        u1 = - M * k2 * e2 - M * e1 + C * x(3:4) + G + M * r2d -u_NN;
        u = u1;
        % u = -u_NN;

    elseif strcmp(ctrlOpt.type, "CTRL3")
        x_in = r1;
        [nn, u_NN, info] = nn_forward(nn, nnOpt, x_in);
        [nn, nnOpt, dot_L, info] = nn_backward(nn, nnOpt, ctrlOpt, e+[0;0;zeta], u_NN);

        % u_NN = zeros(size(u_NN));
        u = - M * k2 * e2 - M * e1 + C * x(3:4) + G + M * r2d -u_NN;
        % u = -u_NN;

    elseif strcmp(ctrlOpt.type, "CTRL4")
        x_in = r1;
        [nn, u_NN, info] = nn_forward(nn, nnOpt, x_in);
        [nn, nnOpt, dot_L, info] = nn_backward(nn, nnOpt, ctrlOpt, e, u_NN);
        
        u = -u_NN;
        
    elseif strcmp(ctrlOpt.type, "CTRL5")
        % it use other type of auxiliary control input
        l1 = 1; l2 = 1; l3 = 13; l4 = 13;
        % p1 = 11/9; p2 = 9/11;
        p1 = 1.25; p2 = 0.75;
        b1 = p1^sign(1-abs(e1(1)));
        b2 = p2^sign(1-abs(e1(2)));

        em = 0.01;
        lbd = 0.01;
        
        k_l = 3; l_l = 3;
        k_g = 1; l_g = 1; 

        % for aux
        sigma_g = [.5 .5];

        r2 = r1d - l1*sig_alp(z1,b1) - l2*sig_alp(z1,b2) - z1;

        e2 = x(3:4) - r2;
        e = [e1;e2];
        z1 = e1; z2 = e2 - zeta;

        b1 = p1^sign(1-abs(z1(1)));
        b2 = p2^sign(1-abs(z2(2)));
        b3 = p1^sign(1-abs(z2(1)));
        b4 = p2^sign(1-abs(z2(2)));
        a1_g = p1^sign(1-abs(zeta(1)));
        a2_g = p2^sign(1-abs(zeta(2)));

        x_in = r1;
        [nn, u_NN, info] = nn_forward(nn, nnOpt, x_in);
        [nn, nnOpt, dot_L, info] = nn_backward(nn, nnOpt, [0;0;z2], u_NN);
       
        % u_NN = zeros(size(u_NN));
        % u = -u_NN;
        u = M * ( -z1 -l3*sig_alp(z2,b3) -l4*sig_alp(z2,b4) - z2 -u_NN ...
            -sigma_g*zeta -k_g*sig_alp(zeta,a1_g) -l_g*sig_alp(zeta,a2_g) ...
            -em*sign(z2) -lbd*sign(z2));
        
    end

    % control input saturation
    p = 100;
    u_bar = ctrlOpt.u_ball;
    u_sat = u / norm(u) * u_bar * ( norm(u)/u_bar / (1+(norm(u)/u_bar)^p)^(1/p) ); 

    % Auxiliary Systems
    if strcmp(ctrlOpt.type, "CTRL3")
        u_max = (1/sqrt(2)*u_bar + u_bar) / 2;
        u_sat_CTRL3 = min(max(u, -u_max), u_max);

        u_del = u - u_sat_CTRL3;

        grad_zeta = -Azeta * zeta + Bzeta * u_del;

        zeta = zeta + grad_zeta * dt;

    elseif strcmp(ctrlOpt.type, "CTRL5")
        u_del = u_sat - u;

        grad_zeta = -k_g*sig_alp(zeta,a1_g) - l_g*sig_alp(zeta,a2_g) - sigma_g*zeta + inv_M * u_del;

        zeta = zeta + grad_zeta * dt;
    end

    % gradient calculation
    grad = grad_x(x, u_sat, t(t_idx));

    % error check
    assert(~isnan(norm(u)));
    assert(~isnan(norm(grad)));

    % step forward
    x = x + grad * dt;

    x_hist(:, t_idx) = x;
    r_hist(:, t_idx) = r1;
    rd_hist(:, t_idx) = r2;
    u_hist(:, t_idx) = u;
    uSat_hist(:, t_idx) = u_sat;
    if strcmp(nnOpt.alg,"Proposed")
        L_hist(:, t_idx) = nnOpt.Lambda;
    end
    V_hist(:, t_idx) = nn_V_norm_cal(nn.V, nnOpt);
    dot_L_hist(:, t_idx) = dot_L;

    aux_hist(:, t_idx) = zeta;

    % simulation report
    if rem(t(t_idx)/dt, rpt_dt/dt) == 0
        fprintf("[INFO] Time Step %.2f/%.2fs (%.2f%%)\r", ...
            t(t_idx), T, t(t_idx)/T*100);
    end
    
end
fprintf("[INFO] Simulation End\n");
fprintf("\n");

if FIGURE_PLOT_FLAG
    %% PLOT
    fprintf("[INFO] Plotting...\n");
    plot_wrapper

    if FIGURE_SAVE_FLAG
        fprintf("[INFO] Figure Saving...\n");

        [~, ~] = mkdir("sim_result/"+ctrl_name);
        
        for idx = 1:1:fig_len
            f_name = "sim_result/"+ctrl_name+"/"+ctrl_name+"_fig"+string(idx);

            saveas(figure(idx), f_name+".png");
            exportgraphics(figure(idx), f_name+'.eps')
        end
    end

    figure(7);clf;
    tiledlayout(2,1);

    nexttile
    plot(t, aux_hist(1,:), 'r', 'LineWidth', 2); hold on;
    grid on

    nexttile
    plot(t, aux_hist(2,:), 'b', 'LineWidth', 2); hold on;
    grid on
end


%% ADDITAIONAL FUNCTIONS
if ANIMATION_FLAG
    fprintf("[INFO] Animation Generating...\n");

    animate
end

if RESULT_SAVE_FLAG
    fprintf("[INFO] Result Saving...\n");

    [~,~] = mkdir("sim_result/"+ctrl_name);
    save("sim_result/"+ctrl_name+"/"+ctrl_name+"_result.mat", ...
        "t", "x_hist", "r_hist", "rd_hist", "nnOpt", "ctrlOpt",...
        "u_hist", "uSat_hist", "L_hist", "V_hist", "dot_L_hist" ...
        );
    clear("nnOpt", "ctrlOpt");
end

% save_weights

% LOCAL FUNCTIONS
function out = sig_alp(in, alp)
    out = abs(in).^alp .* sign(in);
end

