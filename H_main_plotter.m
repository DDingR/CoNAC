clear
clc

%%
FIGURE_PLOT = 1;
% SAVE_FLAG = 1;
SAVE_FLAG = 1;
POSITION_FLAG = 1; % it will plot fiugures in the same position

addpath("InpSatFunc/")

%% DATA LOAD
% group1 = 1;
% group2 = 2;
% group3 = 3:7;
% group4 = 8:12;

for c_idx = 1:1:4
    fprintf("[INFO] Loading CTRL%d\n", c_idx);
    ctrl_results.("ctrl"+string(c_idx)) = load( ...
        "sim_result/CTRL"+string(c_idx)+"/CTRL" + string(c_idx) + "_result.mat" ...
    );
end
c_list = get(gca,'colororder');

group = [1 2 3 4];
quantitive_report(ctrl_results, group);

if ~FIGURE_PLOT
    return
end

%% FIGURE NAMES
fig_names = [
    "fig1" %  state vs ref; CTRL1
    "fig2" %  state vs ref; CTRL2
    "fig3" %  state vs ref; CTRL3
    "fig4" %  state vs ref; CTRL4
    "fig5" %  control input; CTRL1
    "fig6" %  control input; CTRL2
    "fig7" %  control input; CTRL3
    "fig8" %  control input; CTRL4
    "fig9" %  weight ball; CTRL2
    "fig10" %  weight ball; CTRL3
    "fig11" %  weight ball; CTRL4
    "fig12" %  multipliers; CTRL4
    "fig13" %  saturation function
    "fig14" %  total controls
    "fig15" %  total trajectories
    "fig16" %  control ball (focus)
];

% REPRESENTIVE CONTROLLER
ctrl1_rep = ctrl_results.("ctrl"+string(group(1)));
ctrl2_rep = ctrl_results.("ctrl"+string(group(2)));
ctrl3_rep = ctrl_results.("ctrl"+string(group(3)));
ctrl4_rep = ctrl_results.("ctrl"+string(group(4)));

%% FIGURE 1~4: STATE vs REFERENCE
figure(1); clf; 
plot_ref(ctrl1_rep, POSITION_FLAG); 
figure(2); clf;
plot_ref(ctrl2_rep, POSITION_FLAG); 
figure(3); clf; 
plot_ref(ctrl3_rep, POSITION_FLAG); 
figure(4); clf; 
plot_ref(ctrl4_rep, POSITION_FLAG); 

%% FIGURE 5~8: CONTROL INPUT
figure(5); clf;
plot_control(ctrl1_rep, POSITION_FLAG, 1);
figure(6); clf;
plot_control(ctrl2_rep, POSITION_FLAG, 2);
figure(7); clf; 
plot_control(ctrl3_rep, POSITION_FLAG, 3);
figure(8); clf; 
plot_control(ctrl4_rep, POSITION_FLAG, 4);

%% FIGURE 3: WEIGHT BALL
figure(9); clf
plot_weight_norm(ctrl2_rep, POSITION_FLAG, 2);

figure(10); clf
plot_weight_norm(ctrl3_rep, POSITION_FLAG, 3);

figure(11); clf
plot_weight_norm(ctrl4_rep, POSITION_FLAG, 4);

%% FIGURE 4: MULTIPLIERS
figure(12); clf
plot_multipliers(ctrl4_rep, POSITION_FLAG);

%% FIGURE 5: SATURATION FUNCTION
figure(13); clf
plot_saturation(POSITION_FLAG);

%% FIGURE 6: TOTAL CONTROLS
figure(14); clf
tl = tiledlayout(2,1);
plot_total_control(ctrl_results, group);

figure(15); clf
tl = tiledlayout(2,1);
plot_err(ctrl1_rep, group, c_list(1,:));
plot_err(ctrl2_rep, group, c_list(2,:));
plot_err(ctrl3_rep, group, c_list(3,:));
plot_err(ctrl4_rep, group, c_list(4,:));

figure(16);clf
plot_ball(ctrl_results, group, c_list);


%% SAVE FIGURES
if SAVE_FLAG
    for idx = 1:1:length(fig_names)
        fig_name = fig_names(idx);
        f_name = "figures/main_plot/" + fig_name;

        saveas(figure(idx), f_name + ".png")
        exportgraphics(figure(idx), f_name+'.eps')
    end
end

beep()

%% LOCAL FUNCTIONS ************************************
%           *** state vs ref ***
function [] = plot_ref(result, POSITION_FLAG)
    tl = tiledlayout(2, 1);

    t = result.t;

    r1 = result.r_hist(1,:);
    r2 = result.r_hist(2,:);
    x1 = result.x_hist(1,:);
    x2 = result.x_hist(2,:);

    font_size = 16;
    line_width = 2;
    lgd_size = 12;

    nexttile(1)
    plot(t, r1, "Color", "green", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "${q_d}$"); hold on
    plot(t, x1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$q$"); hold on
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("$q_{1}\ (\rm rad)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    % ylim([min(r1) * 0.75, max(r1) * 1.25])
    ylim([-0.5, 2.5])
    % ylim([-1.5, 4.5])
    lgd = legend;
    lgd.Orientation = 'Horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 

    nexttile(2)
    plot(t, r2, "Color", "green", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "${q_d}_2$"); hold on
    plot(t, x2, "Color", "blue", "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "${q}_2$"); hold on
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("$q_{2}\ (\rm rad)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    % ylim([min(r2) * 1.25, max(r2) * 0.75])
    ylim([-2.5, .5])
    % ylim([-2.5, 2.5])
    % lgd = legend;
    % lgd.Orientation = 'Horizontal';
    % lgd.Location = 'northeast';
    % lgd.Interpreter = 'latex';
    % lgd.FontSize = lgd_size; 

    if POSITION_FLAG
        set(gcf, 'Position',  [100, 100, 560, 280])
    end
end

%           *** state vs ref ***
function [] = plot_err(result, POSITION_FLAG, c)
    tl = tiledlayout(2, 1);

    t = result.t;

    x1 = result.x_hist(1,:);
    x2 = result.x_hist(2,:);
    r1 = result.r_hist(1,:);
    r2 = result.r_hist(2,:);

    font_size = 16;
    line_width = 2;
    lgd_size = 12;

    nexttile(1)
    plot(t, x1-r1, "Color", c, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$q$"); hold on
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("${e_1}_{(1)}\ (\rm rad)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    xlim([5.5 8.5])
    ylim([-1, 1] * 1e-1);
    lgd = legend;
    lgd.Orientation = 'Horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 

    nexttile(2)
    plot(t, x2-r2, "Color", c, "LineWidth", line_width, "LineStyle", "-."); hold on
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("${e_1}_{(2)}\ (\rm rad)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    % ylim([min(r2) * 1.25, max(r2) * 0.75])
    xlim([5.5 8.5])
    ylim([-1, 1] * 1e-1);
    % ylim([-2.5, 2.5])

    % if POSITION_FLAG
    %     set(gcf, 'Position',  [100, 100, 560, 250])
    % end
end

%           *** control input ***
function [] = plot_control(result, POSITION_FLAG, CTRL_NUM)
    tl = tiledlayout(2, 1);

    cstr = result.nnOpt.cstr;
    t = result.t;

    u1 = result.u_hist(1,:);
    u2 = result.u_hist(2,:);
    u1_sat = result.uSat_hist(1,:);
    u2_sat = result.uSat_hist(2,:);

    font_size = 16;
    line_width = 2;
    lgd_size = 12;
    
    nexttile
    plot(t, u1, "Color", "red", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\tau$"); hold on
    plot(t, u1_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", " $h(\tau)$"); hold on
    % plot(t(actset1), u1(actset1), "Color", "red", "LineWidth", line_width, "LineStyle", "."); hold on
    % scatter(t(actset1), u1(actset1), "Color", "red", "Marker","."); hold on
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("$\tau_{1}\ (\rm Nm)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    ylim([-50 * 1.25, 60])
%     ylim([-50 * 1.25, 50 * 1.25])
    % ylim([min(u1_sat) * 1.25, max(u1_sat) * 1.25])
    lgd = legend;
    lgd.Orientation = 'Horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size;    
    
    nexttile
    plot(t, u2, "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(t, u2_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("$\tau_{2}\ (\rm Nm)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    ylim([-4 6])
    % ylim([-10 * 1.25, 10 * 1.25])
    % ylim([-50 * 1.25, 50 * 1.25])
    % ylim([min(u2_sat) * 1.25, max(u2_sat) * 1.25])
    
    if CTRL_NUM ~= 1
    nexttile(1)
    rectangle('Position', [2, 0, 1, 200], 'EdgeColor', [0.8500    0.3250    0.0980], 'LineWidth', 2, 'LineStyle', '-.'); hold on
    rectangle('Position', [6.5, 0, 1, 70], 'EdgeColor', [0.8500    0.3250    0.0980], 'LineWidth', 2, 'LineStyle', '-.'); hold on
    
        axes('Position',[.78 .78 .2 .2])
        ang = 0:0.01:2*pi;
        plot(cstr.u_ball*cos(ang), cstr.u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
        plot(u1, u2, "color", 'red', "LineWidth", 2, "LineStyle", "-"); hold on
        plot(u1_sat, u2_sat, "color", 'blue', "LineWidth", 2, "LineStyle", "-"); hold on
        xlabel("$\tau_{1}$", "Interpreter", "latex")
        ylabel("$\tau_{2}$", "Interpreter", "latex")
        set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
        grid on 
        xlim([-cstr.u_ball*1.5, cstr.u_ball*1.5])
        ylim([-cstr.u_ball*1.5, cstr.u_ball*1.5])
        pbaspect([1 1 1])
        
        % if max(u1) > 51
            
            % axes('Position',[.44,.64,.13,.11])
            axes('Position',[.34,.64,.13,.11])
            plot(t, u1, "Color", "red", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\tau$"); hold on
            plot(t, u1_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$h(\tau)$"); hold on
            % plot(t(actset1), u1(actset1), "Color", "red", "LineWidth", line_width, "LineStyle", "."); hold on
            % scatter(t(actset1), u1(actset1), "Color", "red", "Marker","."); hold on
            set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
            grid on
            ylim([00, 200])
            xlim([2 3])
            % ylim([min(u1_sat) * 1.25, max(u1_sat) * 1.25])

            % axes('Position',[.75,.64,.13,.11])
            axes('Position',[0.68,0.64,0.13,0.11])
            plot(t, u1, "Color", "red", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\tau$"); hold on
            plot(t, u1_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$h(\tau)$"); hold on
            % plot(t(actset1), u1(actset1), "Color", "red", "LineWidth", line_width, "LineStyle", "."); hold on
            % scatter(t(actset1), u1(actset1), "Color", "red", "Marker","."); hold on
            set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
            grid on
            ylim([00, 70])
            xlim([6.5 7.5])
            % ylim([min(u1_sat) * 1.25, max(u1_sat) * 1.25])
        % end
    end

    if POSITION_FLAG
        set(gcf, 'Position',  [100, 100, 560, 280])
    end
end

%               *** weight ball ***
function [] = plot_weight_norm(result, POSITION_FLAG, CTRL_NUM)
    % plot(randn(100,1));return;
    tl = tiledlayout(2, 1);

    nnOpt = result.nnOpt;
    cstr = result.nnOpt.cstr;
    t = result.t;

    th = result.V_hist;

    font_size = 16;
    line_width = 2;
    lgd_size = 12;

    l_len = size(th, 1);
% c_list = rand(l_len, 3);
    c_list = eye(3);

    for l_idx = 1:1:l_len
        c = c_list(l_idx, :);

        plot(t, th(l_idx, :), 'color', c, 'DisplayName',"$\Vert\hat\theta_"+string(l_idx-1)+"\Vert$" ...
            , "LineWidth", line_width, "LineStyle", "-"); hold on
        plot(t, ones(size(t)) * cstr.V_max(l_idx), "color", c, 'DisplayName',"$\bar \theta_"+string(l_idx-1)+"$", ...
            "LineWidth", line_width, "LineStyle", "-."); hold on
    end
    lgd = legend;
    lgd.Orientation = 'Vertical';
    if CTRL_NUM == 4
        lgd.Location = 'southeast';
    else
        lgd.Location = 'northwest';
    end
    lgd.Interpreter = 'latex';
    lgd.NumColumns = nnOpt.l_size-1;
    lgd.FontSize = lgd_size;

    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("Weights Norm $ $", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on 

    ylim([0, max(cstr.V_max) * 1.25])

    % ylim([0, max(th, [], 'all') * 1.25])
    % ylim([0, max(th, [], 'all')+20])

    if CTRL_NUM == 4
    rectangle('Position', [6, 38, 7-6, 3], 'EdgeColor', [0.8500    0.3250    0.0980], 'LineWidth', 2, 'LineStyle', '-.'); hold on
    rectangle('Position', [2, 38, 3-2, 3], 'EdgeColor', [0.8500    0.3250    0.0980], 'LineWidth', 2, 'LineStyle', '-.'); hold on
    
    % (a) text 배치
    x_limits = xlim;
    y_limits = ylim;
    text_x_ratio = 0.23;
    text_y_ratio = 0.92;
    text_x = x_limits(1) + text_x_ratio * (x_limits(2) - x_limits(1));
    text_y = y_limits(1) + text_y_ratio * (y_limits(2) - y_limits(1));  
    text(text_x, text_y, {'(a)'}, 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
    
    % (b) text 배치
    x_limits = xlim;
    y_limits = ylim;
    text_x_ratio = 0.63;
    text_y_ratio = 0.92;
    text_x = x_limits(1) + text_x_ratio * (x_limits(2) - x_limits(1));
    text_y = y_limits(1) + text_y_ratio * (y_limits(2) - y_limits(1));  
    text(text_x, text_y, {'(b)'}, 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');


        axes("Position",[.65 .55 0.2,0.15])
        foc_ti = 6;
        foc_tf = 7;
        dt = 1e-4;
        foc_t = foc_ti/dt:foc_tf/dt;
        foc_t = uint32(foc_t);    

        plot(t(foc_t), ones(size(foc_t)) * cstr.V_max(3), "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
        plot(t(foc_t), th(3, foc_t), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

        % xlabel("$t\ [\rm s]$", "Interpreter", "latex")
        % ylabel("$\Vert\tau\Vert$", "Interpreter", "latex")
        set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
        grid on 
        xlim([foc_ti foc_tf])
        ylim([38 41])
        % pbaspect([1 1 1])
        grid on

        axes("Position",[.35 .55 0.2,0.15])
        foc_ti = 2;
        foc_tf = 3;
        dt = 1e-4;
        foc_t = foc_ti/dt:foc_tf/dt;
        foc_t = uint32(foc_t);    

        plot(t(foc_t), ones(size(foc_t)) * cstr.V_max(3), "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
        plot(t(foc_t), th(3, foc_t), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

        % xlabel("$t\ [\rm s]$", "Interpreter", "latex")
        % ylabel("$\Vert\tau\Vert$", "Interpreter", "latex")
        set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
        grid on 
        xlim([foc_ti foc_tf])
        ylim([38 41])
        % pbaspect([1 1 1])
        grid on

    end

    % if CTRL_NUM == 2
    %     axes('Position',[.7 .75 0.2,0.15])

    %     x2 = result.x_hist(3:4,:);
    %     r2 = result.rd_hist;
    %     e2 = x2 - r2;

    %     % foc_ti = 1e-4;
    %     foc_ti = 1.5;
    %     foc_tf = 3.5;
    %     dt = 1e-4;
    %     foc_t = foc_ti/dt:foc_tf/dt;
    %     foc_t = uint32(foc_t);        
        
    %     plot(t(foc_t), e2(1, foc_t), "Color", "magenta", "LineWidth", line_width, "LineStyle", "-"); hold on
    %     plot([foc_ti foc_tf], [0 0], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
    %     % xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    %     set(gca,'xtick',[])
    %     set(gca,'xticklabel',[])
    %     ylabel("${e_{2}}_{(1)}$", "Interpreter", "latex")
    %     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
    %     grid on 
    %     xlim([foc_ti foc_tf])
    %     ylim([-.5 .5])
    %     % pbaspect([1 1 1])

    %     axes('Position',[.7 .57 0.2,0.15])  
        
    %     plot(t(foc_t), e2(2, foc_t), "Color", "magenta", "LineWidth", line_width, "LineStyle", "-"); hold on
    %     plot([foc_ti foc_tf], [0 0], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
    %     xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    %     ylabel("${e_{2}}_{(2)}$", "Interpreter", "latex")
    %     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
    %     grid on 
    %     xlim([foc_ti foc_tf])
    %     % ylim([-cstr.u_ball*1.25, cstr.u_ball*1.25])
    %     ylim([-.02 .02])
    %     % pbaspect([1 1 1])
    % end

    if POSITION_FLAG
        set(gcf, 'Position',  [100, 100, 560, 280])
    end    
end

%           *** multipliers ***
function [] = plot_multipliers(ctrl, POSITION_FLAG)
    t = ctrl.t;
    L = ctrl.L_hist;
    L = L([1 2 3 8], :);

    font_size = 16;
    line_width = 1;
    lgd_size = 12;

    c_list = [eye(3), [.88,.57,.39]'];

    l_list = ["\lambda_{b_0}","\lambda_{b_1}","\lambda_{b_2}","\lambda_{u_b}"];

    for l_idx = 1:1:length(l_list)
        % c = rand(1,3);
        c = c_list(:, l_idx);

        semilogy(t, L(l_idx, :), 'color', c, 'DisplayName',"$"+l_list(l_idx)+"$" ...
        , "LineWidth", line_width, "LineStyle", "-"); hold on
    end
    grid on
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("$\lambda_j$ (log scale)", "Interpreter", "latex")
   
    lgd = legend;
    lgd.Orientation = 'horizontal';
    lgd.Location = 'southeast';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size;
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')

        % axes('Position',[0.78,0.55,0.1,0.2])
        % 
        % u1 = ctrl.u_hist(1,:);
        % u2 = ctrl.u_hist(2,:);
        % uball = ctrl.nnOpt.cstr.u_ball;
        % 
        % % foc_ti = 1e-4;
        % foc_ti = 2.62;
        % foc_tf = 2.625;
        % dt = 1e-4;
        % foc_t = foc_ti/dt:foc_tf/dt;
        % foc_t = uint32(foc_t);        
        % 
        % plot([foc_ti foc_tf], [1 1]*uball, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
        % plot(t(foc_t), sqrt(u1(foc_t).^2+u2(foc_t).^2), "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
        % xlabel("$t\ (\rm s)$", "Interpreter", "latex")
        % ylabel("$\Vert\tau\Vert$", "Interpreter", "latex")
        % set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
        % grid on 
        % xlim([foc_ti foc_tf])
        % ylim([49.7 50.3])
        % % pbaspect([1 1 1])
        % grid on
        

    if POSITION_FLAG
        set(gcf, 'Position',  [100, 100, 560, 280])
    end
end

function [] = plot_saturation(POSITION_FLAG)
    tiledlayout(1,2);
    font_size = 16;
    line_width = 2;
    lgd_size = 12;
    
    nexttile;
    p = 40:20:100;
    tt = 0:0.1:750;
    
    M = 50; m = -50;
    alp = (M+m)/2;
    mu = (M-m)/2;
    
    for p_idx = 1:1:length(p)
        tanh_t = mu * ( (tt/mu) ./ ((1+(tt/mu).^p(p_idx)).^(1/p(p_idx))) );
    
        plot(tt, tanh_t, "DisplayName", "$p = "+string(p(p_idx))+"$", "LineWidth", line_width); hold on
    end
    
    plot(tt, min(tt, M), 'HandleVisibility','off', "LineWidth", line_width, 'color', 'k', 'LineStyle', '-.'); hold on
    xlim([49, 52])
    ylim([49, 50.2])
    
    X = [0.3 0.25];
    Y = [0.6   0.75];
    ano = annotation('textarrow', X,Y, 'String', '$p \ \uparrow$', 'Interpreter', 'latex');
    ano.FontSize = font_size;
    ano.FontName = 'Times New Roman';
    ano.LineWidth = line_width;
    ano.HorizontalAlignment = 'right';
    
    xlabel("$\Vert\tau\Vert$", "Interpreter", "latex")
    ylabel("SSF$_L^U(\Vert\tau\Vert)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    
    grid on
    lgd = legend;
    lgd.Orientation = 'Vertical';
    lgd.Location = 'southeast';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 
    
    nexttile;
    tt = (-150:2:150);
    
    dhdu_F_hist = zeros(length(tt), length(tt));
    
    dhdu_11 = @(t1,t2) dfdx_11(t1, t2);
    dhdu_12 = @(t1,t2) dfdx_12(t1, t2);
    dhdu_21 = @(t1,t2) dfdx_21(t1, t2);
    dhdu_22 = @(t1,t2) dfdx_22(t1, t2);
    dhdu_F = @(t1, t2) norm([dhdu_11(t1, t2), dhdu_12(t1, t2); dhdu_21(t1, t2), dhdu_22(t1, t2)], 'fro');
    
    for t1_idx = 1:1:length(tt)
        for t2_idx = 1:1:length(tt)
            t1 = tt(t1_idx);
            t2 = tt(t2_idx);
            dhdu_F_hist(t1_idx, t2_idx) = dhdu_F(t1, t2);
        end
    end
    
    surfc(tt, tt, dhdu_F_hist, "EdgeAlpha", 0.1, "FaceAlpha", 0.8); hold on
    
    xlabel("$\tau_{(1)}$", "Interpreter", "latex")
    ylabel("$\tau_{(2)}$", "Interpreter","latex")
    zlabel("$\Vert\partial h/\partial \tau\Vert_F$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    zlim([-1 1.5])
    colorbar
    
    if POSITION_FLAG
        set(gcf, 'Position',  [100, 100, 560, 280])
    end

end

function [] = plot_total_control(ctrl_results, group)

    font_size = 16;
    line_width = 2;
    lgd_size = 12;

    % foc_ti = 1e-4;
    foc_ti = 5.5;
    foc_tf = 8.5;
    dt = 1e-4;
    foc_t = foc_ti/dt:foc_tf/dt;
    foc_t = uint32(foc_t);   

    nexttile(1)
    % ctrl_list = [ctrl1, ctrl2, ctrl3, ctrl4];
    ctrl_name = ["BSC" "DNN-BSC" "DNN-BSC-A" "CoNAC"];
    
    for c_idx = 1:1:length(group)
        ctrl_idx = group(c_idx);

        ctrl = ctrl_results.("ctrl"+string(ctrl_idx));
        t = ctrl.t;
        u1 = ctrl.u_hist(1,:);
    
        plot(t(foc_t), u1(foc_t), "DisplayName", ctrl_name(c_idx), "LineWidth", line_width); hold on
    end

    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("$\tau_{(1)}\ (\rm Nm)$", "Interpreter","latex")
    xlim([foc_ti foc_tf])
    ylim([-50 200])
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    lgd = legend;
    lgd.Orientation = 'Horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 
    
    nexttile(2)
    for c_idx = 1:1:length(group)
        ctrl_idx = group(c_idx);

        ctrl = ctrl_results.("ctrl"+string(ctrl_idx));
        t = ctrl.t;
        u2 = ctrl.u_hist(2,:);
    
        plot(t(foc_t), u2(foc_t), "DisplayName", ctrl_name(c_idx), "LineWidth", line_width); hold on
    end
    
    xlabel("$t\ (\rm s)$", "Interpreter", "latex")
    ylabel("$\tau_{(2)}\ (\rm Nm)$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    ylim([-10 10])
    xlim([foc_ti foc_tf])

        axes('Position',[.75 .75 .2 .15])
        % foc_ti = 1e-4;
        foc_ti = 6.1;
        foc_tf = 7.2;
        dt = 1e-4;
        foc_t = foc_ti/dt:foc_tf/dt;
        foc_t = uint32(foc_t);   


        for c_idx = 1:1:length(group)
            ctrl_idx = group(c_idx);
    
            ctrl = ctrl_results.("ctrl"+string(ctrl_idx));
            t = ctrl.t;
            u1 = ctrl.u_hist(1,:);
        
            plot(t(foc_t), u1(foc_t), "DisplayName", "CTRL"+string(ctrl_idx), "LineWidth", line_width); hold on
        end

        xlabel("$t\ (\rm s)$", "Interpreter", "latex")
        ylabel("$\tau_{(1)}\ (\rm Nm)$", "Interpreter","latex")
        set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
        grid on 
        xlim([foc_ti foc_tf])
        % ylim([49.5 52])
        ylim([40 55])
        % pbaspect([1 1 1])

    
    % if POSITION_FLAG
    %     set(gcf, 'Position',  [100, 100, 560, 250])
    % end
end

function [] = plot_ball(ctrl_results, group, c)
    c = c(3:4,:);

    font_size = 16;
    line_width = 2 ;
    lgd_size = 12;

    foc_ti = 5;
    foc_tf = 8;
    dt = 1e-4;
    foc_t = foc_ti/dt:foc_tf/dt;
    foc_t = uint32(foc_t);  
    
    u_ball = 50;
    u_max = (1/sqrt(2)*u_ball+u_ball) /2 ;
    ang = 0:0.01:2*pi;
    plot(u_ball*cos(ang), u_ball*sin(ang), "color", c(2,:), "LineWidth", line_width, "LineStyle", "-.",'HandleVisibility','off'); hold on
    % plot([u_max u_max], [-100 100], "color", 'black', "LineWidth", line_width, "LineStyle", "-.",'HandleVisibility','off'); hold on
    rectangle('Position', [-u_max -u_max 2*u_max 2*u_max], "EdgeColor", c(1,:), "LineWidth", line_width, "LineStyle", "-.",'HandleVisibility','off')

    group = [3 4];
    ctrl_name = [ "DNN-BSC-A", "CoNAC"];
    % ctrl_name = ["CM1", "CM2", "CM3", "CoNAC"];

    for c_idx = 1:1:length(group)
        ctrl_idx = group(c_idx);
        
        ctrl = ctrl_results.("ctrl"+string(ctrl_idx));
        u1 = ctrl.u_hist(1,:);
        u2 = ctrl.u_hist(2,:);

        plot(u1(foc_t), u2(foc_t), "color", c(c_idx,:), "LineWidth", line_width, "LineStyle", "-", "DisplayName", ctrl_name(c_idx)); hold on
    % plot(u1_sat, u2_sat, "color", 'blue', "LineWidth", 2, "LineStyle", "-"); hold on
    end
    xlabel("$\tau_{1}\ (\rm Nm)$", "Interpreter", "latex")
    ylabel("$\tau_{2}\ (\rm Nm)$", "Interpreter", "latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on 
    % xlim([-cstr.u_ball*1.5, cstr.u_ball*1.5])
    % ylim([-cstr.u_ball*1.5, cstr.u_ball*1.5])
    xlim([40 55])
    ylim([-4 4])
    % pbaspect([1 1 1])

    % DNN-BSC-A constraint
    x_limits = xlim;
    y_limits = ylim;
    text_x_ratio = 0.2;
    text_y_ratio = 0.15;
    text_x = x_limits(1) + text_x_ratio * (x_limits(2) - x_limits(1));
    text_y = y_limits(1) + text_y_ratio * (y_limits(2) - y_limits(1));  
    text(text_x, text_y, {'\leftarrow Input Bound', '     Constraint of','     DNN-BSC-A'}, 'FontName', 'Times New Roman', 'FontSize', font_size);
    
    % CoNAC constraint
    text_x_ratio = 0.68;
    text_y_ratio = 0.15;
    text_x = x_limits(1) + text_x_ratio * (x_limits(2) - x_limits(1));
    text_y = y_limits(1) + text_y_ratio * (y_limits(2) - y_limits(1));  
    text(text_x, text_y, {'\leftarrow Input Norm', '    Constraint of', '       CoNAC'}, 'FontName', 'Times New Roman', 'FontSize', font_size);

    lgd = legend;
    lgd.Orientation = 'Horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 

        axes('Position',[.7 .7 .25 .25])
        foc_ti = 1e-4;
        foc_ti = 6.1;
        foc_tf = 10;
        dt = 1e-4;
        foc_t = foc_ti/dt:foc_tf/dt;
        foc_t = uint32(foc_t);   

        plot(u_ball*cos(ang), u_ball*sin(ang), "color", c(2,:), "LineWidth", line_width, "LineStyle", "-.",'HandleVisibility','off'); hold on
        % plot([u_max u_max], [-100 100], "color", 'black', "LineWidth", line_width, "LineStyle", "-.",'HandleVisibility','off'); hold on
        rectangle('Position', [-u_max -u_max 2*u_max 2*u_max], "EdgeColor", c(1,:), "LineWidth", line_width, "LineStyle", "-.",'HandleVisibility','off')

        for c_idx = 1:1:length(group)
            ctrl_idx = group(c_idx);
    
            ctrl = ctrl_results.("ctrl"+string(ctrl_idx));
           
            u1 = ctrl.u_hist(1,:);
            u2 = ctrl.u_hist(2,:);
        
            plot(u1(foc_t), u2(foc_t), "color", c(c_idx,:), "LineWidth", line_width); hold on
        end

        xlabel("$\tau_{1}$", "Interpreter", "latex")
        ylabel("$\tau_{2}$", "Interpreter","latex")
        set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
        grid on 
        xlim([-u_ball*1.1, u_ball*1.1])
        ylim([-u_ball*1.1, u_ball*1.1])
        % xlim([foc_ti foc_tf])
        % ylim([49.5 52])
        % ylim([40 55])
        pbaspect([1 1 1])

end


%               *** Quantitive Report ***       
function [] = quantitive_report(ctrl_results, group)
    ctrl1 = ctrl_results.("ctrl"+string(group(1)));
    ctrl2 = ctrl_results.("ctrl"+string(group(2)));
    ctrl3 = ctrl_results.("ctrl"+string(group(3)));
    ctrl4 = ctrl_results.("ctrl"+string(group(4)));

    %% QUANTITIVE RESULTS
    rst1 = quantitive_calc(ctrl1);
    rst2 = quantitive_calc(ctrl2);
    rst3 = quantitive_calc(ctrl3);
    rst4 = quantitive_calc(ctrl4);

    %% PRINT RESULTS
    fprintf("      |  CTRL1      |  CTRL2      |  CTRL5      |  CTRL4\n");
    fprintf("RMSE  |  %.3f      |  %.3f      |  %.3f      |  %.3f\n", rst1.e, rst2.e, rst3.e, rst4.e);
    fprintf("Jeck  |  %.3e  |  %.3e  |  %.3e  |  %.3e\n", rst1.j, rst2.j, rst3.j, rst4.j);
    fprintf("      |  %.2f%%      |  %.2f%%      |  %.2f%%      |  %.2f%%\n", ...
        (rst1.j / rst4.j), (rst2.j / rst4.j), (rst3.j / rst4.j), 1.0);
    fprintf("Jeck E|  %.3e  |  %.3e  |  %.3e  |  %.3e\n", rst1.j_e, rst2.j_e, rst3.j_e, rst4.j_e);
    fprintf("      |  %.2f%%      |  %.2f%%      |  %.2f%%      |  %.2f%%\n", ...
        (rst1.j_e / rst4.j_e), (rst2.j_e / rst4.j_e), (rst3.j_e / rst4.j_e), 1.0);
end

function out = quantitive_calc(result)

    foc_ti = 4;
    foc_tf = 8;
    dt = 1e-4;
    foc_t = foc_ti/dt:foc_tf/dt;
    foc_t = uint32(foc_t);        
    
    x1 = result.x_hist(1,:);
    x2 = result.x_hist(2,:);

    r1 = result.r_hist(1,:);
    r2 = result.r_hist(2,:);

    e1 = x1 - r1; e2 = x2 - r2;
    e1 = e1(foc_t); e2 = e2(foc_t);
    e = e1.^2 + e2.^2;

    out.e = sqrt(mean(e)); % RMSE

    u1 = result.u_hist(1,:);
    u2 = result.u_hist(2,:);
    u1 = u1(foc_t); u2 = u2(foc_t);

    j1 = jeck_calc(u1);
    j2 = jeck_calc(u2);

    out.j = sqrt(mean(j1.^2 + j2.^2)); % Jeck

    j_e1 = jeck_calc(e1);
    j_e2 = jeck_calc(e2);

    out.j_e = sqrt(mean(j_e1.^2 + j_e2.^2));

end

function out = jeck_calc(in)

    dt = 1e-4;

    v = diff(in) / dt;
    a  = diff(v) / dt;

    out = diff(a) / dt;

end