clear

%%
SAVE_FLAG = 1;

%% DATA LOAD
ctrl1 = "sim_result/ctrl1/ctrl1_result.mat";
ctrl2 = "sim_result/ctrl2/ctrl2_result.mat";
ctrl3 = "sim_result/ctrl3/ctrl3_result.mat";
ctrl4 = "sim_result/ctrl4/ctrl4_result.mat";

ctrl1 = load(ctrl1);
ctrl2 = load(ctrl2);
ctrl3 = load(ctrl3);
ctrl4 = load(ctrl4);

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
    "fig12" %  multipliers; CTRL3vs4
];

%% FIGURE 1~4: STATE vs REFERENCE
figure(1); clf
figure1_plot(ctrl1); 
figure(2); clf
figure1_plot(ctrl2);
figure(3); clf
figure1_plot(ctrl3);
figure(4); clf
figure1_plot(ctrl4);

%% FIGURE 5~8: CONTROL INPUT
figure(5); clf
figure3_plot(ctrl1); 
figure(6); clf
figure3_plot(ctrl2);
figure(7); clf
figure3_plot(ctrl3);
figure(8); clf
figure3_plot(ctrl4);

%% FIGURE 3: WEIGHT BALL
figure(9); clf
figure4_plot(ctrl2);
figure(10); clf
figure4_plot(ctrl3);
figure(11); clf
figure4_plot(ctrl4);

%% FIGURE 4: MULTIPLIERS
figure(12); clf
figure6_plot(ctrl3, ctrl4);


%% SAVE FIGURES
if SAVE_FLAG
    for idx = 1:1:length(fig_names)
        fig_name = fig_names(idx);
        f_name = "sim_result/main_plot/" + fig_name;

        saveas(figure(idx), f_name + ".png")
        exportgraphics(figure(idx), f_name+'.eps')
    end
end

beep()


%% LOCAL FUNCTIONS
function [] = figure1_plot(result)
    tl = tiledlayout(2, 1);

    t = result.t;

    x1 = result.x_hist(1,:);
    x2 = result.x_hist(2,:);

    r1 = result.r_hist(1,:);
    r2 = result.r_hist(2,:);

    font_size = 16;
    line_width = 2;
    lgd_size = 12;

    nexttile
    plot(t, r1, "Color", "green", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$r_1$"); hold on
    plot(t, x1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$q$"); hold on
    xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    ylabel("$q_{(1)}\ [\rm rad]$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    % ylim([min(r1) * 0.75, max(r1) * 1.25])
    % ylim([-0.5, 2.5])
    ylim([-1.5, 4.5])
    lgd = legend;
    lgd.Orientation = 'Horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 

    nexttile
    plot(t, r2, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(t, x2, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
    xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    ylabel("$q_{(2)}\ [\rm rad]$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    % ylim([min(r2) * 1.25, max(r2) * 0.75])
    ylim([-2.5, .5])
    % ylim([-2.5, 2.5])

    set(gcf, 'Position',  [100, 100, 560, 350])

end

function [] = figure3_plot(result)
    % plot(randn(100,1));return;
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
    plot(t, u1_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "Saturated $\tau$"); hold on
    % plot(t(actset1), u1(actset1), "Color", "red", "LineWidth", line_width, "LineStyle", "."); hold on
    % scatter(t(actset1), u1(actset1), "Color", "red", "Marker","."); hold on
    xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    ylabel("$\tau_{(1)}\ [\rm Nm]$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    ylim([-50 * 1.25, 50 * 1.25])
    % ylim([min(u1_sat) * 1.25, max(u1_sat) * 1.25])
    lgd = legend;
    lgd.Orientation = 'Horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size;    
    
    nexttile
    plot(t, u2, "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(t, u2_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
    xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    ylabel("$\tau_{(2)}\ [\rm Nm]$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on
    ylim([-50 * 1.25, 50 * 1.25])
    % ylim([min(u2_sat) * 1.25, max(u2_sat) * 1.25])
    
        axes('Position',[.75 .75 .2 .2])
        ang = 0:0.01:2*pi;
        plot(cstr.u_ball*cos(ang), cstr.u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
        plot(u1, u2, "color", 'red', "LineWidth", 2, "LineStyle", "-"); hold on
        plot(u1_sat, u2_sat, "color", 'blue', "LineWidth", 2, "LineStyle", "-"); hold on
        xlabel("$\tau_{(1)}$", "Interpreter", "latex")
        ylabel("$\tau_{(2)}$", "Interpreter", "latex")
        set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
        grid on 
        xlim([-cstr.u_ball*1.25, cstr.u_ball*1.25])
        ylim([-cstr.u_ball*1.25, cstr.u_ball*1.25])
        pbaspect([1 1 1])
   
    set(gcf, 'Position',  [100, 100, 560, 350])

end

function [] = figure4_plot(result)
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
    lgd.Location = 'northwest';
    lgd.Interpreter = 'latex';
    lgd.NumColumns = nnOpt.l_size-1;
    lgd.FontSize = lgd_size;

    xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    ylabel("Weights Norm $ $", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    grid on 
    ylim([0, max(cstr.V_max) * 1.25])
    % ylim([0, max(th, [], 'all') * 1.25])
    % ylim([0, max(th, [], 'all')+20])

    set(gcf, 'Position',  [100, 100, 560, 350])
    
end

function [] = figure6_plot(ctrl3, ctrl4)
    % plot(randn(100,1));return;
    tl = tiledlayout(2, 1);

    t_ctrl3 = ctrl3.t;
    L_ctrl3 = ctrl3.L_hist;
    L_ctrl3 = L_ctrl3([1 2 3 8], :);

    t_ctrl4 = ctrl4.t;
    L_ctrl4 = ctrl4.L_hist;
    L_ctrl4 = L_ctrl4([1 2 3 8], :);

    font_size = 16;
    line_width = 1;
    lgd_size = 12;

    c_list = [eye(3), [.88,.57,.39]'];

    l_list = ["\lambda_{b,0}","\lambda_{b,1}","\lambda_{b,2}","\lambda_{u_b}"];

    nexttile
    for l_idx = 1:1:size(L_ctrl3, 1)
        % c = rand(1,3);
        c = c_list(:, l_idx);
    
        semilogy(t_ctrl3, L_ctrl3(l_idx, :), 'color', c, 'DisplayName',"$"+l_list(l_idx)+"$" ...
        , "LineWidth", line_width, "LineStyle", "-"); hold on
    end
    grid on
    xlabel("$t\ [\rm s]$", "Interpreter", "latex")
    ylabel("$\lambda_i$ (log scale)", "Interpreter", "latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
    lgd = legend;
    lgd.Orientation = 'horizontal';
    lgd.Location = 'northoutside';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size;

    nexttile
    for l_idx = 1:1:size(L_ctrl4, 1)
        % c = rand(1,3);
        c = c_list(:, l_idx);

        semilogy(t_ctrl4, L_ctrl4(l_idx, :), 'color', c, 'DisplayName',"$\lambda_"+string(l_idx)+"$" ...
        , "LineWidth", line_width, "LineStyle", "-"); hold on
    end
    grid on
    xlabel("Time $[\rm s]$", "Interpreter", "latex")
    ylabel("$\lambda_i$ (log scale)", "Interpreter", "latex")
   
    set(gcf, 'Position',  [100, 100, 560, 350])

end


