
%% FIGURE SETTING
% ********************************************************
fig_names = [
    "fig1 state vs ref"
    "fig2 dot state vs ref"
    "fig3 control input"
    "fig4 weight ball"
    "fig5 input ball"
    "fig6 multipliers & dot Lagrangian"        
];

fig_len = length(fig_names);

%% FIG.1: STATE vs REFERENCE
% ********************************************************
figure(1); clf
tl = tiledlayout(2, 1);

nexttile
p1 = plot(t, r1, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
p2 = plot(t, x1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("State 1 $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(r1) * 0.75, max(r1) * 1.25])
% ylim([-0.5, 2.5])
legend([p1, p2], ["$r_1$","$x_1$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

nexttile
p1 = plot(t, r2, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
p2 = plot(t, x2, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("State 2 $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(r2) * 1.25, max(r2) * 0.75])
% ylim([-2.5, .5])
legend([p1, p2], ["$r_2$","$x_2$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

%% FIG.2: DOT STATE vs REFERENCE
% ********************************************************

figure(2); clf
tl = tiledlayout(2, 1);

nexttile
p1 = plot(t, rd1, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
p2 = plot(t, xd1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("State Dot 1 $[\rm rad/s]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
ylim([min(rd1) * 1.25, max(rd1) * 1.25])
legend([p1, p2], ["$\dot r_1$","$\dot x_1$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

nexttile
p1 = plot(t, rd2, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
p2 = plot(t, xd2, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("State Dot 2 $[\rm rad/s]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
ylim([min(rd2) * 1.25, max(rd2) * 1.25])
legend([p1, p2], ["$\dot r_2$","$\dot x_2$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

%% FIG.3: CONTROL INPUT 
% ********************************************************
figure(3); clf
tl = tiledlayout(2, 1);

% actset1 = find(u1 > cstr.U_max(1) | u1 < cstr.U_min(1));
% actset2 = find(u2 > cstr.U_max(2) | u2 < cstr.U_min(2));

nexttile
plot(t, u1, "Color", "red", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\tau$"); hold on
plot(t, u1_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "Saturated $\tau$"); hold on
% plot(t(actset1), u1(actset1), "Color", "red", "LineWidth", line_width, "LineStyle", "."); hold on
% scatter(t(actset1), u1(actset1), "Color", "red", "Marker","."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\tau_1\ [\rm Nm]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
ylim([min(u1_sat) * 1.25, max(u1_sat) * 1.25])
lgd = legend;
lgd.Orientation = 'Horizontal';
lgd.Location = 'northoutside';
lgd.Interpreter = 'latex';
lgd.FontSize = lgd_size;    

nexttile
plot(t, u2, "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(t, u2_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\tau_2\ [\rm Nm]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
ylim([min(u2_sat) * 1.25, max(u2_sat) * 1.25])

    axes('Position',[.75 .75 .2 .2])
    ang = 0:0.01:2*pi;
    plot(cstr.u_ball*cos(ang), cstr.u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
    plot(u1, u2, "color", 'red', "LineWidth", 2, "LineStyle", "-"); hold on
    plot(u1_sat, u2_sat, "color", 'blue', "LineWidth", 2, "LineStyle", "-"); hold on
    xlabel("$\tau_1$", "Interpreter", "latex")
    ylabel("$\tau_2$", "Interpreter", "latex")
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')
    grid on 
    xlim([-cstr.u_ball*1.25, cstr.u_ball*1.25])
    ylim([-cstr.u_ball*1.25, cstr.u_ball*1.25])
    pbaspect([1 1 1])

%% FIG.4: WEIGHTS BALL CONTSRAINT
% ********************************************************
figure(4); clf

l_len = size(th, 1);
% c_list = rand(l_len, 3);
c_list = eye(3);

for l_idx = flip(1:1:l_len)
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

xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("Weights Norm $ $", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on 
ylim([0, max(cstr.V_max) * 1.25])
% ylim([0, max(th, [], 'all') * 1.25])
% ylim([0, max(th, [], 'all')+20])

%% FIG.5: CONTROL BALL CONSTRAINT
% ********************************************************
figure(5); clf

ang = 0:0.01:2*pi;

plot(cstr.u_ball*cos(ang), cstr.u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
p1 = plot(u1, u2, "color", 'red', "LineWidth", 2, "LineStyle", "-"); hold on
p2 = plot(u1_sat, u2_sat, "color", 'blue', "LineWidth", 2, "LineStyle", "-"); hold on
xlabel("$\tau_1\ [\rm Nm]$", "Interpreter", "latex")
ylabel("$\tau_2\ [\rm Nm]$", "Interpreter", "latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on 
xlim([-cstr.u_ball*1.25, cstr.u_ball*1.25])
ylim([-cstr.u_ball*1.25, cstr.u_ball*1.25])
pbaspect([1 1 1])
legend([p1, p2], ["$\tau$", "Saturated $\tau$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

%% FIG.6: MULTIPLIERS & DOT LAGRANGIAN
% ********************************************************
figure(6); clf
tl = tiledlayout(2, 1);

nexttile
for l_idx = 1:1:size(L, 1)
    c = rand(1,3);

    % plot(t, L(l_idx, :), 'color', c, 'DisplayName',"$\lambda_"+string(l_idx)+"$" ...
    % , "LineWidth", line_width, "LineStyle", "-"); hold on
    semilogy(t, L(l_idx, :), 'color', c, 'DisplayName',"$\lambda_"+string(l_idx)+"$" ...
    , "LineWidth", line_width, "LineStyle", "-"); hold on
end
grid on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\lambda_i$ (log scale)", "Interpreter", "latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
lgd = legend;
lgd.Orientation = 'horizontal';
lgd.Location = 'northoutside';
lgd.Interpreter = 'latex';
lgd.FontSize = lgd_size;

nexttile
semilogy(t, dot_L, 'color', 'blue', 'DisplayName', "$\dot L$" ...
        , "LineWidth", line_width, "LineStyle", "-"); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\dot L$", "Interpreter", "latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
lgd = legend;
lgd.Interpreter = 'latex';
lgd.FontSize = lgd_size;
lgd.Location = 'northwest';