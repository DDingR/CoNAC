
%% FIG.1: STATE vs REFERENCE
% ********************************************************
figure(1); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

p1 = plot(t, r1, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
p2 = plot(t, x1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("State 1 $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(r1) * 0.75, max(r1) * 1.25])
% ylim([-0.5, 2.5])
legend([p1, p2], ["$r_1$","$x_1$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

figure(2); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

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

figure(3); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

p1 = plot(t, rd1, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
p2 = plot(t, xd1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("State Dot 1 $[\rm rad/s]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(rd1) * 1.25, max(rd1) * 1.25])
legend([p1, p2], ["$\dot r_1$","$\dot x_1$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

figure(4); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

p1 = plot(t, rd2, "Color", "green", "LineWidth", line_width, "LineStyle", "-"); hold on
p2 = plot(t, xd2, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("State Dot 2 $[\rm rad/s]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(rd2) * 1.25, max(rd2) * 1.25])
legend([p1, p2], ["$\dot r_2$","$\dot x_2$"], "Interpreter","latex", "FontSize", lgd_size, "FontWeight", "bold", "Location", "northwest")

%% FIG.3: CONTROL INPUT 
% ********************************************************
figure(5); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

plot(t, ones(size(t)) * cstr.uMax1, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
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

figure(6); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

plot(t, ones(size(t)) * cstr.uMax2, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(t, u2, "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(t, u2_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\tau_2\ [\rm Nm]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
ylim([min(u2_sat) * 1.25, max(u2_sat) * 1.25])

%% FIG.4: WEIGHTS BALL CONTSRAINT
% ********************************************************
figure(7); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

l_len = size(th, 1);
% c_list = rand(l_len, 3);
c_list = eye(3);

for l_idx = flip(1:1:l_len)
    c = c_list(l_idx, :);

    plot(t, th(l_idx, :), 'color', c, 'DisplayName',"$\Vert\hat\theta_"+string(l_idx-1)+"\Vert$" ...
        , "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(t, ones(size(t)) * cstr.th_max(l_idx), "color", c, 'DisplayName',"$\bar \theta_"+string(l_idx-1)+"$", ...
         "LineWidth", line_width, "LineStyle", "-."); hold on
end
lgd = legend;
lgd.Orientation = 'Vertical';
lgd.Location = 'northwest';
lgd.Interpreter = 'latex';
lgd.NumColumns = opt.l_size-1;
lgd.FontSize = lgd_size;

xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("Weights Norm $ $", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on 
ylim([0, max(cstr.th_max) * 1.25])
% ylim([0, max(th, [], 'all') * 1.25])
% ylim([0, max(th, [], 'all')+20])

%% FIG.5: CONTROL TOP VIEW
% ********************************************************
figure(8); clf
hF = gcf; 
% hF.Position(3:4) = [fig_width, fig_height];

ang = 0:0.01:2*pi;

plot(cstr.u_ball*cos(ang), cstr.u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([-100, 100], [1, 1] * cstr.uMax2, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([-100, 100], [-1, -1] * cstr.uMax2, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([1, 1] * cstr.uMax1, [-100, 100], "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([-1, -1] * cstr.uMax1, [-100, 100], "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
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
figure(9); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];
tl = tiledlayout(2, 1);

nexttile(1);
p2 = plot(t, z1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("Zeta 1 $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on

nexttile(2);
p2 = plot(t, z2, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("Zeta 2 $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on

%% FIG.10: CONTROL INPUT NORM
figure(10); clf
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

u_norm = sqrt(u1.^2 + u2.^2);

plot(t, ones(size(t)) * cstr.u_ball, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(t, u_norm, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\Vert \tau \Vert[\rm Nm]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on