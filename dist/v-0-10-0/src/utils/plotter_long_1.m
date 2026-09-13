

figure(1); clf
tl = tiledlayout(6, 2);
% hF = gcf; 
% hF.Position(3:4) = [fig_width, fig_height];

start_t = 7.5;

%% [1, 1:2] State vs Reference (Joint 1)
nexttile([1, 2])
plot(t, x1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$x_1$"); hold on
plot(t, r1, "Color", "red", "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$r_1$"); hold on
maxVal = max([x1; r1], [], 'all');
minVal = min([x1; r1], [], 'all');
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel(" $q_1$ $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(r1) * 0.75, max(r1) * 1.25])
% ylim([-0.5, 2.5])
% lgd = legend;
% lgd.Interpreter = 'latex';
% lgd.FontSize = lgd_size;
% lgd.FontWeight = "bold";
% lgd.Location = "northwest";


xlim([start_t, t(end)])
len = maxVal - minVal;
if len ~= 0; ratio = 0.25; ylim([minVal-len*ratio, maxVal+len*ratio]); end

%% [2, 1:2] State vs Reference (Joint 2)
nexttile([1, 2])
plot(t, x2, "Color", "blue", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(t, r2, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
maxVal = max([x2; r2], [], 'all');
minVal = min([x2; r2], [], 'all');
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel(" $q_2$ $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(r2) * 1.25, max(r2) * 0.75])
% ylim([-2.5, .5])
% lgd = legend;
% lgd.Interpreter = 'latex';
% lgd.FontSize = lgd_size;
% lgd.FontWeight = "bold";
% lgd.Location = "northwest";

xlim([start_t, t(end)])
len = maxVal - minVal;
if len ~= 0; ratio = 0.25; ylim([minVal-len*ratio, maxVal+len*ratio]); end

%% [3, 1:2] Control Input (Joint 1)
nexttile([1, 2])
plot(t, ones(size(t)) * cstr.uMax1, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(t, u1, "Color", "red", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\tau$"); hold on
plot(t, u1_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "Saturated $\tau$"); hold on
maxVal = max([u1; u1_sat], [], 'all');
minVal = min([u1; u1_sat], [], 'all');
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\tau_1\ [\rm Nm]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% lgd = legend;
% lgd.Orientation = 'Horizontal';
% lgd.Location = 'northoutside';
% lgd.Interpreter = 'latex';
% lgd.FontSize = lgd_size;    

xlim([start_t, t(end)])
len = maxVal - minVal;
if len ~= 0; ratio = 0.25; ylim([minVal-len*ratio, maxVal+len*ratio]); end

%% [4, 1:2] Control Input (Joint 2)
nexttile([1, 2])
plot(t, ones(size(t)) * cstr.uMax2, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(t, u2, "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(t, u2_sat, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
maxVal = max([u2; u2_sat], [], 'all');
minVal = min([u2; u2_sat], [], 'all');
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\tau_2\ [\rm Nm]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on

xlim([start_t, t(end)])
len = maxVal - minVal;
if len ~= 0; ratio = 0.25; ylim([minVal-len*ratio, maxVal+len*ratio]); end

%% [5, 1:2] Control Input Norm
nexttile([1, 2])
u_norm = sqrt(u1.^2 + u2.^2);
plot(t, ones(size(t)) * cstr.u_ball, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(t, u_norm, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel("$\Vert \tau \Vert[\rm Nm]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on

xlim([start_t, t(end)])
len = max(u_norm) - min(u_norm);
if len ~= 0; ratio = 0.25; ylim([min(u_norm)-len*ratio, max(u_norm)+len*ratio]); end

%% FIG.4: WEIGHTS BALL CONTSRAINT
nexttile([1, 2])

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
xlim([start_t, t(end)])
ylim([0, max(cstr.th_max) * 1.25])
% ylim([0, max(th, [], 'all') * 1.25])
% ylim([0, max(th, [], 'all')+20])

%% Aux Sys
figure(2);clf
plot(t, zeta_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$x_1$"); hold on
plot(t, zeta_hist(2,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$r_1$"); hold on
maxVal = max([x1; r1], [], 'all');
minVal = min([x1; r1], [], 'all');
xlabel("Time $[\rm s]$", "Interpreter", "latex")
ylabel(" $q_1$ $[\rm rad]$", "Interpreter","latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on
% ylim([min(r1) * 0.75, max(r1) * 1.25])
% ylim([-0.5, 2.5])
% lgd = legend;
% lgd.Interpreter = 'latex';
% lgd.FontSize = lgd_size;
% lgd.FontWeight = "bold";
% lgd.Location = "northwest";


xlim([start_t, t(end)])
len = maxVal - minVal;
if len ~= 0; ratio = 0.25; ylim([minVal-len*ratio, maxVal+len*ratio]); end
