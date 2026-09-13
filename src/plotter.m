

%%
% SIM = true;
SIM = false;
SAVE_FLAG = 1;

uMax2 = 3.8;
u_ball = 11;
uMax1 = sqrt(u_ball^2 - uMax2^2);
th_max = [6,6,6]*10;

warmup_time = 11+3;
ep_time = 2*4;

global_start_t = warmup_time;
% global_end_t = T;
global_end_t = warmup_time + 2*ep_time;

start_t = 26.6; 
end_t   = 27.3;

% start_t = 32.4; 
% end_t   = 34.9;


%% VARIABLE MAPPING
if SIM
    load("RESULT.mat")
    T = dataSet{1}.T;
    t = dataSet{1}.t;
    obs_t = 1:length(t);

    warmup_start_t = find(t >= warmup_time, 1);
    epi_end_t = find(t >= warmup_time + 2*ep_time, 1);
    epi_idx = find(t >= warmup_time & t <= warmup_time + 2*ep_time);
    ctrl_obs_idx = find(t >= start_t & t <= end_t);
    ctrl_obs_idx2 = find(t >= start_t-0.5 & t <= end_t+0.5); %little longer
else
    dataSet = cell(4,1);

    data_path_4 = "20260728_161911"; % CONAC high - 1
    data_path_3 = "20260728_162454"; % CONAC low - 2
    data_path_2 = "20260728_151903"; % AUX - 3
    data_path_1 = "20260728_141636"; % NAC - 4

    dataSet{1} = loadFromMeas("meas_result/"+data_path_1+"/"+data_path_1+".mat", 4);
    dataSet{2} = loadFromMeas("meas_result/"+data_path_2+"/"+data_path_2+".mat", 3);
    dataSet{3} = loadFromMeas("meas_result/"+data_path_3+"/"+data_path_3+".mat", 2);
    dataSet{4} = loadFromMeas("meas_result/"+data_path_4+"/"+data_path_4+".mat", 1);

    [t_end_max, t_max_idx] = max([length(dataSet{1}.t), length(dataSet{2}.t), length(dataSet{3}.t), length(dataSet{4}.t)]);
    % fill nan values


    t = dataSet{t_max_idx}.t;
    T = dataSet{t_max_idx}.T;
    obs_t = 1:length(t);
    warmup_start_t = find(t >= warmup_time, 1);
    epi_end_t = find(t >= warmup_time + 2*ep_time, 1);
    epi_idx = find(t >= warmup_time & t <= warmup_time + 2*ep_time);
    ctrl_obs_idx = find(t >= start_t & t <= end_t);
    ctrl_obs_idx2 = find(t >= start_t-.5 & t <= end_t+.5); %little longer
    % ctrl_obs_idx2 = ctrl_obs_idx
end


CTRL_NUM = 4;
color_list = [ ...
    "#808080";
    "magenta";
    "cyan";
    "blue";
];
tex_name_list = {
    "NAC";
    "AUX";
    "CONAClow";
    "CONAChigh";
};
name_list = {
    "(C$_4$)";
    "(C$_3$)";
    "(C$_2$)";
    "(C$_1$)";
};

%%
font_size = 12;
ax_font_size = 12;
line_width = 1;
lgd_size = 16;

fig_height = 200;
fig_width = 1000;
ax_height = fig_height * .9;
ax_width = fig_width * .9;

fig_unit = 'pixels';

norm_func = @(x) sqrt(sum(x.^2, 1));

%% FIG. 1-4
figW = 17;   % cm, one-column width 정도
figH = 5;
axSize = [2.0 2 14.3 5.1];

ax_list = {};
% maxVals = ones(5,1)*-inf; minVals = ones(5,1)*inf;
maxVals = [
    50
    60
    rad2deg(3)
    rad2deg(3)
    10.5
    4
    10.5
];
minVals = [
    -60
    -90
    rad2deg(-3)
    rad2deg(-3)
    6.5
    0
    % -10.4
    % -3.5
    6.5
];

for f_idx = [1:1:7, 13, 14]
    fig = figure(f_idx); clf;
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, 'on');
    grid(ax, 'minor');
    box(ax, 'on');

    set(fig, 'Units', 'centimeters');
    fig.Position(3:4) = [figW figH];
    % set(ax, 'Units', 'centimeters');
    % set(ax, 'Position', axSize);

    set(ax, 'FontName', 'Times New Roman');
    set(ax, 'FontSize', ax_font_size);
    set(ax, 'LineWidth', 1.1);
    set(ax, 'TickLabelInterpreter', 'latex');

    % set(fig, 'PaperUnits', 'centimeters');
    % set(fig, 'PaperSize', [figW figH]);
    % set(fig, 'PaperPosition', [0 0 figW figH]);
    % print(fig, 'fig_error.pdf', '-dpdf', '-painters');

    ax_list{f_idx} = ax;
end

% ============================
% before all data


% ============================
% all data
for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    CTRL_INFO = data.CTRL_INFO;

    x1_hist = rad2deg(data.x1_hist);
    x2_hist = rad2deg(data.x2_hist);
    xd1_hist = rad2deg(data.xd1_hist);
    xd2_hist = rad2deg(data.xd2_hist);
    u_hist = data.u_hist;
    uSat_hist = data.uSat_hist;
    color = color_list(ctrl_idx);

    % fill nan values for same length
    if length(t) > length(x1_hist)
        x1_hist = [x1_hist, nan(2, length(t)-length(x1_hist))];
        x2_hist = [x2_hist, nan(2, length(t)-length(x2_hist))];
        u_hist = [u_hist, nan(2, length(t)-length(u_hist))];
        uSat_hist = [uSat_hist, nan(2, length(t)-length(uSat_hist))];
    end
    
    plot(ax_list{1}, t, x1_hist(1,:),       "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_list{2}, t, x1_hist(2,:),       "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_list{3}, t, x2_hist(1,:),       "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_list{4}, t, x2_hist(2,:),       "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_list{5}, t, u_hist(1,:),        "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_list{5}, t, uSat_hist(1,:),     "Color", color, "LineWidth", line_width, "LineStyle", "-."); 
    plot(ax_list{6}, t, u_hist(2,:),        "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_list{6}, t, uSat_hist(2,:),     "Color", color, "LineWidth", line_width, "LineStyle", "-.");
    plot(ax_list{7}, t, norm_func(u_hist),  "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    
    plot(ax_list{13}, t, x1_hist(1,:),       "Color", color, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_list{14}, t, x1_hist(2,:),       "Color", color, "LineWidth", line_width, "LineStyle", "-");  
end

% ============================
% after all data
plot(ax_list{1}, t, rad2deg(dataSet{t_max_idx}.xd1_hist(1,:)), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); 
plot(ax_list{2}, t, rad2deg(dataSet{t_max_idx}.xd1_hist(2,:)), "Color", "red", "LineWidth", line_width, "LineStyle", "--");
plot(ax_list{3}, t, rad2deg(dataSet{t_max_idx}.xd2_hist(1,:)), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); 
plot(ax_list{4}, t, rad2deg(dataSet{t_max_idx}.xd2_hist(2,:)), "Color", "red", "LineWidth", line_width, "LineStyle", "--");

plot(ax_list{5}, [t(1) t(end)], [+1 +1]*uMax1, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax_list{5}, [t(1) t(end)], [-1 -1]*uMax1, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax_list{6}, [t(1) t(end)], [+1 +1]*uMax2, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax_list{6}, [t(1) t(end)], [-1 -1]*uMax2, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax_list{7}, [t(1) t(end)], [+1 +1]*u_ball, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax_list{7}, [t(1) t(end)], [-1 -1]*u_ball, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on

yticks(ax_list{2}, [-90 0 90]);  % 실제 0, 90 위치에 눈금

plot(ax_list{13}, t, rad2deg(dataSet{t_max_idx}.xd1_hist(1,:)), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); 
plot(ax_list{14}, t, rad2deg(dataSet{t_max_idx}.xd1_hist(2,:)), "Color", "red", "LineWidth", line_width, "LineStyle", "--");


% 
zoom_color = "#008602"; % orange
for f_idx = 1:1:7
    % plot(ax_list{f_idx}, [warmup_time, warmup_time], [-5e2 5e2], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
    plot(ax_list{f_idx}, [warmup_time+ep_time, warmup_time+ep_time], [-5e2 5e2], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
    % plot(ax_list{f_idx}, [warmup_time+2*ep_time, warmup_time+2*ep_time], [-5e2 5e2], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
    
    plot(ax_list{f_idx}, [start_t, start_t], [-5e2 5e2], "Color", zoom_color, "LineWidth", line_width, "LineStyle", "--"); hold on
    plot(ax_list{f_idx}, [end_t, end_t], [-5e2 5e2], "Color", zoom_color, "LineWidth", line_width, "LineStyle", "--"); hold on

    if f_idx <= 4
        text(ax_list{f_idx}, 0.02, 0.9, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')
        text(ax_list{f_idx}, 0.52, 0.9, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')
        text(ax_list{f_idx}, 0.845, 0.9, "see, Fig. 9", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized', 'Color', zoom_color)
    else
        text(ax_list{f_idx}, 0.02, 0.1, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')
        text(ax_list{f_idx}, 0.52, 0.1, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')
        text(ax_list{f_idx}, 0.845, 0.1, "see, Fig. 9", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized', 'Color', zoom_color)
    end

    len = maxVals(f_idx)-minVals(f_idx); ratio = .3;
    if len~=0
        ax_list{f_idx}.YLim = [minVals(f_idx)-len*ratio maxVals(f_idx)+len*ratio];
    end
    ax_list{f_idx}.XLim = [global_start_t global_end_t];
end

% label
ax_list{1}.XLabel.String = 'Time / s';
ax_list{1}.YLabel.String = '$q_1$ / deg';
ax_list{2}.XLabel.String = 'Time / s';
ax_list{2}.YLabel.String = '$q_2$ / deg';
ax_list{3}.XLabel.String = 'Time / s';
ax_list{3}.YLabel.String = '$\dot{q}_1$ / deg/s';
ax_list{4}.XLabel.String = 'Time / s';
ax_list{4}.YLabel.String = '$\dot{q}_2$ / deg/s';
ax_list{5}.XLabel.String = 'Time / s';
ax_list{5}.YLabel.String = '$\tau_1$ / Nm';
ax_list{6}.XLabel.String = 'Time / s';
ax_list{6}.YLabel.String = '$\tau_2$ / Nm';
ax_list{7}.XLabel.String = 'Time / s';
ax_list{7}.YLabel.String = '$\Vert\mbox{\boldmath $\tau$}\Vert$ / Nm';

ax_list{13}.XLabel.String = 'Time / s';
ax_list{13}.YLabel.String = '$q_1$ / deg';
ax_list{14}.XLabel.String = 'Time / s';
ax_list{14}.YLabel.String = '$q_2$ / deg';

for ax_idx = [1:1:7, 13, 14]
    ax_list{ax_idx}.XLabel.Interpreter = 'latex';
    ax_list{ax_idx}.YLabel.Interpreter = 'latex';
end

%% ============================
% Top view of the control input
% ============================
fig = figure(8); clf;
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
grid(ax, 'minor');
box(ax, 'on');

set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [figW figH*1.5];
set(ax, 'FontName', 'Times New Roman');
set(ax, 'FontSize', ax_font_size);
set(ax, 'LineWidth', 1.1);
set(ax, 'TickLabelInterpreter', 'latex');

p = ax.Position;
axInset = axes( ...
    'Parent', fig, ...
    'Units', 'normalized', ...
    'Position', ...
    [p(1) + 0.50*p(3), ...
    p(2) + 0.50*p(4), ...
    0.45*p(3), ...
    0.43*p(4)]);

hold(axInset, 'on');
grid(axInset, 'on');
grid(axInset, 'minor');
box(axInset, 'on');

set(axInset, 'FontName', 'Times New Roman');
set(axInset, 'FontSize', ax_font_size);
set(axInset, 'LineWidth', 1.1);
set(axInset, 'TickLabelInterpreter', 'latex');

maxminX = [-inf inf]; maxminY = [-inf inf];
% for ctrl_idx = 1:1:length(dataSet)
for ctrl_idx = [2, 3, 4, 1]
    data = dataSet{ctrl_idx};
    CTRL_INFO = data.CTRL_INFO;

    u_hist = data.u_hist;
    u_hist = [u_hist, nan(2, length(t)-length(u_hist))];
    u_hist = u_hist(:, ctrl_obs_idx);
    color = color_list(ctrl_idx);

    if ctrl_idx == 1
        plot(axInset, u_hist(1,:), u_hist(2,:), "Color", color, "LineWidth", line_width, "LineStyle", "-");

        maxValX = [min(u_hist(1,:)) max(u_hist(1,:))];
        maxValY = [min(u_hist(2,:)) max(u_hist(2,:))];
        lenX = maxValX(2)-maxValX(1); lenY = maxValY(2)-maxValY(1);
        ratio = 0.1;
        axInset.XLim = [maxValX(1)-lenX*ratio maxValX(2)+lenX*ratio];
        axInset.YLim = [maxValY(1)-lenY*ratio maxValY(2)+lenY*ratio];


    else
        plot(ax, u_hist(1,:), u_hist(2,:), "Color", color, "LineWidth", line_width, "LineStyle", "-");

        maxminX = [min(u_hist(1,:)) max(u_hist(1,:))];
        maxminY = [min(u_hist(2,:)) max(u_hist(2,:))];
    end
end
plot(ax, u_ball*cos(0:0.01:2*pi), u_ball*sin(0:0.01:2*pi), "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax, [-1000 1000], [+1 +1]*uMax2, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax, [-1000 1000], [-1 -1]*uMax2, "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
plot(ax, [1 1]*uMax1, [-10000, 10000], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
ax.XLabel.String = '$\tau_1$ / Nm'; ax.XLabel.Interpreter = 'latex';
ax.YLabel.String = '$\tau_2$ / Nm'; ax.YLabel.Interpreter = 'latex';
maxminX = [-u_ball, u_ball];
maxminY = [-u_ball, u_ball];
ax.XLim = [9.5 14.5];
ax.XLim = [10 13];
ax.YLim = [0 4.2];

axInset.XLabel.String = '$\tau_1$ / Nm'; axInset.XLabel.Interpreter = 'latex';
axInset.YLabel.String = '$\tau_2$ / Nm'; axInset.YLabel.Interpreter = 'latex';

%% ============================
% [CONAC] Multipliers
% ============================

fig = figure(9); clf;
ax = axes(fig);

set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [figW figH];

C1_u_ball_lbd = dataSet{4}.lbd_hist(4,:);
C2_u_ball_lbd = dataSet{3}.lbd_hist(4,:);
C1_u2_max_lbd = dataSet{4}.lbd_hist(7,:);
C2_u2_max_lbd = dataSet{3}.lbd_hist(7,:);

semilogy(ax, dataSet{4}.t, C1_u_ball_lbd, "Color", color_list(4), "LineWidth", line_width, "LineStyle", "-", "DisplayName", '$\lambda_{\overline{\tau}}$'); hold on
semilogy(ax, dataSet{3}.t, C2_u_ball_lbd, "Color", color_list(3), "LineWidth", line_width, "LineStyle", "-", "DisplayName", '$\lambda_{\overline{\tau}}$'); hold on
% semilogy(ax, dataSet{4}.t, C1_u2_max_lbd, "Color", color_list(4), "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\lambda_{\overline{\tau}_2}$"); hold on
% semilogy(ax, dataSet{3}.t, C2_u2_max_lbd, "Color", color_list(3), "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\lambda_{\overline{\tau}_2}$"); hold on

hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
grid(ax, 'minor');

% lgd = legend(ax);
% lgd.Location = 'southeast';
% lgd.Interpreter = 'latex';
% lgd.FontSize = 10;
% lgd.NumColumns = 2;
set(ax, 'FontName', 'Times New Roman');
set(ax, 'FontSize', ax_font_size);
set(ax, 'LineWidth', 1.1);
set(ax, 'TickLabelInterpreter', 'latex');
ax.XLabel.String = 'Time / s';
ax.YLabel.String = '$\lambda_{\overline{\tau}}$';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';

ax.XLim = [start_t end_t];
ax.YLim = [0.000078389151053,7.12744126403738];

%% ============================
%   Weight Norms
% ============================
fig = figure(10); clf;
ax = axes(fig);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [figW figH];
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
grid(ax, 'minor');
set(ax, 'FontName', 'Times New Roman');
set(ax, 'FontSize', ax_font_size);
set(ax, 'LineWidth', 1.1);
set(ax, 'TickLabelInterpreter', 'latex');

for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    CTRL_INFO = data.CTRL_INFO;
    color = color_list(ctrl_idx);
    
    t = data.t;

    th_hist = data.th_hist;
    th0 = th_hist(1,:);
    th1 = th_hist(2,:);
    th2 = th_hist(3,:);

    plot(ax, t, th0, "Color", color, "LineWidth", line_width, "LineStyle", "-", "HandleVisibility", "off"); hold on
    plot(ax, t, th1, "Color", color, "LineWidth", line_width, "LineStyle", "-.", "HandleVisibility", "off"); hold on
    plot(ax, t, th2, "Color", color, "LineWidth", line_width, "LineStyle", "--", "HandleVisibility", "off"); hold on
    
end
% dummy for legend (not plotted)
plot(ax, NaN, NaN, "Color", 'k', "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\Vert\widehat{\mbox{\boldmath $\theta$}}_0\Vert$"); hold on
plot(ax, NaN, NaN, "Color", 'k', "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\Vert\widehat{\mbox{\boldmath $\theta$}}_1\Vert$"); hold on
plot(ax, NaN, NaN, "Color", 'k', "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\Vert\widehat{\mbox{\boldmath $\theta$}}_2\Vert$"); hold on
plot(ax, [0 T], [+1 +1]*th_max(1), "Color", "black", "LineWidth", line_width, "LineStyle", "-.", "HandleVisibility", "off"); hold on

lgd = legend(ax);
lgd.Location = 'northwest';
lgd.Interpreter = 'latex';
lgd.FontSize = 10;
lgd.NumColumns = 3;

ax.XLim = [global_start_t global_end_t];
ax.YLim = [0 14];
ax.XLabel.String = 'Time / s';
ax.YLabel.String = '$\Vert\widehat{\mbox{\boldmath $\theta$}}_i\Vert$';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';

%% ============================
% [AUX] Zeta
% ============================
% zeta 1
fig = figure(11); clf;
ax = axes(fig);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [figW figH];
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
grid(ax, 'minor');
set(ax, 'FontName', 'Times New Roman');
set(ax, 'FontSize', ax_font_size);
set(ax, 'LineWidth', 1.1);
set(ax, 'TickLabelInterpreter', 'latex');

for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    color = color_list(ctrl_idx);

    if data.CTRL_INFO.CTRL_NUM ~= 3
        continue
    end
    
    t = data.t;
    zeta_hist = data.zeta_hist;

    plot(ax, t, zeta_hist(1,:), "Color", color, "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(ax, t, zeta_hist(2,:), "Color", color, "LineWidth", line_width, "LineStyle", "-."); hold on
end

%% ============================
% Computation Time
% ============================
fig = figure(12); clf;
ax = axes(fig);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [figW figH];
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');
grid(ax, 'minor');
set(ax, 'FontName', 'Times New Roman');
set(ax, 'FontSize', ax_font_size);
set(ax, 'LineWidth', 1.1);
set(ax, 'TickLabelInterpreter', 'latex');

for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    color = color_list(ctrl_idx);
    
    t = data.t;
    comp_time_hist = data.comp_time_hist;

    plot(ax, t, comp_time_hist, "Color", color, "LineWidth", line_width, "LineStyle", "-"); hold on
end

ax.XLim = [global_start_t global_end_t];
ax.XLabel.String = 'Time / s';
ax.YLabel.String = 'Comp. Time / $\mu$s';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';

%%============================
% % zoom (ctrl_obs_idx); q1, q2
% %  ===========================
ax_list{13}.XLim = [start_t end_t];
ax_list{14}.XLim = [start_t end_t];

maxVal1 = -inf; minVal1 = inf;
maxVal2 = -inf; minVal2 = inf;
for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    x1_hist = rad2deg(data.x1_hist);
    xd1_hist = rad2deg(data.xd1_hist);
    maxVal1 = max(maxVal1, max(x1_hist(1,ctrl_obs_idx))); maxVal1 = max(maxVal1, max(xd1_hist(1,ctrl_obs_idx)));
    minVal1 = min(minVal1, min(x1_hist(1,ctrl_obs_idx))); minVal1 = min(minVal1, min(xd1_hist(1,ctrl_obs_idx)));
    maxVal2 = max(maxVal2, max(x1_hist(2,ctrl_obs_idx))); maxVal2 = max(maxVal2, max(xd1_hist(2,ctrl_obs_idx)));
    minVal2 = min(minVal2, min(x1_hist(2,ctrl_obs_idx))); minVal2 = min(minVal2, min(xd1_hist(2,ctrl_obs_idx)));
end
ax_list{13}.YLim = [minVal1-5 maxVal1+5];
ax_list{14}.YLim = [minVal2-5 maxVal2+5];

%% ============================
% % zoom (ctrl_obs_idx); e1, e2
% %  ===========================
fig1 = figure(15); clf;
ax1 = axes(fig1);
set(fig1, 'Units', 'centimeters');
fig1.Position(3:4) = [figW figH];
hold(ax1, 'on');
grid(ax1, 'on');
box(ax1, 'on');
grid(ax1, 'minor');
set(ax1, 'FontName', 'Times New Roman');
set(ax1, 'FontSize', ax_font_size);
set(ax1, 'LineWidth', 1.1);
set(ax1, 'TickLabelInterpreter', 'latex');
fig2 = figure(16); clf;
ax2 = axes(fig2);
set(fig2, 'Units', 'centimeters');
fig2.Position(3:4) = [figW figH];
hold(ax2, 'on');
grid(ax2, 'on');
box(ax2, 'on');
grid(ax2, 'minor');
set(ax2, 'FontName', 'Times New Roman');
set(ax2, 'FontSize', ax_font_size);
set(ax2, 'LineWidth', 1.1);
set(ax2, 'TickLabelInterpreter', 'latex');

maxVal1 = -inf; minVal1 = inf;
maxVal2 = -inf; minVal2 = inf;
for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    x1_hist = rad2deg(data.x1_hist);
    xd1_hist = rad2deg(data.xd1_hist);
    e1_hist = x1_hist - xd1_hist;
    e1_hist = abs(e1_hist);

    t = data.t;

    plot(ax1, t, e1_hist(1,:), "Color", color_list(ctrl_idx), "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(ax2, t, e1_hist(2,:), "Color", color_list(ctrl_idx), "LineWidth", line_width, "LineStyle", "-"); hold on


    maxVal1 = max(maxVal1, max(e1_hist(1,ctrl_obs_idx))); 
    minVal1 = min(minVal1, min(e1_hist(1,ctrl_obs_idx))); 
    maxVal2 = max(maxVal2, max(e1_hist(2,ctrl_obs_idx))); 
    minVal2 = min(minVal2, min(e1_hist(2,ctrl_obs_idx))); 
end

yline(ax1, 0, 'k--', 'LineWidth', line_width)
yline(ax2, 0, 'k--', 'LineWidth', line_width)
xline(ax1, warmup_time+ep_time, 'k--', 'LineWidth', line_width)
xline(ax2, warmup_time+ep_time, 'k--', 'LineWidth', line_width)
minVal1 = 0; minVal2 = 0;

ax1.XLabel.String = 'Time / s';
ax1.YLabel.String = '$\vert q_1-{q_d}_1\vert$ / deg';
ax2.XLabel.String = 'Time / s';
ax2.YLabel.String = '$\vert q_2-{q_d}_2\vert$ / deg';
ax1.XLabel.Interpreter = 'latex';
ax1.YLabel.Interpreter = 'latex';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';

% ax1.XLim = [start_t end_t];
% ax2.XLim = [start_t end_t];
ax1.XLim = [global_start_t global_end_t];
ax2.XLim = [global_start_t global_end_t];

len1 = maxVal1 - minVal1; len2 = maxVal2 - minVal2;
ax1.YLim = [minVal1-len1*.1 maxVal1+len1*.1];
ax2.YLim = [minVal2-len2*.1 maxVal2+len2*.1];

text(ax1, 0.02, 0.9, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')
text(ax2, 0.02, 0.9, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')
text(ax1, 0.52, 0.9, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')
text(ax2, 0.52, 0.9, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman','Units','normalized')


%% ============================
% zoom Top view of the q
% ============================
fig = figure(17); clf;
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
grid(ax, 'minor');
box(ax, 'on');

set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [figW figH*1.5];
set(ax, 'FontName', 'Times New Roman');
set(ax, 'FontSize', ax_font_size);
set(ax, 'LineWidth', 1.1);
set(ax, 'TickLabelInterpreter', 'latex');
ax.XLabel.String = '$q_1-{q_d}_1$ / deg'; ax.XLabel.Interpreter = 'latex';
ax.YLabel.String = '$q_2-{q_d}_2$ / deg'; ax.YLabel.Interpreter = 'latex';

maxminX = [inf 0]; maxminY = [inf -inf];
for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    CTRL_INFO = data.CTRL_INFO;

    e1_hist = data.x1_hist-data.xd1_hist;
    e1_hist = rad2deg(e1_hist);
    e1_hist = [e1_hist, nan(2, length(t)-length(e1_hist))];
    e1_hist = e1_hist(:, ctrl_obs_idx);
    color = color_list(ctrl_idx);

    plot(ax, e1_hist(1,:), e1_hist(2,:), "Color", color, "LineWidth", line_width, "LineStyle", "-");
    % marker_idx = round(linspace(1, length(x1_hist), 5));
    % plot(ax, x1_hist(1,marker_idx), x1_hist(2,marker_idx), "Color", color, "Marker", "o", "MarkerSize", 6, "LineStyle", "none");

    maxminX = [min(maxminX(1), min(e1_hist(1,:))) max(maxminX(2), max(e1_hist(1,:)))];
    maxminY = [min(maxminY(1), min(e1_hist(2,:))) max(maxminY(2), max(e1_hist(2,:)))];
end
len = maxminX(2)-maxminX(1); ratio = .1;
ax.XLim = [maxminX(1)-len*ratio maxminX(2)+len*ratio];
len = maxminY(2)-maxminY(1); ratio = .1;
ax.YLim = [maxminY(1)-len*ratio maxminY(2)+len*ratio];

backstep_list = [7,60,10,10];

for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    CTRL_INFO = data.CTRL_INFO;

    e1_hist = data.x1_hist-data.xd1_hist;
    e1_hist = rad2deg(e1_hist);
    e1_hist = [e1_hist, nan(2, length(t)-length(e1_hist))];
    e1_hist = e1_hist(:, ctrl_obs_idx);
    color = color_list(ctrl_idx);

    backstep_for_arrow = backstep_list(ctrl_idx);
    smoothen_traj = smoothdata(e1_hist, 2, 'movmean', 10);
    % e1_hist(1,end-backstep_for_arrow), e1_hist(2,end-backstep_for_arrow),...
    drawArrow(ax,...
        smoothen_traj(1,end-backstep_for_arrow), smoothen_traj(2,end-backstep_for_arrow),...
        smoothen_traj(1,end), smoothen_traj(2,end),...
        color);
    hold on
end

yline(ax, 0, 'k--', 'LineWidth', line_width);
xline(ax, 0, 'k--', 'LineWidth', line_width);

% %% ============================
% % [Zoom] weight norm
% %   ============================
% fig = figure(17); clf;
% ax = axes(fig);
% set(fig, 'Units', 'centimeters');
% fig.Position(3:4) = [figW figH];
% hold(ax, 'on');
% grid(ax, 'on');
% box(ax, 'on');
% grid(ax, 'minor');
% set(ax, 'FontName', 'Times New Roman');
% set(ax, 'FontSize', ax_font_size);
% set(ax, 'LineWidth', 1.1);
% set(ax, 'TickLabelInterpreter', 'latex');   

% for ctrl_idx = 1:1:length(dataSet)
%     data = dataSet{ctrl_idx};
%     CTRL_INFO = data.CTRL_INFO;
%     color = color_list(ctrl_idx);
    
%     t = data.t;
%     th_hist = data.th_hist;
%     th0 = th_hist(1,:);
%     th1 = th_hist(2,:);
%     th2 = th_hist(3,:);

%     % th = [th0; th1; th2];
%     % del_th = th(2:end,:) - th(1:end-1,:);
%     % norm_th = vecnorm(del_th, 2, 1);
%     % plot(ax, t, norm_th, "Color", color, "LineWidth", line_width, "LineStyle", "-"); hold on

%     plot(ax, t, th0, "Color", color, "LineWidth", line_width, "LineStyle", "-", "HandleVisibility", "off"); hold on
%     plot(ax, t, th1, "Color", color, "LineWidth", line_width, "LineStyle", "-.", "HandleVisibility", "off"); hold on
%     plot(ax, t, th2, "Color", color, "LineWidth", line_width, "LineStyle", "--", "HandleVisibility", "off"); hold on
% end
% % dummy for legend (not plotted)
% % plot(ax, NaN, NaN, "Color", 'k', "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\Vert\widehat{\mbox{\boldmath $\theta$}}_0\Vert$"); hold on
% % plot(ax, NaN, NaN, "Color", 'k', "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\Vert\widehat{\mbox{\boldmath $\theta$}}_1\Vert$"); hold on
% % plot(ax, NaN, NaN, "Color", 'k', "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\Vert\widehat{\mbox{\boldmath $\theta$}}_2\Vert$"); hold on
% % plot(ax, [0 T], [+1 +1]*th_max(1), "Color", "black", "LineWidth", line_width, "LineStyle", "-.", "HandleVisibility", "off"); hold on
% % lgd = legend(ax);
% % lgd.Location = 'northwest';
% % lgd.Interpreter = 'latex';
% % lgd.FontSize = 10;
% % lgd.NumColumns = 3;

% ax.XLim = [start_t end_t];
% ax.YLim = [2 4];
% ax.XLabel.String = 'Time / s';
% ax.YLabel.String = '$\Vert\widehat{\mbox{\boldmath $\theta$}}\Vert$';
% ax.XLabel.Interpreter = 'latex';
% ax.YLabel.Interpreter = 'latex';


%%
bar_plotter

%% SAVE FIGURES
if SAVE_FLAG
    FIG_SAVE_PATH = "figures";
    [~,~] = mkdir(FIG_SAVE_PATH);

    for idx = [1:1:17, 21]

        f_name = FIG_SAVE_PATH + "/Fig" + string(idx);

        saveas(figure(idx), f_name + ".png")

        figure(idx);
        % set(gcf, 'Position', [0, 0, fig_width, fig_height]); % [left, bottom, width, height] 
        exportgraphics(gcf, f_name+'.eps', 'ContentType', 'vector')
        % exportgraphics(figure(idx), f_name+'.eps',"Padding","figure")
    
        % matlab2tikz(char(f_name+".tex"))

        fprintf("Saved Figure %d\n", idx)
    end
end


%% NUMERICAL ANALYSIS
ctrl_dt = 1/250;
sim_dt = ctrl_dt / 1000;

ep1_idx = find(t >= warmup_time & t <= warmup_time + ep_time);
ep2_idx = find(t >= warmup_time + ep_time & t <= warmup_time + 2*ep_time);

for c_idx = 1:1:length(dataSet)
    data = dataSet{c_idx};

    x1_hist = data.x1_hist;
    x2_hist = data.x2_hist;
    xd1_hist = data.xd1_hist;
    xd2_hist = data.xd2_hist;
    e1_hist = x1_hist - xd1_hist;
    e2_hist = x2_hist - xd2_hist;

    Lambda = diag([5 15]);
    r = e2_hist + Lambda*e1_hist;

    dataSet{c_idx}.r_hist = r(1,:);
end

RMSE = @(e) sqrt(mean(e.^2, 2));

fprintf("%% RMSE: \n")
for c_idx = 1:1:length(dataSet)
    data = dataSet{c_idx};

    if data.t(end) < warmup_time + 2*ep_time
        fprintf("%% C %s : failed to control\n", tex_name_list{c_idx})
        continue
    end

    e_hist = rad2deg(data.x1_hist - data.xd1_hist);
    r_hist = rad2deg(data.r_hist);

    e1_ep1 = e_hist(1, ep1_idx);
    e2_ep1 = e_hist(2, ep1_idx);
    e1_ep2 = e_hist(1, ep2_idx);
    e2_ep2 = e_hist(2, ep2_idx);
    r_ep1  = r_hist(ep1_idx);
    r_ep2  = r_hist(ep2_idx);

    ctrl_name = tex_name_list{c_idx};

    fprintf("%% -------------------------------\n")
    fprintf("%% C %s :\n", ctrl_name)
    fprintf("%% error q1 (1) %.3f, (2) %.3f [imp: %.3f]\n", RMSE(e1_ep1), RMSE(e1_ep2), 1-RMSE(e1_ep2)/RMSE(e1_ep1))
    fprintf("%% error q2 (1) %.3f, (2) %.3f [imp: %.3f]\n", RMSE(e2_ep1), RMSE(e2_ep2), 1-RMSE(e2_ep2)/RMSE(e2_ep1))
    fprintf("%% filter error   (1) %.3f, (2) %.3f [imp: %.3f]\n", RMSE(r_ep1), RMSE(r_ep2), 1-RMSE(r_ep2)/RMSE(r_ep1))
    fprintf("\n")

    fprintf("\\newcommand{\\%sEpiOneAngOne}{%.3f}\n", ctrl_name, RMSE(e1_ep1))
    fprintf("\\newcommand{\\%sEpiOneAngTwo}{%.3f}\n", ctrl_name, RMSE(e2_ep1))
    fprintf("\\newcommand{\\%sEpiTwoAngOne}{%.3f}\n", ctrl_name, RMSE(e1_ep2))
    fprintf("\\newcommand{\\%sEpiTwoAngTwo}{%.3f}\n", ctrl_name, RMSE(e2_ep2))
    fprintf("\\newcommand{\\%sImpAngOne}{%.3f}\n", ctrl_name, (RMSE(e1_ep2)/RMSE(e1_ep1)-1)*100)
    fprintf("\\newcommand{\\%sImpAngTwo}{%.3f}\n", ctrl_name, (RMSE(e2_ep2)/RMSE(e2_ep1)-1)*100)
end

    fprintf("%% -------------------------------\n")
    fprintf("\\newcommand{\\zoomStartTime}{%.1f}\n", start_t)
    fprintf("\\newcommand{\\zoomEndTime}{%.1f}\n", end_t)
beep()

function data = loadFromMeas(data_path, CTRL_INFO)
    data_orig = load(data_path);

    
    data = struct();
    data.CTRL_INFO.CTRL_NUM = CTRL_INFO;
    
    % ref start to end idx
    xd1_hist = data_orig.r(1:2, :);
    ref_start_idx = find(isnan(xd1_hist(1,:))==0, 1);
    ref_end_idx = find(isnan(xd1_hist(1,:))==0, 1, 'last');
    
    TRUE_CTRL_NUM = data_orig.conac_controller_id(ref_start_idx);
    fprintf("TRUE_CTRL_NUM: %d\n", TRUE_CTRL_NUM)
    
    % data mapping
    data.x1_hist = data_orig.q(1:2, ref_start_idx:ref_end_idx);
    data.x2_hist = data_orig.qdot(1:2, ref_start_idx:ref_end_idx);
    data.xd1_hist = data_orig.r(1:2, ref_start_idx:ref_end_idx);
    data.xd2_hist = data_orig.rdot(1:2, ref_start_idx:ref_end_idx);
    data.u_hist = data_orig.u_command(1:2, ref_start_idx:ref_end_idx);
    data.uSat_hist = data_orig.u_saturated(1:2, ref_start_idx:ref_end_idx);
    % data.xd3_hist = data_orig.rddot(1:2, :);

    data.lbd_hist = data_orig.lbd(:, ref_start_idx:ref_end_idx);
    data.th_hist = data_orig.Vn(:, ref_start_idx:ref_end_idx);
    data.zeta_hist = data_orig.zeta(:, ref_start_idx:ref_end_idx);
    
    if isfield(data_orig, 'conac_computation_time_us')
        data.comp_time_hist = data_orig.conac_computation_time_us(ref_start_idx:ref_end_idx);
    else
        data.comp_time_hist = zeros(1, ref_end_idx-ref_start_idx+1);
    end

    t = data_orig.t(ref_start_idx:ref_end_idx);
    data.t = t - t(1);
    data.T = data.t(end);
end

function drawArrow(ax, x1, y1, x2, y2, color)

    % axes position (normalized in figure)
    axPos = ax.Position;

    % axis limits
    xl = xlim(ax);
    yl = ylim(ax);

    % data -> normalized axes
    xn1 = (x1-xl(1))/(xl(2)-xl(1));
    yn1 = (y1-yl(1))/(yl(2)-yl(1));

    xn2 = (x2-xl(1))/(xl(2)-xl(1));
    yn2 = (y2-yl(1))/(yl(2)-yl(1));

    % axes -> figure
    xf1 = axPos(1) + xn1*axPos(3);
    yf1 = axPos(2) + yn1*axPos(4);

    xf2 = axPos(1) + xn2*axPos(3);
    yf2 = axPos(2) + yn2*axPos(4);

    annotation(gcf,'arrow',...
        [xf1 xf2],...
        [yf1 yf2],...
        'Color',color,...
        'LineWidth',1.5,...
        'HeadLength',14,...
        'HeadWidth',14);

end