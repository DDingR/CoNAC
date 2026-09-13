%% NUMERICAL ANALYSIS
ctrl_dt = 1/250;
sim_dt  = ctrl_dt/1000;

ep1_idx = find(t >= warmup_time & ...
               t <= warmup_time + ep_time);

ep2_idx = find(t >= warmup_time + ep_time & ...
               t <= warmup_time + 2*ep_time);

Lambda = diag([5 15]);

%% FILTERED ERROR CALCULATION
for c_idx = 1:length(dataSet)
    data = dataSet{c_idx};

    x1_hist  = data.x1_hist;
    x2_hist  = data.x2_hist;
    xd1_hist = data.xd1_hist;
    xd2_hist = data.xd2_hist;

    e1_hist = x1_hist - xd1_hist;
    e2_hist = x2_hist - xd2_hist;

    r = e2_hist + Lambda*e1_hist;

    dataSet{c_idx}.r_hist = r(1,:);
end

%% RMSE CALCULATION
RMSE = @(e) sqrt(mean(e.^2, 2));

num_ctrl = length(dataSet);

rmse_q1 = nan(num_ctrl, 2);
rmse_q2 = nan(num_ctrl, 2);

valid_ctrl = false(num_ctrl, 1);

fprintf("%% RMSE:\n");

for c_idx = 1:num_ctrl
    data = dataSet{c_idx};

    if data.t(end) < warmup_time + 2*ep_time
        fprintf("%% C %s : failed to control\n", name_list{c_idx});
        continue
    end

    valid_ctrl(c_idx) = true;

    e_hist = rad2deg(data.x1_hist - data.xd1_hist);
    r_hist = rad2deg(data.r_hist);

    e1_ep1 = e_hist(1, ep1_idx);
    e2_ep1 = e_hist(2, ep1_idx);

    e1_ep2 = e_hist(1, ep2_idx);
    e2_ep2 = e_hist(2, ep2_idx);

    r_ep1 = r_hist(ep1_idx);
    r_ep2 = r_hist(ep2_idx);

    ctrl_name = name_list{c_idx};

    rmse_q1(c_idx,1) = RMSE(e1_ep1);
    rmse_q1(c_idx,2) = RMSE(e1_ep2);

    rmse_q2(c_idx,1) = RMSE(e2_ep1);
    rmse_q2(c_idx,2) = RMSE(e2_ep2);
end

%% REMOVE FAILED CONTROLLERS
rmse_q1 = rmse_q1(valid_ctrl,:);
rmse_q2 = rmse_q2(valid_ctrl,:);
ctrl_names = name_list(valid_ctrl);

%% REVERSE ORDER: CONAC_HIGH -> ... -> NAC
rmse_q1 = flipud(rmse_q1);
rmse_q2 = flipud(rmse_q2);
ctrl_names = flip(ctrl_names);

%% FIGURE
fig = figure(21);
clf(fig);

set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [17 10];

ax1 = subplot(2,1,1, 'Parent', fig);
plotGroupedRMSE( ...
    ax1, ...
    rmse_q1, ...
    ctrl_names, ...
    '$\mathrm{RMSE}(q_1-{q_d}_1)$ / deg');

ax2 = subplot(2,1,2, 'Parent', fig);
plotGroupedRMSE( ...
    ax2, ...
    rmse_q2, ...
    ctrl_names, ...
    '$\mathrm{RMSE}(q_2-{q_d}_2)$ / deg');


%% LOCAL FUNCTION
function plotGroupedRMSE(ax, rmse_data, ctrl_names, y_label_text)

    cla(ax);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    ax.XGrid = 'off';
    ax.YGrid = 'on';

    num_ctrl = size(rmse_data, 1);

    %% BAR POSITIONS
    bar_width = 0.4;
    pair_gap  = 0.00;
    ctrl_gap  = 0.2;

    pair_width = 2*bar_width + pair_gap;

    pair_center = ...
        (0:num_ctrl-1)*(pair_width + ctrl_gap) + 1;

    x_ep1 = pair_center - (bar_width + pair_gap)/2;
    x_ep2 = pair_center + (bar_width + pair_gap)/2;

    %% BAR PLOTS
    b1 = bar( ...
        ax, ...
        x_ep1, ...
        rmse_data(:,1), ...
        bar_width, ...
        'FaceColor', [0.7 0.7 0.7], ...
        'DisplayName', 'Episode 1');

    b2 = bar( ...
        ax, ...
        x_ep2, ...
        rmse_data(:,2), ...
        bar_width, ...
        'FaceColor', [0.3 0.3 0.3], ...
        'DisplayName', 'Episode 2');

    %% AXIS SETTINGS
    ax.XTick = pair_center;
    ax.XTickLabel = ctrl_names;

    ax.FontName = 'Times New Roman';
    ax.FontSize = 10;
    ax.LineWidth = 1.0;
    ax.TickLabelInterpreter = 'latex';

    ylabel( ...
        ax, ...
        y_label_text, ...
        'Interpreter', 'latex', ...
        'FontSize', 11);

    %% Y-LIMIT
    valid_values = rmse_data(~isnan(rmse_data));

    if isempty(valid_values)
        max_val = 1;
    else
        max_val = max(valid_values);
    end

    if max_val <= 0
        max_val = 1;
    end

    ylim(ax, [0, 1.5*max_val]);

    value_offset  = 0.025*max_val;
    change_offset = 0.14*max_val;

    %% VALUE LABELS
    for c_idx = 1:num_ctrl
        ep1 = rmse_data(c_idx,1);
        ep2 = rmse_data(c_idx,2);

        text( ...
            ax, ...
            x_ep1(c_idx), ...
            ep1 + value_offset, ...
            sprintf('%.3f', ep1), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 9);

        text( ...
            ax, ...
            x_ep2(c_idx), ...
            ep2 + value_offset, ...
            sprintf('%.3f', ep2), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 9);

        %% IMPROVEMENT RATE
        if ep1 ~= 0
            improvement = ...
                -(1 - ep2/ep1)*100;
        else
            improvement = NaN;
        end

        pair_top = max(ep1, ep2);

        % pair_center(c_idx), ...
        % pair_top + change_offset, ...
        text( ...
            ax, ...
            x_ep2(c_idx), ...
            0, ...
            sprintf('%+.1f\\%%', improvement), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'Interpreter', 'latex', ...
            'FontName', 'Times New Roman', ...
            'Color', [1 1 1], ...
            'FontSize', 9);
    end

    %% LEGEND
    legend( ...
        ax, ...
        [b1 b2], ...
        'Location', 'northwest', ...
        'Interpreter', 'latex', ...
        'Box', 'off');

    %% X-LIMIT
    xlim( ...
        ax, ...
        [x_ep1(1) - bar_width, ...
         x_ep2(end) + bar_width]);
end