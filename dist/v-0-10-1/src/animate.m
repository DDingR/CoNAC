
figure(30); close
fig = figure(7);
% fig.Units = 'Pixels';
fig.Color = 'black';
% fig.Renderer = 'opengl';
fig.Position(3:4) = [800 600];
ax1         = subplot(3,4, [1,2,5,6], 'Parent', fig); 
ax_traj1    = subplot(3,4, [9,10], 'Parent', fig);
ax_traj2    = subplot(3,4, [11,12], 'Parent', fig);
ax_input    = subplot(3,4, [3,4,7,8], 'Parent', fig);

ax_traj1.Position(1) = ax_traj1.Position(1) + 0.04;
ax_traj2.Position(1) = ax_traj2.Position(1) + 0.04;
ax_traj1.Position(3) = ax_traj1.Position(3) - 0.04;
ax_traj2.Position(3) = ax_traj2.Position(3) - 0.04;
ax_input.Position(1) = ax_input.Position(1) + 0.02;
ax_input.Position(3) = ax_input.Position(3) - 0.03;

AINMATION_SAVE_FLAG = 1;
video_name = "sample";

line_width = 1.5;
font_size = 16;
% accel = 1;
accel = 20;
dt = t(2) - t(1);


%%
WARMUP_AND_IDLE_TIME = 14; % warmup and idle time
EPISODE_TIME = 8; % duration of each episode (after warmup and cool down)

%%
% T = dataSet{1}.T;
T = 34; % s
t = dataSet{1}.t;
obs_t = 1:length(t);

animation_end_t = WARMUP_AND_IDLE_TIME + EPISODE_TIME*3;
% assert(T == animation_end_t, "T does not match the expected animation end time.");

PhaseList = {
    "Warm-up",
    "Episode 1",
    "Episode 2",
    "Ending"
};

%%
if AINMATION_SAVE_FLAG
    v = VideoWriter(video_name, 'MPEG-4');
    % v = VideoWriter(video_name, 'Uncompressed AVI');
    % v.Quality = 100;
    v.FrameRate = 1/dt/accel; 
    open(v);
end

% fig = figure(7); clf   

REPRESENTATIVE_CTRL_IDX = 4; % index of the representative control method for animation
data = dataSet{REPRESENTATIVE_CTRL_IDX};
CTRL_INFO = data.CTRL_INFO;

% th_max = opt.cstr.th_max;
% u_ball = opt.cstr.u_ball;
u_max2 = uMax2;
u_max1 = uMax1;

%%
PLOT_WINDOW_STATE = 2; % s
PLOT_WINDOW_INPUT = 1; % s

%%
max_t_num = 0;
for ctrl_idx = 1:1:length(dataSet)
    data = dataSet{ctrl_idx};
    x1_hist = data.x1_hist;
    c = color_list{ctrl_idx};
    t = data.t;
    hold(ax_traj1, "on"); hold(ax_traj2, "on");
    plot(ax_traj1, t, x1_hist(1,:), "color", c, "LineWidth", line_width, "LineStyle", "-"); 
    plot(ax_traj2, t, x1_hist(2,:), "color", c, "LineWidth", line_width, "LineStyle", "-"); 

    if ctrl_idx == REPRESENTATIVE_CTRL_IDX
        xd1_hist = data.xd1_hist;
        plot(ax_traj1, t, xd1_hist(1,:), "color", "red", "LineWidth", line_width, "LineStyle", "--"); 
        plot(ax_traj2, t, xd1_hist(2,:), "color", "red", "LineWidth", line_width, "LineStyle", "--"); 
    end

    max_t_num = max(max_t_num, length(t));
end

%%
L1 = .45;            
L2 = .45;

MAXMIN_VALS = [-inf(1, 3); inf(1, 3)];
for cur_t_idx = 1:accel:max_t_num
% for cur_t_idx = 1

    cur_t = t(cur_t_idx);
    
    % index for state
    PLOT_WINDOW_STATE_END = min(cur_t + PLOT_WINDOW_STATE, T);
    plot_t_idx_state = find(t >= cur_t - PLOT_WINDOW_STATE, 1):find(t >= PLOT_WINDOW_STATE_END, 1);
    plot_t_idx_input = find(t >= cur_t - PLOT_WINDOW_INPUT, 1):find(t >= cur_t, 1);

    t_idx = find(t >= cur_t, 1);
    if isempty(t_idx) || isempty(plot_t_idx_state) || isempty(plot_t_idx_input) || cur_t > T
        break
    end

    % main plot
    for ctrl_idx = 1:1:length(dataSet)

        data = dataSet{ctrl_idx};
        CTRL_INFO = data.CTRL_INFO;
        
        c = color_list{ctrl_idx};

        x1_hist = data.x1_hist; 
        u_hist = data.u_hist;
        uSat_hist = data.uSat_hist;

        q1 = x1_hist(1, t_idx); q2 = x1_hist(2, t_idx);
        
        plot(ax1, ...
            L1* [0 cos(q1), cos(q1)+cos(q1+q2)], ...
            L1* [0 sin(q1), sin(q1)+sin(q1+q2)], ...
            "Color", c, ...
            "LineStyle","-", ...
            "Marker","o", ...
            "LineWidth", line_width ...
            ); hold on

        if CTRL_INFO.CTRL_NUM == 1
            r_hist = data.r_hist;

            r1 = xd1_hist(1, t_idx); r2 = xd1_hist(2, t_idx);
            plot(ax1, ...
                L1* [0 cos(r1), cos(r1)+cos(r1+r2)], ...
                L1* [0 sin(r1), sin(r1)+sin(r1+r2)], ...
                "Color", "red", ...
                "LineStyle","--", ...
                "Marker","o", ...
                "LineWidth", line_width ...
                ); hold on

            % time 
            text(ax1, .1,.8,0, ...
                sprintf( ...
                    'Time: %.2f/%.0f s (%.1f%%) ', t(t_idx), T, round(t(t_idx)/T*100, 3) ...
                ), ...
                "FontSize", font_size, ...
                "FontName", "Times New Roman", ...
                "Units", "normalized" ...
                );
            if cur_t < WARMUP_AND_IDLE_TIME
                phase_num = 1;
            elseif cur_t < WARMUP_AND_IDLE_TIME + EPISODE_TIME
                phase_num = 2;
            elseif cur_t < WARMUP_AND_IDLE_TIME + EPISODE_TIME*2
                phase_num = 3;
            else
                phase_num = 4;
            end
            text(ax1, .1,.75,0, ...
                sprintf( ...
                    'Phase: %s', PhaseList{phase_num} ...
                ), ...
                "FontSize", font_size, ...
                "FontName", "Times New Roman", ...
                "Units", "normalized" ...
                );

            xlabel("$x$", "Interpreter", "latex")
            ylabel("$y$", "Interpreter","latex")
            set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
        end

        % mini plots for state
        MAXMIN_VALS(1, 1) = max(MAXMIN_VALS(1,1), max(x1_hist(1, plot_t_idx_state)));
        MAXMIN_VALS(2, 1) = min(MAXMIN_VALS(2,1), min(x1_hist(1, plot_t_idx_state)));
        MAXMIN_VALS(1, 2) = max(MAXMIN_VALS(1,2), max(x1_hist(2, plot_t_idx_state)));
        MAXMIN_VALS(2, 2) = min(MAXMIN_VALS(2,2), min(x1_hist(2, plot_t_idx_state)));
        MAXMIN_VALS(1, 1) = max(MAXMIN_VALS(1,1), max(xd1_hist(1, plot_t_idx_state)));
        MAXMIN_VALS(2, 1) = min(MAXMIN_VALS(2,1), min(xd1_hist(1, plot_t_idx_state)));
        MAXMIN_VALS(1, 2) = max(MAXMIN_VALS(1,2), max(xd1_hist(2, plot_t_idx_state)));
        MAXMIN_VALS(2, 2) = min(MAXMIN_VALS(2,2), min(xd1_hist(2, plot_t_idx_state)));
        xline(ax_traj1, cur_t, "color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
        xline(ax_traj2, cur_t, "color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on

        % mini plot for input ball
        ang = 0:0.01:2*pi;    
        cur_bar_q1 = plot(ax_input, u_hist(1,plot_t_idx_input ), u_hist(2,plot_t_idx_input ), "color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
        cur_bar_q2 = plot(ax_input, u_hist(1,t_idx), u_hist(2,t_idx), "color", c,"Marker","o"); hold on

        plot(ax_input, u_ball*cos(ang), u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
        yline(ax_input, uMax2, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
        yline(ax_input, -uMax2, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
        xline(ax_input, uMax1, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
        xline(ax_input, -uMax1, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on

    end

    for ax_ = [ax1, ax_traj1, ax_traj2, ax_input]
        grid(ax_, "on"); box(ax_, "on"); hold(ax_, "on");
    end

    % axis settings
    ax1.XLim = [-1.2 1.2];
    ax1.YLim = [-1.2 1.2];
    % pbaspect(ax1, [1 1 1])
    ax1.XTickLabel = [];
    ax1.YTickLabel = [];

    ax_traj1.XLim = [cur_t-PLOT_WINDOW_STATE, cur_t+PLOT_WINDOW_STATE];
    ax_traj2.XLim = [cur_t-PLOT_WINDOW_STATE, cur_t+PLOT_WINDOW_STATE];
    len_x1 = MAXMIN_VALS(1,1) - MAXMIN_VALS(2,1);
    len_x2 = MAXMIN_VALS(1,2) - MAXMIN_VALS(2,2);
    ax_traj1.YLim = [MAXMIN_VALS(2,1)-0.2*len_x1, MAXMIN_VALS(1,1)+0.2*len_x1];
    ax_traj2.YLim = [MAXMIN_VALS(2,2)-0.2*len_x2, MAXMIN_VALS(1,2)+0.2*len_x2];
    ax_traj1.XLabel.String = "Time / s";
    ax_traj2.XLabel.String = "Time / s";
    ax_traj1.YLabel.String = "$q_1$ / rad";
    ax_traj2.YLabel.String = "$q_2$ / rad";
    ax_traj1.XLabel.Interpreter = "latex";
    ax_traj2.XLabel.Interpreter = "latex";
    ax_traj1.YLabel.Interpreter = "latex";
    ax_traj2.YLabel.Interpreter = "latex";
    ax_traj1.XLabel.Color = 'white';
    ax_traj1.YLabel.Color = 'white';
    ax_traj2.XLabel.Color = 'white';
    ax_traj2.YLabel.Color = 'white';
    ax_traj1.XColor = 'white';
    ax_traj1.YColor = 'white';
    ax_traj2.XColor = 'white';
    ax_traj2.YColor = 'white';
    ax_traj1.GridColor = 'black'; ax_traj1.MinorGridColor = 'black';
    ax_traj2.GridColor = 'black'; ax_traj2.MinorGridColor = 'black';
    set(ax_traj1, 'FontSize', font_size, 'FontName', 'Times New Roman')
    set(ax_traj2, 'FontSize', font_size, 'FontName', 'Times New Roman')

    ax_input.XLabel.String = "$\tau_1$ / Nm";
    ax_input.YLabel.String = "$\tau_2$ / Nm";
    ax_input.XLabel.Interpreter = "latex";
    ax_input.YLabel.Interpreter = "latex";
    ax_input.XLabel.Color = 'white';
    ax_input.YLabel.Color = 'white';
    ax_input.XColor = 'white';
    ax_input.YColor = 'white';  
    ax_input.GridColor = 'black'; ax_input.MinorGridColor = 'black';
    set(ax_input, 'FontSize', font_size, 'FontName', 'Times New Roman')
    mul_range = 1.2;
    % ax_input.XLim = mul_range* [-u_ball, u_ball];
    % ax_input.YLim = mul_range* [-u_ball, u_ball];
    ax_input.XLim = [-2 13];
    ax_input.YLim = [-3.5 13.5];
    % pbaspect(ax_input, [1 1 1])

    % get frame
    drawnow
    pause(0.01)

    if AINMATION_SAVE_FLAG
        f = getframe(gcf);
        writeVideo(v, f);
    end

    cla(ax1); cla(ax_input);
    delete(findall(ax_traj1, 'Type', 'ConstantLine'));
    delete(findall(ax_traj2, 'Type', 'ConstantLine'));
end

if AINMATION_SAVE_FLAG
    close(v);
end
