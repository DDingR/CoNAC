
figure(7); clf

L1 = .45;            
L2 = .45;

accel = 100;

if AINMATION_SAVE_FLAG
    v = VideoWriter("sim_result/"+ctrl_name, 'MPEG-4');
    % v.Quality = 100;
    v.FrameRate = 1/dt/accel; 
    open(v);
end

for t_idx = 1:accel:num_t
    clf

    q1 = x_hist(1, t_idx); q2 = x_hist(2, t_idx);
    r1 = r_hist(1, t_idx); r2 = r_hist(2, t_idx);

    plot( ...
        L1* [0 cos(r1), cos(r1)+cos(r1+r2)], ...
        L1* [0 sin(r1), sin(r1)+sin(r1+r2)], ...
        "Color", "red", ...
        "LineStyle","-", ...
        "Marker","o", ...
        "LineWidth", line_width ...
        ); hold on

    plot( ...
        L1* [0 cos(q1), cos(q1)+cos(q1+q2)], ...
        L1* [0 sin(q1), sin(q1)+sin(q1+q2)], ...
        "Color", "blue", ...
        "LineStyle","-", ...
        "Marker","o", ...
        "LineWidth", line_width ...
        ); hold on

    % time 
    text(-.95,+.8,0, ...
        sprintf( ...
            't [s]: %.2f\t/%.2f (%.1f%%) ', t(t_idx), T, round(t(t_idx)/T*100, 3) ...
        ), ...
        "FontSize", font_size, ...
        "FontName", "Times New Roman" ...
        );

    xlabel("$x$", "Interpreter", "latex")
    ylabel("$y$", "Interpreter","latex")
    set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')

    grid on
    xlim([-1.2 1.2]);
    ylim([-1.2 1.2]);
    pbaspect([1 1 1])

    % input ball
    axes('Position',[.6 .22 .2 .2])
    box on
    ang = 0:0.01:2*pi;
    shadow = max(1, t_idx-1e3:1:t_idx);
    
    plot(u_hist(1,shadow ), u_hist(2,shadow ), "color", 'red', "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(u_hist(1,t_idx), u_hist(2,t_idx), "color", 'red',"Marker","o"); hold on

    plot(uSat_hist(1,shadow ), uSat_hist(2,shadow ), "color", 'blue', "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(uSat_hist(1,t_idx), uSat_hist(2,t_idx), "color", 'blue',"Marker","o"); hold on
    
    plot(nnOpt.cstr.u_ball*cos(ang), nnOpt.cstr.u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
    
    xlabel("$\tau_1$", "Interpreter", "latex")
    ylabel("$\tau_2$", "Interpreter", "latex")
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')

    grid on
    mul_range = 2;
    xlim(mul_range* [-nnOpt.cstr.u_ball, nnOpt.cstr.u_ball])
    ylim(mul_range* [-nnOpt.cstr.u_ball, nnOpt.cstr.u_ball])
    pbaspect([1 1 1])

    % get frame
    drawnow

    if AINMATION_SAVE_FLAG
        f = getframe(gcf);
        writeVideo(v, f);
    end

end

if AINMATION_SAVE_FLAG
    close(v);
end
