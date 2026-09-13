function r_func = ref4_load()

r_func = @(t) epi_gen2(t);

end

function [r1, r2, r3] = epi_gen1(t)
    T = 3; % traj. duration
    POINT_NUM = 4; % number of points in one cycle
    CYCLE_TIME = T*POINT_NUM;
    t = mod(t, CYCLE_TIME);

    x0 = deg2rad([-60;60]);  
    xd1 = deg2rad([45;-90]);
    xd2 = deg2rad([-45;45]);
    % xd3 = deg2rad([-50;-145]);
    xd3 = xd1;

    if t < T
        [r1,r2,r3] = poly_filter(x0, xd1, T, t);
    elseif t < 2*T
        [r1,r2,r3] = poly_filter(xd1, xd2, T, (t-T));
    elseif t < 3*T
        [r1,r2,r3] = poly_filter(xd2, xd3, T, (t-2*T));
    else
        [r1,r2,r3] = poly_filter(xd3, x0, T, (t-3*T));
    end
end

function [r1, r2, r3] = epi_gen2(t)
    WARMUP_T = 3; % warmup duration
    COOL_T = 8; % cool down duration after warmup (for smooth transition)
    UNREACHABLE_T = 3; % duration of unreachable desired state (for showing KKT satisfaction at s.s.)

    INIT_X = [-1/2*pi; 0]; % initial state for warmup (link angle)
    WARMUP_X = [-1/3*pi; 1/3*pi]; % same as x0
    UNREACHABLE_X = deg2rad([45; -90]);

    if t < WARMUP_T
        [r1, r2, r3] = poly_filter(INIT_X, WARMUP_X, WARMUP_T, t);
    elseif t < WARMUP_T + COOL_T
        r1 = WARMUP_X;
        r2 = [0;0];
        r3 = [0;0];
    elseif t < WARMUP_T + COOL_T + 24 % (24s normal reference duration)
        [r1, r2, r3] = epi_gen1(t-WARMUP_T-COOL_T);
    elseif t < WARMUP_T + COOL_T + 24 + UNREACHABLE_T
        [r1, r2, r3] = poly_filter( ...
            WARMUP_X, ...
            UNREACHABLE_X, ...
            UNREACHABLE_T, t-WARMUP_T-COOL_T-24);
    else 
        r1 = UNREACHABLE_X;
        r2 = [0;0];
        r3 = [0;0];
    end
end

function [r1, r2, r3] = poly_filter(x0, xd, T, t)

a0 = x0;
a1 = zeros(2,1);
a2 = zeros(2,1);
a3 = 10*(xd-a0)/T^3;
a4 = -15*(xd-x0)/T^4;
a5 = 6*(xd-a0)/T^5;

r1 = [a0 a1 a2 a3 a4 a5] * [1 t t^2 t^3 t^4 t^5]';
r2 = [a1 a2 a3 a4 a5] * [1 2*t 3*t^2 4*t^3 5*t^4]';
r3 = [a2 a3 a4 a5] * [2 6*t 12*t^2 20*t^3]';
end