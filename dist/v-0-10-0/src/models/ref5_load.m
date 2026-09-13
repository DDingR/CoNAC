function r_func = ref5_load()
% ref5_load - Reference signal generator for simulation
% Syntax:  r_func = ref5_load()
% Inputs:
%
% Outputs:
%   r_func   - (function handle) Reference signal generator function
%

r_func = @(t) ref_gen(t);

end

function [r1, r2] = ref_gen(t)
% 

    % warmup duration and cool down duration
    WARMUP_T = 3; % warmup duration
    COOL_T = 8; % cool down duration after warmup (for smooth transition)
    EPR_T = 12; % duration of each episode (after warmup and cool down)

    if t < WARMUP_T + COOL_T
        [r1, r2] = warmup_gen(t, WARMUP_T);
    elseif t < WARMUP_T + COOL_T + EPR_T
        [r1, r2] = epi1_gen(t-WARMUP_T-COOL_T);
    elseif t < WARMUP_T + COOL_T + EPR_T*2
        [r1, r2] = epi1_gen(t-WARMUP_T-COOL_T-EPR_T);
    else
        [r1, r2] = epi1_gen(t-WARMUP_T-COOL_T-EPR_T*2);
    end

end

function [r1, r2] = warmup_gen(t, WARMUP_T)
    INIT_X      = deg2rad([-90; 0]); % initial state for warmup (link angle)
    WARMUP_X    = deg2rad([-60; 60]); % same as x0

    if t < WARMUP_T
        [r1, r2] = poly_filter(INIT_X, WARMUP_X, WARMUP_T, t);
    else
        r1 = WARMUP_X;
        r2 = [0;0];
    end
end

function [r1, r2] = epi1_gen(t)
    T_list = [3, 3, 3, 3];  % duration for each segment
    pt_list = {             % points for each segment
        deg2rad([-60; 60]), ...
            deg2rad([30; -90]) ...
            deg2rad([60; -60]) ...
            deg2rad([-60; 45+60]), ...
    };

    if t < T_list(1)
        [r1,r2] = poly_filter(pt_list{1}, pt_list{2}, T_list(1), t);
    elseif t < sum(T_list(1:2))
        [r1,r2] = poly_filter(pt_list{2}, pt_list{3}, T_list(2), t - T_list(1));
    elseif t < sum(T_list(1:3))
        [r1,r2] = poly_filter(pt_list{3}, pt_list{4}, T_list(3), t - sum(T_list(1:2)));
    else
        [r1,r2] = poly_filter(pt_list{4}, pt_list{1}, T_list(4), t - sum(T_list(1:3)));
    end
end


function [r1, r2] = poly_filter(x0, xd, T, t)

a0 = x0;
a1 = zeros(2,1);
a2 = zeros(2,1);
a3 = 10*(xd-a0)/T^3;
a4 = -15*(xd-x0)/T^4;
a5 = 6*(xd-a0)/T^5;

r1 = [a0 a1 a2 a3 a4 a5] * [1 t t^2 t^3 t^4 t^5]';
r2 = [a1 a2 a3 a4 a5] * [1 2*t 3*t^2 4*t^3 5*t^4]';

end