function [r1, r2] = ref3_load()

max_ang = [ % rad
    3;
    3
];

r1 = @(x, t) ref1_gen(x, t) .* max_ang - [deg2rad(90); 0];
r2 = @(x, t) ref2_gen(x, t) .* max_ang;

end

function r1 = ref1_gen(~, t)
    cycle_time = 4;
    freq = 4/4/2*pi;

    rep_num = 2;
    episode_time = rep_num * cycle_time;

    % k = fix( t/cycle_time );
    % resetIndex = fix(episode_time/cycle_time);
    % k = mod(k, rep_num);

    t = mod(t, episode_time);

    r1 = [
        -1*sin(freq*t) 
        +1*sin(freq*t) 
    ] * t/episode_time;
end

function r2 = ref2_gen(~, t)
    cycle_time = 4;
    freq = 4/4/2*pi;

    rep_num = 3;
    episode_time = rep_num * cycle_time;

    % k = fix( t/cycle_time );
    % resetIndex = fix(episode_time/cycle_time);
    % k = mod(k, rep_num);

    t = mod(t, episode_time);

    r2 = [
        -1*cos(freq*t) * freq
        +1*cos(freq*t) * freq
    ] * t/episode_time;
end