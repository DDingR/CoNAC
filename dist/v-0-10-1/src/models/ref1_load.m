function [r1, r2] = ref1_load()
%% REFERENCE SIGNAL
% ********************************************************
freq = .5*pi;

r1 = @(x, t)[ % desired trajectory
    +1*cos(freq*t) + 1
    -1*cos(freq*t) - 1
    ] ;

r2 = @(x,t) [ 
    -1*sin(freq*t)*freq
    +1*sin(freq*t)*freq
];

end
