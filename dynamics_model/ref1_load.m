function [r, rd, rdd] = ref1_load()
%% REFERENCE SIGNAL
% ********************************************************

r = @(x, t)[ % desired trajectory
    +1*cos(0.5*pi*t) + 1
    -1*cos(0.5*pi*t) - 1
    ] ;

rd = @(x,t) [ 
    -1*sin(0.5*pi*t)*0.5*pi
    +1*sin(0.5*pi*t)*0.5*pi
];

rdd = @(x, t) [
    -1*cos(0.5*pi*t)*(0.5*pi)^2
    +1*cos(0.5*pi*t)*(0.5*pi)^2
];

end
