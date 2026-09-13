function [r, rd, rdd] = ref2_load(x, t)

dt = 1e-4;
r = @(x, t)[ % desired trajectory
    10
    10
    ] ;

rd = @(x,t) [ 
    0
    0
];

rdd = @(x, t) [
    0;0
];

end