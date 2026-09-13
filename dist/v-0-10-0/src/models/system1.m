function [grad_x, xd, IC] = system1()

%% SYSTEM DYNAMICS
% ********************************************************
grad_x = @(x, u) [
    x(2)
    % x(1)*x(2)*tanh(x(2))+sech(x(1)) + sum(u)
    -2*x(1)-3*x(2) + sum(u)
];

%% REFERENCE SIGNAL
% ********************************************************
xd = @(x, t) [
    sin(2*t) + sin(3*t)
    2*cos(2*t) + 3*cos(3*t)
];

%% INITIAL CONDITION
% ********************************************************
IC.x = [0; 0];
IC.u = [0; 0];

end