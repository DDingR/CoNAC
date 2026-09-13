
x_in = [x1;x2;r];                           % NN input vector 

fil_x = x2 + opt.Lambda*x1;
fil_xd = xd2 + opt.Lambda*xd1;
fil_xdd = xd3 + opt.Lambda*xd2;

% x_in = [fil_xd; fil_xdd; fil_x];
x_in = [fil_x; fil_xd; xd3];

% fil_xd1: -5.309 ~ 4.000, fil_xd2: -23.577 ~ 15.723
% xd3_1: -1.176 ~ 1.176, xd3_2: -1.679 ~ 1.679

x_in_MIN = [
    -5.309; 
    -23.577; 
    -5.309; 
    -23.577; 
    -1.176; 
    -1.679
];
x_in_MAX = [
    4.000; 
    15.723; 
    4.000; 
    15.723; 
    1.176; 
    1.679
];
x_in = normalize_vertor(x_in, x_in_MIN, x_in_MAX);  % normalize the input vector

[nn, u_NN] = nnForward(nn, opt, x_in);      % NN forward propagation
u = u_NN;                           % control input (before saturation)   



function [x] = normalize_vertor(x, x_min, x_max)
    % Normalize the input vector x to the range [-1, 1]
    x = 2 * (x - x_min) ./ (x_max - x_min) - 1;
end