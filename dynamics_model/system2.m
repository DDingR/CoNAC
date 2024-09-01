function [grad_x, r, rd, rdd, IC] = system2()

% x = [q1 q2 \dot{q1} \dot{q2}]

%% INITIAL CONDITION
% ********************************************************
IC.x = [1;-1;0;0];  
IC.u = [0; 0];

%% SYSTEM DYNAMICS
% ********************************************************

env.v = 13.33;
env.q01 = 8.98;
env.q02 = 8.75;
env.g = 9.8;
env.alpha = 3;
env.gamma = 20;
env.k1 = 0.001;

M = @(x) [
    env.v+env.q01+2*env.gamma*cos(x(2))     env.q01+env.q02*cos(x(2))
    env.q01+env.q02*cos(x(2))               env.q01
    ];

C = @(x) [
    -env.q02*x(4)*sin(x(2))                 -env.q02*(x(3)+x(4))*sin(x(2))
    env.q02*x(3)*sin(x(2))                  0
    ];
    

G = @(x) [
    15*env.g*cos(x(1))+8.75*env.g*cos(x(1)+x(2))
    8.75*env.g*cos(x(1)+x(2))
    ];

% d = @(x) zeros(2,1);
d = @(x, t) [
    tanh(x(1))*20 + sin(5*t)*3e1
    -tanh(x(2))*20 + cos(5*t)*3e1
];

grad_x = @(x,u, t) blkdiag(eye(2), M(x)) \ ...
    ([zeros(2) eye(2); zeros(2) -C(x)] * x + [zeros(2);eye(2)] * (-G(x)+u+d(x, t)));

%% REFERENCE SIGNAL
% ********************************************************


% dt = 1e-4;
% r = @(x, t)[ % desired trajectory
%     heaviside(t - 3) + 1
%     heaviside(t - 3) - 1
%     ] ;
% 
% rd = @(x,t) [ 
%     1/dt * (heaviside(t - 3) - heaviside(t - 3-dt))
%     1/dt * (heaviside(t - 3) - heaviside(t - 3-dt))
% ];
% 
% rdd = @(x, t) [
%     0;0
% ];

r = @(x, t)[ % desired trajectory
    1+1*sin(0.5*pi*t)
    -1-1*sin(0.5*pi*t)
    ] ;

rd = @(x,t) [ 
    1*cos(0.5*pi*t)*0.5*pi
    -1*cos(0.5*pi*t)*0.5*pi
];

rdd = @(x, t) [
    -1*sin(0.5*pi*t)*(0.5*pi)^2
    +1*sin(0.5*pi*t)*(0.5*pi)^2
];

end