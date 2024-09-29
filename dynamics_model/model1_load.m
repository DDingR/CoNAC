function [grad_x, IC] = model1_load()
% The Cartesian Coordinates are ignored
% x = [q1 q2 \dot{q1} \dot{q2}]

%% INITIAL CONDITION
% ********************************************************
% IC.x = [1;-1;0;0];  
IC.x = [2;-2;0;0];  
IC.u = [0; 0];

%% SYSTEM DYNAMICS
% ********************************************************
grad_x = @(q, u, t) model_wrapper(q, u, t);

end

%% Local function
% ********************************************************
function f = model_wrapper(q, u, t)
    [M, C, G, F ] = model1(q, t);
    
    del_M = zeros(2,2); 
    % del_C = [
    %     tanh(q(1))*100
    %     -tanh(q(2))*100
    % ]; 
    del_C = zeros(2,2); 

    % del_G = rand(2,1)*1e0;
    del_G = zeros(2,1);

    M = M + del_M;
    C = C + del_C;
    G = G + del_G;

    % d = [
    %     tanh(q(1))*20 + sin(5*t)*3e1
    %     -tanh(q(2))*20 + cos(5*t)*3e1
    % ];
    d = zeros(2,1);

    f1 = q(3:4);
    f2 = M \ (-C*q(3:4) -G -F +u +d);
    
    f = [f1 ; f2];
end

% function x = forward_kinematics(q, L1, L2)

%     J = calculate_jacobian(q(1), q(2), L1, L2);
%     x = [L1*cos(q(1)) + L2*cos(q(1) + q(2));
%          L1*sin(q(1)) + L2*sin(q(1) + q(2));
%          J * q(3:4)];
% end

% function J = calculate_jacobian(q1, q2, L1, L2)

%     J = [(-L1*sin(q1) - L2*sin(q1 + q2)), -L2*sin(q1 + q2);
%           (L1*cos(q1) + L2*cos(q1 + q2)),  L2*cos(q1 + q2)];
% end



