function grad_x = model1_load()
% The Cartesian Coordinates are ignored
% x = [q1 q2 \dot{q1} \dot{q2}]

%% SYSTEM DYNAMICS
% ********************************************************
grad_x = @(q, u, t) model_wrapper(q, u, t);

end

%% Local function


% ********************************************************
function f = model_wrapper(q, u, t)
    % [M, C, G, F ] = model1(q, u, t);

    dof = 2;
    p = config_robot(dof);
    [M, C, G, F] = robot_dyn(q(1:2), q(3:4), p, dof);
    
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



