function [c, cd] = nn_cstr(nn, nnOpt, u_NN, nnGrad)
    %% CONSTRAINTS PARAMETERS
    V_max = nnOpt.cstr.V_max;
    U_max = nnOpt.cstr.U_max;
    U_min = nnOpt.cstr.U_min;
    u_ball = nnOpt.cstr.u_ball;
    
    %% ERROR CHECK
    assert(length(V_max) == nnOpt.l_size-1)
    
    %% ACTIVE SET CHECK
    % ball condition
    c_b = zeros(nnOpt.l_size-1, 1);
    
    cumsum_V = [0;cumsum(nnOpt.v_size_list)];
    for l_idx = 1:1:nnOpt.l_size-1
        start_pt = cumsum_V(l_idx)+1;
        end_pt = cumsum_V(l_idx+1);   

        c_b(l_idx) = norm(nn.V(start_pt: end_pt))^2 ...
            - V_max(l_idx)^2;
    end
    
    % input saturation 
    u = -u_NN;
    c_M = u - U_max;
    c_m = U_min - u;
    
    % input ball
    c_ub = norm(u)^2 - u_ball^2;

    % summary
    c = [c_b; c_M; c_m; c_ub];
    
    %% ACTIVE SET CONSTRAINT GRADIENT CALC
    cd = zeros(length(c), nnOpt.v_size);
    for l_idx = 1:1:length(V_max)
        start_pt = cumsum_V(l_idx)+1;
        end_pt = cumsum_V(l_idx+1);
        cd(l_idx, start_pt:end_pt) = 2 * nn.V(start_pt:end_pt,1);
    end
    
    cd(l_idx+1:end, :) = [
        -nnGrad;
        +nnGrad;
        2*u(1)*-nnGrad(1,:) + 2*u(2)*-nnGrad(2,:)
        ];
end









