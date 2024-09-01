function [nn, nnOpt, dot_L, info] = nn_backward(nn, nnOpt, e, u_NN)

    %% NN GRADIENT CALCULATION
    nnGrad = nn_nnGradient(nn, nnOpt);

    %% UPDATE
    k1 = nnOpt.ctrl_param.k1; k2 = nnOpt.ctrl_param.k2;
    % M = nnOpt.ctrl_param.M;
    A = [-k1*eye(2) eye(2); -eye(2) -k2*eye(2)];
    B = [zeros(2); eye(2)];

    if strcmp(nnOpt.alg, "Dixon")

        [c, cd] = nn_cstr(nn, nnOpt, u_NN, nnGrad);

        ActSet  = double(c>0);
        c       = double(c .* ActSet);
        cd      = double(cd .* ActSet);

        % V_grad =  - nnOpt.alpha * -nnGrad' * B' * e;
        V_grad =  - nnOpt.alpha * -nnGrad' * e(1:2);
        V_grad = V_grad - nnOpt.rho * norm(e) * nn.V;
        V_grad = V_grad * nnOpt.dt;

        act_chk = double(cd*V_grad > 0) & c;

        for idx = 1:1:nnOpt.l_size-1
            
            cumsum_V = [0;cumsum(nnOpt.v_size_list)];
            if act_chk(idx)
                start_pt = cumsum_V(idx)+1;
                end_pt = cumsum_V(idx+1);       

                cd_focus = cd(idx, (start_pt:end_pt))';

                V_grad(start_pt:end_pt) = (eye(end_pt-start_pt+1) - ...
                    nnOpt.alpha * (cd_focus*cd_focus') / ...
                    (cd_focus'*nnOpt.alpha*cd_focus)) * V_grad(start_pt:end_pt);
            end
        end

    elseif strcmp(nnOpt.alg, "Kasra")
        V_grad = - nnOpt.alpha * (A\B * nnGrad)'  * e;
        V_grad = V_grad - nnOpt.rho * norm(e) * nn.V;
        V_grad = V_grad * nnOpt.dt;

    elseif strcmp(nnOpt.alg, "Proposed")
        % dynamical gradient
        nn.eta = nn.eta + (A*nn.eta+B*-nnGrad) * nnOpt.dt;

        % active set check
        [c, cd] = nn_cstr(nn, nnOpt, u_NN, nnGrad);

        ActSet  = double(c>0);
        Lambda  = nnOpt.Lambda .* ActSet;
        c       = double(c .* ActSet);
        cd      = double(cd .* ActSet);

        % find gradient; theta, lambda
        V_grad = - nnOpt.alpha * (nn.eta'*nnOpt.W*e + cd' * Lambda);
        L_grad = diag(nnOpt.Beta) * c;

        V_grad = V_grad * nnOpt.dt;
        L_grad = L_grad * nnOpt.dt;

        % update mupliers
        nnOpt.Lambda = Lambda + L_grad;
        nnOpt.Lambda = max(nnOpt.Lambda, 0);

    elseif strcmp(nnOpt.alg, "ALM")
        error("not ready")

    else
        error("wrong NN learning algorithm; check in [nn_opt]")

    end

    % dead-zone
    if norm(e) > nnOpt.e_tol
        nn.V = nn.V + V_grad;
    else
        info = NaN;
        dot_L = norm(V_grad);
        return;
    end

    %% TERMINATION
    dot_L = norm(V_grad);
    info = NaN;
    
end