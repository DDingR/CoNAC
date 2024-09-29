function [nn, nnOpt, dot_L, info] = nn_backward(nn, nnOpt, ctrlOpt, e, u_NN)

    %% NN GRADIENT CALCULATION
    nnGrad = nn_nnGradient(nn, nnOpt);

    % dead-zone
    if norm(e) < nnOpt.e_tol
        info = NaN;
        % dot_L = norm(V_grad);
        dot_L = norm(0);
        return;
    end

    %% UPDATE
    k1 = ctrlOpt.k1; k2 = ctrlOpt.k2;
    inv_M = ctrlOpt.inv_M;
    % M = ctrlOpt.M;
    A = [-k1*eye(2) eye(2); -eye(2) -k2*eye(2)];
    B = [zeros(2); inv_M];

    if strcmp(nnOpt.alg, "Dixon")
        % V_grad =  - nnOpt.alpha * -nnGrad' * B' * e;
        V_grad =  - nnOpt.alpha * -nnGrad' * inv_M * e(3:4);
        V_grad = V_grad - nnOpt.rho * norm(e(3:4)) * nn.V;
        V_grad = V_grad * nnOpt.dt;

        [c, cd] = nn_cstr(nn, nnOpt, u_NN, nnGrad);

        ActSet  = double(c>0);
        c       = double(c .* ActSet);
        cd      = double(cd .* ActSet);

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
        dhdu_11 = dfdx_11(u_NN(1), u_NN(2));
        dhdu_12 = dfdx_12(u_NN(1), u_NN(2));
        dhdu_21 = dfdx_21(u_NN(1), u_NN(2));
        dhdu_22 = dfdx_22(u_NN(1), u_NN(2));
        dhdu = [dhdu_11 dhdu_12; dhdu_21 dhdu_22];

        nn.eta = nn.eta + (A*nn.eta+B*dhdu*-nnGrad) * nnOpt.dt;

        % active set check
        [c, cd] = nn_cstr(nn, nnOpt, u_NN, nnGrad);

        ActSet  = double(c>0);
        Lambda  = nnOpt.Lambda .* ActSet;
        % Lambda  = nnOpt.Lambda;
        c       = double(c .* ActSet);
        cd      = double(cd .* ActSet);

        % find gradient; theta, lambda
        V_grad = - nnOpt.alpha * (nn.eta'*nnOpt.W*e + cd' * Lambda);
        % V_grad = - nnOpt.alpha * ((A\B*nnGrad)'*e + cd' * Lambda);
        % V_grad = - nnOpt.alpha * (-nnGrad' *e(3:4)+ cd' * Lambda);
        % V_grad = - nnOpt.alpha * (-nnGrad' *nnOpt.W(3:4,3:4)*e(3:4)+ cd' * Lambda);
        L_grad = diag(nnOpt.beta) * c;
        % L_grad = -1e2 * Lambda.*~ActSet + L_grad;
 
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

    nn.V = nn.V + V_grad;

    %% TERMINATION
    dot_L = norm(V_grad);
    info = NaN;
    
end