function [nn, opt] = nnBackward(nn, opt, e, u_NN)

    %% NN GRADIENT CALCULATION
    nnGrad = nnGradient(nn, opt);

    % dead-zone
    if norm(e) < opt.e_tol
        return;
    end

    %% UPDATE
    % active set check
    [c, cd] = nnCstr(nn, opt, nnGrad);

    % find gradient; theta, lambda
    th_grad = - opt.alpha * (nnGrad'*opt.W*e);
    th_grad = th_grad + - opt.rho * nn.th;

    th_grad = th_grad * opt.dt;
    % update
    nn.th = nn.th + th_grad;

    %% TERMINATION
    
end