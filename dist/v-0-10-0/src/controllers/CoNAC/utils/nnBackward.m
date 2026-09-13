function [nn, opt] = nnBackward(nn, opt, e, u_NN)

    %% NN GRADIENT CALCULATION
    nnGrad = nnGradient(nn, opt);

    % dead-zone
    if norm(e) < opt.e_tol
        return;
    end

    %% BLF??
    % max_x = 2;
    % e(1) = e(1)/(max_x^2 - e(1)^2);

    %% UPDATE
    % active set check
    [c, cd] = nnCstr(nn, opt, u_NN, nnGrad);

    ActSet  = double(c>0);
    ActSet = ones(size(ActSet)); % ignore constraints for now
    % ActSet = zeros(size(ActSet)); % ignore constraints for now
    lbd  = opt.lbd .* ActSet;

    % find gradient; theta, lambda
    th_grad = - opt.alpha * (nnGrad'*opt.W*e + cd'*lbd);
    th_grad = th_grad + - opt.rho * nn.th;
    lbd_grad = diag(opt.beta) * c;

    th_grad = th_grad * opt.dt;
    lbd_grad = lbd_grad * opt.dt;

    % update
    nn.th = nn.th + th_grad;
    opt.lbd = lbd + lbd_grad;
    opt.lbd = max(opt.lbd, 0);

    %% TERMINATION
    
end