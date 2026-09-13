% -------------------------------
% not used, see loadGlobalOpts.m
% -------------------------------

function opt = loadOpts(dt)
    opt.dt = dt;

    opt.e_tol = 0e-3;
    opt.init_range = 1e-3;
    
    opt.Lambda = diag([1 1]) * 5e0;

    opt.alpha = 8e-1;
    opt.rho = opt.alpha*0e-2;
    opt.NN_size = [6,4,4,2];
    opt.W = 1;
    opt.e_size = 2;

    % opt.beta = [1 1 1] * 1e-1;
    opt.beta(1:3) = [1 1 1] * 0e-3;
    opt.beta(4) = 0e-2; % control ipnut
    opt.cstr.th_max = [11;12;13] * 1e0;
    opt.cstr.u_ball = 20;

    c_num = length(opt.beta);
    opt.lbd = zeros(c_num,1);
    
%% PASSIVE NUMBER
    % layer number
    opt.l_size = length(opt.NN_size);
    % total tape number
    opt.tp_size = sum(opt.NN_size(1:end-1));
    % total weight number (will be calc-ed)
    opt.th_size_list = zeros(opt.l_size-1 ,1);    
    for idx = 1:1:opt.l_size-1
        opt.th_size_list(idx) = (opt.NN_size(idx)+1) * opt.NN_size(idx+1);
    end

    opt.th_size = sum(opt.th_size_list);    

end
