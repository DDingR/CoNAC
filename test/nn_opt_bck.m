function nnOpt = nn_opt(ctrl_name, dt)
    %% DEFAULT PARAMETERS
    nnOpt.dt = dt; % sampling time

    nnOpt.e_tol = 0e-3; % error tolerance

    nnOpt.init_range = 1e-3; % weight init range

    % nnOpt.W = eye(4);
    nnOpt.W = diag([1 1 15 15]);

    %% CONTROLLER PARAMETERS
    nnOpt.ctrl_param.k1 = 1.1;
    nnOpt.ctrl_param.k2 = 10;
    
    nnOpt.ctrl_param.Azeta = diag([1 1]) * 1e1 * 2.0;
    nnOpt.ctrl_param.Bzeta = diag([1 1]) * 1e0 * 10;
    % nnOpt.ctrl_param.delta = 1e-1;
    nnOpt.ctrl_param.delta = 1e3;

    nnOpt.ctrl_param.M = eye(2) * 1;
    nnOpt.ctrl_param.C = eye(2) * 1;
    nnOpt.ctrl_param.G = zeros(2,1);
    nnOpt.ctrl_param.inv_M = inv(nnOpt.ctrl_param.M);

    %% ALGORITHM PARAMETERS
    if strcmp(ctrl_name, "CTRL2")
        nnOpt.alg = "Dixon";

        nnOpt.alpha = 2e3; % learning rate
        % nnOpt.rho = nnOpt.alpha * 1e-2; % e-modification
        nnOpt.rho = 0;      
        
    elseif strcmp(ctrl_name, "CTRL3")
        nnOpt.alg = "Dixon";

        nnOpt.alpha = 2e3; % learning rate
        % nnOpt.rho = nnOpt.alpha * 1e-2; % e-modification
        nnOpt.rho = 0;

%     elseif strcmp(ctrl_name, "Kasra")
%         error("not use now")

    elseif strcmp(ctrl_name, "CTRL4")
        nnOpt.alg = "Proposed";

        nnOpt.alpha = 1e3; % learning rate

        c_num = 8;

        nnOpt.Lambda = zeros(c_num,1);
%         nnOpt.Beta = ones(c_num,1) * 1e+0;
        nnOpt.Beta = [ ...
            1 1 1 ...   % weight ball
            0 0 ...     % maximum input
            0 0 ...     % minimum input
            1 ...       % input ball
            ] * 1e0;

    elseif strcmp(ctrl_name, "CTRL1")
        nnOpt.alg = "Backstepping";

    elseif strcmp(ctrl_name, "CTRL5")
        nnOpt.alg = "Dixon";

        nnOpt.alpha = 1e3; % learning rate
        % nnOpt.rho = nnOpt.alpha * 1e-2; % e-modification
        nnOpt.rho = 0;

    elseif strcmp(ctrl_name, "ALM")
        error("not ready")
    end

    nnOpt.cstr.V_max = [
        20
        30
        50
    ];
    
    nnOpt.cstr.U_max = [
        40
        2e2
    ] * 1e0;
    
    nnOpt.cstr.U_min = [
        1e2
        5e1
    ] * -1e2;

    nnOpt.cstr.u_ball = 5e1;
    
    %% NN STRUCTURE
    nnOpt.NN_size = [
        % layer info (node numbers)
        2 % input layer [r; rd]
        8
        8
        2 % output layer
    ]; 

    %% NUMBERS
    % layer number
    nnOpt.l_size = length(nnOpt.NN_size);
    % total tape number
    nnOpt.t_size = sum(nnOpt.NN_size(1:end-1));
    % total weight number (will be calc-ed)
    nnOpt.v_size_list = zeros(nnOpt.l_size-1 ,1);    
    for idx = 1:1:nnOpt.l_size-1
        nnOpt.v_size_list(idx) = (nnOpt.NN_size(idx)+1) * nnOpt.NN_size(idx+1);
                  
    end

    nnOpt.v_size = sum(nnOpt.v_size_list);    

    %% ETC
    
end