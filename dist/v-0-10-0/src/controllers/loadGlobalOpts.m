function opt = loadGlobalOpts(dt, CONTROL_NUM, OPT_NUM)
    opt.dt = dt;                        % control sampling time

    % opt.Lambda = diag([1 5]) * .1;      % filtered error gain
    opt.Lambda = diag([5 15]);      % filtered error gain
    % r = e2 + Lambda * e1
    
    opt.init_range = 1e0;              % initial NN weight range
    opt.NN_size = [6,4,4,2];            % NN layer size (input, hiddens, output)
                                        %   e.g. [6,4,4,2] 
                                        %   ->  6 input nodes, 
                                        %       2 hidden layers with 4 nodes each, 
                                        %       2 output nodes
                                        %   -> W0, W1, W2 weight matrices with sizes
                                        %       (6+1)x4, (4+1)x4, (4+1)x2 (with bias)

    % opt.alpha = 15;                     % learning rate (for NNs)
    opt.alpha = .2;                     % learning rate (for NNs)
    opt.e_size = 2;                     % error size (number of error signals)

    opt.e_tol = 0e-3;                   % error tolerance (not used)
    opt.rho = opt.alpha*0e-2;           % sigma-modification gain (not used)
    opt.W = diag([1 1]);                % feedback error weighting matrix (not used)

    opt.cstr.th_max = [6;6;6] * 1e1; % NN weight L2 norm constraints 
    opt.cstr.u_ball = 11;             % control input L2 norm constraint (ball constraint)   
    % opt.cstr.uMax1 = 10;                 % control input 2 max constraint (box constraint)
    opt.cstr.uMax2 = 4;                 % control input 2 max constraint (box constraint)
    opt.cstr.uMax1 = sqrt(opt.cstr.u_ball^2 - opt.cstr.uMax2^2);
    % opt.cstr.uMax2 = opt.cstr.u_ball/sqrt(2);
    % opt.cstr.uMax1 = opt.cstr.u_ball/sqrt(2);
                                        % control input 1 max constraint 
                                        %       (box constraint, calculated from ball constraint and control input 2 max constraint)

    if CONTROL_NUM == 1
        % beta is Lagrange multiplier update gain for each constraint;
        %   \dot{lambda} = beta * c, where c is the constraint violation
        opt.beta(1:3) = [1 1 1] * 1e0; % NN weight constraints
        opt.beta(4) = 1e0;              % control input ball
        opt.beta(5) = 0e2;              % control input 1 Max
        opt.beta(6) = 1e2;              % control input 2 Max
        opt.beta(7) = opt.beta(5);      % control input 1 Min
        opt.beta(8) = opt.beta(6);      % control input 2 Min
        % opt.beta = opt.beta/opt.alpha;

        c_num = length(opt.beta);
        opt.lbd = zeros(c_num,1);

    elseif CONTROL_NUM == 2
        % for auxiliary control
        %   \dot{z} = A_zeta * z + B_zeta * \Delta u,
        %   where auxiliary state z and input saturation \Delta u.
        opt.A_zeta = -10 * eye(2);     % stable matrix
        opt.B_zeta = 10000 * eye(2);    % input gain matrix

    elseif CONTROL_NUM == 3
        % for complex auxiliary control
        opt.alp_prev = zeros(2,1);
        opt.fil_alp = zeros(2,1);

    end
    
    %% PASSIVE NUMBER
    opt.l_size = length(opt.NN_size);           % layer number
    opt.tp_size = sum(opt.NN_size(1:end-1));    % total tape number
    opt.th_size_list = zeros(opt.l_size-1 ,1);  % total weight number (will be calc-ed) 
    for idx = 1:1:opt.l_size-1
        % for each layer, weight size = (input size + 1) * output size (with bias)
        opt.th_size_list(idx) = (opt.NN_size(idx)+1) * opt.NN_size(idx+1);
    end
    opt.th_size = sum(opt.th_size_list);        % total weight number

    %% OPTIONAL: OVERWRITE OPT WITH VAR_OPTS
    switch CONTROL_NUM
        case 1
            % CoNAC
            switch OPT_NUM
                case 1
                    opt.beta(4) = opt.beta(4) * 100;
                    opt.beta(6) = opt.beta(6) * 100;
                    opt.beta(8) = opt.beta(8) * 100;
                case 2
                   
                case 3
                    opt.beta(4:end) = 0;
        case 2
            % Auxiliary Control
        case 3
            % Complex Auxiliary Control
    end


end
