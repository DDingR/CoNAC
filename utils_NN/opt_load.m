function [ctrlOpt, nnOpt] = opt_load(c_idx, opts, dt)

    type_idx = find(opts.Properties.RowNames == "type");
    alg_idx = find(opts.Properties.RowNames == "alg");
    e_tol_idx = find(opts.Properties.RowNames == "e_tol");
    init_range_idx = find(opts.Properties.RowNames == "init_range");
    k1_idx = find(opts.Properties.RowNames == "k1");
    k2_idx = find(opts.Properties.RowNames == "k2");
    M_idx = find(opts.Properties.RowNames == "M");
    C_idx = find(opts.Properties.RowNames == "C");  
    G_idx = find(opts.Properties.RowNames == "G");
    NN_size_idx = find(opts.Properties.RowNames == "NN_size");
    alpha_idx = find(opts.Properties.RowNames == "alpha");
    rho_idx = find(opts.Properties.RowNames == "rho");
    W_idx = find(opts.Properties.RowNames == "W");
    beta_idx = find(opts.Properties.RowNames == "beta");
    V_max_idx = find(opts.Properties.RowNames == "V_max");
    U_max_idx = find(opts.Properties.RowNames == "U_max");
    U_min_idx = find(opts.Properties.RowNames == "U_min");
    u_ball_idx = find(opts.Properties.RowNames == "u_ball");
    Azeta_idx = find(opts.Properties.RowNames == "Azeta");
    Bzeta_idx = find(opts.Properties.RowNames == "Bzeta");
    L2_idx = find(opts.Properties.RowNames == "L2");
 
    %% DEFAULT PARAMETERS
    nnOpt.dt = dt; % sampling time

    type = opts{type_idx, c_idx};

    if strcmp(type, "CTRL1") % BSC
        ctrlOpt.type = type;
        nnOpt.alg = opts{alg_idx, c_idx};

        ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
        ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});
    
        ctrlOpt.M = str2num(opts{M_idx, c_idx});
        ctrlOpt.C = str2num(opts{C_idx, c_idx});
        ctrlOpt.G = str2num(opts{G_idx, c_idx});
        ctrlOpt.inv_M = inv(ctrlOpt.M);

        ctrlOpt.u_ball = str2num(opts{u_ball_idx, c_idx});

        % nnOpt.v_size = 1;
        % nnOpt.init_range = 1;
        % nnOpt.NN_size = [2 2];
        % nnOpt.t_size = 1;
        nnOpt.l_size = 2;
        % nnOpt.v_size_list = [2 2];

        nnOpt.cstr.V_max = str2num(opts{V_max_idx, c_idx});
        nnOpt.cstr.u_ball = str2num(opts{u_ball_idx, c_idx});
        nnOpt.cstr.U_max = str2num(opts{U_max_idx, c_idx});
        nnOpt.cstr.U_min = str2num(opts{U_min_idx, c_idx});
    
    elseif strcmp(type, "CTRL2") % BSC + Dixon
        ctrlOpt.type = type;

        nnOpt.alg = opts{alg_idx, c_idx};

        nnOpt.e_tol = str2num(opts{e_tol_idx, c_idx});
        nnOpt.init_range = str2num(opts{init_range_idx, c_idx});
    
        ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
        ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});
    
        ctrlOpt.M = str2num(opts{M_idx, c_idx});
        ctrlOpt.C = str2num(opts{C_idx, c_idx});
        ctrlOpt.G = str2num(opts{G_idx, c_idx});
        ctrlOpt.inv_M = inv(ctrlOpt.M);

        ctrlOpt.u_ball = str2num(opts{u_ball_idx, c_idx});

        nnOpt.alpha = str2num(opts{alpha_idx, c_idx});
        nnOpt.rho = str2num(opts{rho_idx, c_idx});
        nnOpt.NN_size = str2num(opts{NN_size_idx, c_idx});

        nnOpt.cstr.V_max = str2num(opts{V_max_idx, c_idx});   
        nnOpt.cstr.u_ball = str2num(opts{u_ball_idx, c_idx});
        nnOpt.cstr.U_max = str2num(opts{U_max_idx, c_idx});
        nnOpt.cstr.U_min = str2num(opts{U_min_idx, c_idx});

        nnOpt = numbersCalc(nnOpt);
    
    elseif strcmp(type, "CTRL3") % Dixon + Aux. State
        ctrlOpt.type = type;

        nnOpt.alg = opts{alg_idx, c_idx};

        nnOpt.e_tol = str2num(opts{e_tol_idx, c_idx});
        nnOpt.init_range = str2num(opts{init_range_idx, c_idx});
    
        ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
        ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});
    
        ctrlOpt.M = str2num(opts{M_idx, c_idx});
        ctrlOpt.C = str2num(opts{C_idx, c_idx});
        ctrlOpt.G = str2num(opts{G_idx, c_idx});
        ctrlOpt.inv_M = inv(ctrlOpt.M);

        ctrlOpt.u_ball = str2num(opts{u_ball_idx, c_idx});

        nnOpt.alpha = str2num(opts{alpha_idx, c_idx});
        nnOpt.rho = str2num(opts{rho_idx, c_idx});
        nnOpt.NN_size = str2num(opts{NN_size_idx, c_idx});

        nnOpt.cstr.V_max = str2num(opts{V_max_idx, c_idx});
        nnOpt.cstr.u_ball = str2num(opts{u_ball_idx, c_idx});
        nnOpt.cstr.U_max = str2num(opts{U_max_idx, c_idx});
        nnOpt.cstr.U_min = str2num(opts{U_min_idx, c_idx});

        ctrlOpt.Azeta = str2num(opts{Azeta_idx, c_idx});
        ctrlOpt.Bzeta = str2num(opts{Bzeta_idx, c_idx});

        nnOpt = numbersCalc(nnOpt);
    
    elseif strcmp(type, "CTRL4") % Proposed
        ctrlOpt.type = type;

        nnOpt.alg = opts{alg_idx, c_idx};

        nnOpt.e_tol = str2num(opts{e_tol_idx, c_idx});
        nnOpt.init_range = str2num(opts{init_range_idx, c_idx});
        
        ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
        ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});
    
        ctrlOpt.M = str2num(opts{M_idx, c_idx});
        ctrlOpt.inv_M = inv(ctrlOpt.M);

        ctrlOpt.u_ball = str2num(opts{u_ball_idx, c_idx});

        nnOpt.alpha = str2num(opts{alpha_idx, c_idx});
        nnOpt.NN_size = str2num(opts{NN_size_idx, c_idx});
        nnOpt.W = str2num(opts{W_idx, c_idx});
        nnOpt.beta = str2num(opts{beta_idx, c_idx});
        nnOpt.cstr.V_max = str2num(opts{V_max_idx, c_idx});
        nnOpt.cstr.u_ball = str2num(opts{u_ball_idx, c_idx});
        nnOpt.cstr.U_max = str2num(opts{U_max_idx, c_idx});
        nnOpt.cstr.U_min = str2num(opts{U_min_idx, c_idx});

        c_num = length(nnOpt.beta);
        nnOpt.Lambda = zeros(c_num,1);
        
        nnOpt.L2 = str2num(opts{L2_idx, c_idx});

        nnOpt = numbersCalc(nnOpt);

    end
    
end

%% LOCAL FUNCTIONS
function nnOpt = numbersCalc(nnOpt)
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

end