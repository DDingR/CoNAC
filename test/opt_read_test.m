clear

opts = readtable("ctrlOpts.csv",'delimiter', ',', 'TextType', 'string','ReadVariableNames',0, 'ReadRowNames', 1);

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
beta_idx = find(opts.Properties.RowNames == "beta");
V_max_idx = find(opts.Properties.RowNames == "V_max");
U_max_idx = find(opts.Properties.RowNames == "U_max");
U_min_idx = find(opts.Properties.RowNames == "U_min");
u_ball_idx = find(opts.Properties.RowNames == "u_ball");
Azeta_idx = find(opts.Properties.RowNames == "Azeta");
Bzeta_idx = find(opts.Properties.RowNames == "Bzeta");

c_idx = 1;

type = opts{type_idx, c_idx};

if strcmp(type, "CTRL1") % BSC
    ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
    ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});

    ctrlOpt.M = str2num(opts{M_idx, c_idx});
    ctrlOpt.C = str2num(opts{C_idx, c_idx});
    ctrlOpt.G = str2num(opts{G_idx, c_idx});

elseif strcmp(type, "CTRL2") % BSC + Dixon
    nnOpt.e_tol = str2num(opts{e_tol_idx, c_idx});
    nnOpt.init_range = str2num(opts{init_range_idx, c_idx});

    ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
    ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});

    ctrlOpt.M = str2num(opts{M_idx, c_idx});
    ctrlOpt.C = str2num(opts{C_idx, c_idx});
    ctrlOpt.G = str2num(opts{G_idx, c_idx});
    
    nnOpt.alpha = str2num(opts{alpha_idx, c_idx});
    nnOpt.rho = str2num(opts{rho_idx, c_idx});

    nnOpt.cstr.V_max = str2num(opts{V_max_idx, c_idx});   
    nnOpt.cstr.U_max = str2num(opts{U_max_idx, c_idx});
    nnOpt.cstr.U_min = str2num(opts{U_min_idx, c_idx});

elseif strcmp(type, "CTRL3") % Dixon + Aux. State
    nnOpt.e_tol = str2num(opts{e_tol_idx, c_idx});
    nnOpt.init_range = str2num(opts{init_range_idx, c_idx});

    ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
    ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});

    ctrlOpt.M = str2num(opts{M_idx, c_idx});
    ctrlOpt.C = str2num(opts{C_idx, c_idx});
    ctrlOpt.G = str2num(opts{G_idx, c_idx});
    
    nnOpt.alpha = str2num(opts{alpha_idx, c_idx});
    nnOpt.rho = str2num(opts{rho_idx, c_idx});

    nnOpt.cstr.V_max = str2num(opts{V_max_idx, c_idx});
    nnOpt.cstr.u_ball = str2num(opts{u_ball_idx, c_idx});
    nnOpt.cstr.U_max = str2num(opts{U_max_idx, c_idx});
    nnOpt.cstr.U_min = str2num(opts{U_min_idx, c_idx});

    ctrlOpt.A_zeta = str2num(opts{Azeta_idx, c_idx});
    ctrlOpt.B_zeta = str2num(opts{Bzeta_idx, c_idx});

elseif strcmp(type, "CTRL4") % Proposed
    nnOpt.e_tol = str2num(opts{e_tol_idx, c_idx});
    nnOpt.init_range = str2num(opts{init_range_idx, c_idx});
    
    ctrlOpt.k1 = str2num(opts{k1_idx, c_idx});
    ctrlOpt.k2 = str2num(opts{k2_idx, c_idx});

    ctrlOpt.M = str2num(opts{M_idx, c_idx});
    
    nnOpt.alpha = str2num(opts{alpha_idx, c_idx});

    nnOpt.beta = str2num(opts{beta_idx, c_idx});
    nnOpt.cstr.V_max = str2num(opts{V_max_idx, c_idx});
    nnOpt.cstr.u_ball = str2num(opts{u_ball_idx, c_idx});
    nnOpt.cstr.U_max = str2num(opts{U_max_idx, c_idx});
    nnOpt.cstr.U_min = str2num(opts{U_min_idx, c_idx});
end
