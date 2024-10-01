function nn = nn_init(nnOpt)

    if strcmp(nnOpt.alg, "none")
        nn.V = 0;
        return
    end

    %% MATRICES PRE-ALLOCATION
    nn.V = (rand(nnOpt.v_size,1)-1/2)*2*nnOpt.init_range;
    nn.tape = zeros(nnOpt.t_size, 1);
    
    %% SENSITIVY INITIALIZATION (PROPOSED)
    if strcmp(nnOpt.alg, "Proposed") || strcmp(nnOpt.alg, "L2")
        nn.eta = zeros(nnOpt.NN_size(end)*2, nnOpt.v_size);
    end

end