function nn = nn_init(nnOpt)

    %% MATRICES PRE-ALLOCATION
    nn.V = (rand(nnOpt.v_size,1)-1/2)*2*nnOpt.init_range;
    nn.tape = zeros(nnOpt.t_size, 1);
    
    %% SENSITIVY INITIALIZATION (PROPOSED)
    if strcmp(nnOpt.alg, "Proposed")
        nn.eta = zeros(nnOpt.NN_size(end)*2, nnOpt.v_size); % dPhi/dth
        % backsteping check
%         nn.eta = zeros(nnOpt.NN_size(end), nnOpt.v_size); % dPhi/dth
    end

end