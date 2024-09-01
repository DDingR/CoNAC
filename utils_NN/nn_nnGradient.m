function total_grad = nn_nnGradient(nn, nnOpt)

    %% PRE-ALLOCATION
    total_grad = zeros(nnOpt.NN_size(end), nnOpt.v_size);
    
    %% 
    pt_V = nnOpt.v_size;
    pt_tape = nnOpt.t_size;

    grad_to_back = 1;

    %% BACKWARD (TOTAL NN GRDIENT CALCULATION)
    for l_idx = flip(2:1:nnOpt.l_size)
        l2 = nnOpt.NN_size(l_idx);        
        l1 = nnOpt.NN_size(l_idx-1);        
        
        % pre_Phi
        Phi = nn.tape(pt_tape-l1+1: pt_tape);
        aug_phi = [tanh(Phi); 1];

        % (I \otimes phi_k^T)
        grad = kron(eye(l2), aug_phi');
        if l_idx ~= nnOpt.l_size
            n = l2+1; m = nnOpt.NN_size(l_idx+1);
            % get phi_dot
            phi_dot = tanh_dot(pre_Phi);
            % get V weight
            V = reshape(nn.V(pt_V-(n*m)+1: pt_V), n,m);
            grad_to_back = grad_to_back * V'*phi_dot;

            % pointer to next
            pt_V = pt_V - (n*m);
        end
        grad = grad_to_back*grad;

        % save to total grad
        total_grad(:, pt_V-((l1+1)*l2)+1: pt_V) = grad;

        % pointer to next
        pt_tape = pt_tape - l1;

        % for efficiency
        pre_Phi = Phi;

    end

    %% ERROR CHECK
    assert(pt_V == (nnOpt.NN_size(1)+1)*nnOpt.NN_size(2));
    assert(pt_tape == 0);
end

%% LOCAL FUNCTIONS
function phi_dot = tanh_dot(Phi)
    phi_dot = zeros(length(Phi)+1, length(Phi));
    for idx = 1:1:length(Phi)
        phi_dot(idx,idx) = 1-tanh(Phi(idx))^2;
    end
end