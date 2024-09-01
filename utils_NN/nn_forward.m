function [nn, out, info] = nn_forward(nn, nnOpt, in)
    
    %% POINTER TO START POINT
    pt_V = 1;
    pt_tape = 1;

    %% FORWARD
    % save to gradient tape
    nn.tape(1: 1+nnOpt.NN_size(1)-1) = in;
    pt_tape = pt_tape + nnOpt.NN_size(1);

    for l_idx = 1:1:nnOpt.l_size-1
        n = nnOpt.NN_size(l_idx)+1;
        m = nnOpt.NN_size(l_idx+1);

        % get V weight
        V = reshape(nn.V(pt_V: pt_V+(n*m)-1), n,m);
    
        % phi_k = tanh(nn_out); 1]
        if l_idx == 1
            in = [in; 1]; % only one augmentation for 1st layer
        else
            in = [tanh(in); 1];
        end

        % forward
        in = V'*in;
        
        % save to gradient tape
        if l_idx~=nnOpt.l_size-1
            nn.tape(pt_tape: pt_tape+m-1) = in;
        end

        % pointer to next
        pt_V = pt_V + (n*m);
        pt_tape = pt_tape + m;
    end

    %% ERROR CHECK
    assert(pt_V-1 == nnOpt.v_size);
%     assert(pt_tape-1 == nnOpt.t_size);

    %% TERMINATION
    out = in;
    info = NaN;
    

end