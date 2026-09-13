function [nn, out] = nnForward(nn, opt, in)
    
    %% POINTER TO START POINT
    pt_V = 1;
    pt_tape = 1;

    %% FORWARD
    % save to gradient tape
    nn.tape(1: 1+opt.NN_size(1)-1) = in;
    pt_tape = pt_tape + opt.NN_size(1);

    for l_idx = 1:1:opt.l_size-1
        n = opt.NN_size(l_idx)+1;
        m = opt.NN_size(l_idx+1);

        % get V weight
        V = reshape(nn.th(pt_V: pt_V+(n*m)-1), n,m);
    
        % phi_k = tanh(nn_out); 1]
        if l_idx == 1
            in = [in; 1]; % only one augmentation for 1st layer
        else
            in = [tanh(in); 1];
        end

        % forward
        in = V'*in;
        
        % save to gradient tape
        if l_idx~=opt.l_size-1
            nn.tape(pt_tape: pt_tape+m-1) = in;
        end

        % pointer to next
        pt_V = pt_V + (n*m);
        pt_tape = pt_tape + m;
    end

    %% ERROR CHECK
    assert(pt_V-1 == opt.th_size);
%     assert(pt_tape-1 == opt.t_size);

    %% TERMINATION
    out = in;    

end