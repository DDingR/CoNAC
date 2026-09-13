function [c, cd] = nnCstr(nn, opt, nnGrad)
    %% CONSTRAINTS PARAMETERS
    th_max = opt.cstr.th_max;
    
    %% ERROR CHECK
    assert(length(th_max) == opt.l_size-1)
    
    %% ACTIVE SET CHECK
    % ball condition
    c = zeros(opt.l_size-1, 1);
    
    cumsum_th = [0;cumsum(opt.th_size_list)];
    for l_idx = 1:1:opt.l_size-1
        start_pt = cumsum_th(l_idx)+1;
        end_pt = cumsum_th(l_idx+1);   

        c(l_idx) = norm(nn.th(start_pt: end_pt))^2 ...
            - th_max(l_idx)^2;
    end
    
    %% ACTIVE SET CONSTRAINT GRADIENT CALC
    for l_idx = 1:1:length(th_max)
        start_pt = cumsum_th(l_idx)+1;
        end_pt = cumsum_th(l_idx+1);
        cd(l_idx, start_pt:end_pt) = 2 * nn.th(start_pt:end_pt,1);
    end
    
end









