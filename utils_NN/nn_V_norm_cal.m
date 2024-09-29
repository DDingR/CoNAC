function Vn = nn_V_norm_cal(V, nnOpt)

    if strcmp(nnOpt.alg, "none")
        Vn = 0;
        return
    end

    Vn = zeros(nnOpt.l_size-1,1);

    cumsum_V = [0;cumsum(nnOpt.v_size_list)];
    for l_idx = 1:1:nnOpt.l_size-1
        start_pt = cumsum_V(l_idx)+1;
        end_pt = cumsum_V(l_idx+1);   

        Vn(l_idx) = norm(V(start_pt: end_pt));
    end
    
end