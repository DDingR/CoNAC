function Vn = nnWeightNorm(th, nnOpt)

    Vn = zeros(nnOpt.l_size-1,1);

    cumsum_V = [0;cumsum(nnOpt.th_size_list)];
    for l_idx = 1:1:nnOpt.l_size-1
        start_pt = cumsum_V(l_idx)+1;
        end_pt = cumsum_V(l_idx+1);   

        Vn(l_idx) = norm(th(start_pt: end_pt));
    end
    
end