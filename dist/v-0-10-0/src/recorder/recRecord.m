x1_hist(:, t_idx) = x1;
x2_hist(:, t_idx) = x2;
xd1_hist(:, t_idx) = xd1;
xd2_hist(:, t_idx) = xd2;

u_hist(:, t_idx) = u;
uSat_hist(:, t_idx) = u_sat;
th_hist(:, t_idx) = nnWeightNorm(nn.th, opt);

if CONTROL_NUM == 1
    % CoNAC
    lbd_hist(:, t_idx) = opt.lbd;
elseif CONTROL_NUM == 2 || CONTROL_NUM == 3
    % Aux.
    zeta_hist(:, t_idx) = z;
end
