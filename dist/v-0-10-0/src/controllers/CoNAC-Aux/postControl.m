
% auxiliary system update
del_u = zeros(2,1);
if u(1) > opt.cstr.uMax1
    del_u(1) = u(1) - opt.cstr.uMax1;
elseif u(1) < -opt.cstr.uMax1
    del_u(1) = u(1) + opt.cstr.uMax1;
end
if u(2) > opt.cstr.uMax2
    del_u(2) = u(2) - opt.cstr.uMax2;
elseif u(2) < -opt.cstr.uMax2
    del_u(2) = u(2) + opt.cstr.uMax2;
end

z = z + (opt.A_zeta*z + opt.B_zeta*del_u) * opt.dt;

% NN backward propagation, update NN weights
[nn, opt] = nnBackward(nn, opt, r+z, u_NN);
