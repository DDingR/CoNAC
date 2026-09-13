% ---------------------------------------------
% Linear auxiliary system update
% ---------------------------------------------
K1 = diag([9 6]) * 1e0;
K2 = diag([10 5]) * 1e0;
K3 = diag([1 1]) * -diag([000, 3]);
% K3 = K3 * 0;

M0 = diag([1 1]) * 1;

% virtual control
K = diag([10 5]) * 1e-2;
u = u_NN - K * z;

% auxiliary system update
del_u(1,1) = min(max(u(1), -opt.cstr.uMax1), opt.cstr.uMax1) - u(1);
del_u(2,1) = min(max(u(2), -opt.cstr.uMax2), opt.cstr.uMax2) - u(2);
del_u = -del_u;

% NONLINEAR AUXILIARY SYSTEM UPDATE
A_z = 1e1*eye(2); B_z = 1e4*eye(2);
z_grad = ...
    - A_z * z ...
    + B_z * del_u;    

% z = z + z_grad * ctrl_dt;

% exact discretization
Ad_z = expm(-A_z*ctrl_dt);
z = Ad_z * z + -A_z \ (Ad_z - eye(2)) * B_z*del_u;

% NN backward propagation, update NN weights
[nn, opt] = nnBackward(nn, opt, r+1*z, u_NN);


% % ---------------------------------------------
% % Linear auxiliary system update
% % ---------------------------------------------
% K1 = diag([9 6]) * 1e0;
% K2 = diag([10 5]) * 1e0;
% K3 = diag([1 1]) * -diag([000, 3]);
% % K3 = K3 * 0;
% 
% M0 = diag([1 1]) * 1;
% 
% % virtual control
% alp = xd2 - K1 * e1;
% fil_alp_grad = (alp - opt.fil_alp) / 0.5;
% opt.fil_alp = opt.fil_alp + fil_alp_grad * ctrl_dt;
% opt.fil_alp = alp; fil_alp = [0;0];
% BSC_e2 = x2-opt.fil_alp;
% % u = -e1-K2*BSC_e2-K3*z -M0*( - fil_alp_grad);
% u = -e1-K2*BSC_e2-K3*z -M0*( -u_NN- fil_alp_grad);
% 
% % auxiliary system update
% del_u(1,1) = min(max(u(1), -opt.cstr.uMax1), opt.cstr.uMax1) - u(1);
% del_u(2,1) = min(max(u(2), -opt.cstr.uMax2), opt.cstr.uMax2) - u(2);
% 
% % NONLINEAR AUXILIARY SYSTEM UPDATE
% A_z = 1e1*eye(2); B_z = 1.5e2*eye(2);
% 
% z_grad = ...
%     - A_z * z ...
%     + B_z * del_u;    
% 
% % z = z + z_grad * ctrl_dt;
% 
% % exact discretization
% Ad_z = expm(-A_z*ctrl_dt);
% z = Ad_z * z + -A_z \ (Ad_z - eye(2)) * B_z*del_u;
% 
% % NN backward propagation, update NN weights
% [nn, opt] = nnBackward(nn, opt, BSC_e2-1*z, u_NN);

% ---------------------------------------------
% FIXED TIME
% ---------------------------------------------
% % z = zeros(size(z));
% 
% z1 = e1;
% 
% % p1 = 11/9;
% p1 = 1.01;
% p2 = 1/p1;
% 
% beta_1i     = p1 .^ sign( abs(z1) - 1 );
% beta_2i     = p2 .^ sign( 1 - abs(z1) );
% alp_gamma_1 = p1 .^ sign( abs(z) - 1 );
% alp_gamma_2 = p2 .^ sign( 1 - abs(z) );
% 
% l1 = 3; 
% l3 = 8;
% 
% sig_z = 1.1;
% k_z = 0.1;
% 
% l2 = l1;
% l4 = l3;
% 
% l_z = k_z;
% 
% sig_alp = @(x, alp) abs(x).^alp .* sign(x);
% 
% alp = xd2 - l1*sig_alp(z1, beta_1i) - l2*sig_alp(z1, beta_2i) - z1;
% 
% 
% 
% z2 = x2 - alp - z;
% 
% beta_3i     = p1 .^ sign( abs(z2) - 1 );
% beta_4i     = p2 .^ sign( 1 - abs(z2) );
% 
% % M0 = diag([1 1]) * 2.465;
% M0 = eye(2)*1;
% 
% u_NN = zeros(size(u));
% u_fix = M0 * ( -z1 - l3*sig_alp(z2, beta_3i) - l4*sig_alp(z2, beta_4i) - z2);
% u_aux = M0 * ( - sig_z * z - k_z * sig_alp(z, alp_gamma_1) - l_z * sig_alp(z, alp_gamma_2) );
% u = u_fix + u_aux - M0 * u_NN;
% 
% % % del_u(1,1) = min(max(u(1), -opt.cstr.uMax1), opt.cstr.uMax1) - u(1);
% del_u(1,1) = 0;
% del_u(2,1) = min(max(u(2), -opt.cstr.uMax2), opt.cstr.uMax2) - u(2);
% % del_u = -del_u;
% % del_u = del_u * 10;
% 
% z_grad = -k_z * sig_alp(z, alp_gamma_1) - l_z * sig_alp(z, alp_gamma_2) - sig_z * z + M0\del_u;
% z = z + z_grad * ctrl_dt;
% 
% % NN backward propagation, update NN weights
% [nn, opt] = nnBackward(nn, opt, -z2, u_NN);
% 
% % ---------------------------------------------
% % NONLINEAR TIME
% % ---------------------------------------------
% K1 = diag([5 1]) * 5e0;
% K2 = diag([1 1]) * 15e0; 
% K3 = diag([1 1]) * -1.5e3;
% B1 = diag([1 1]) * 1e2;
% B2 = diag([1 1]) * 5e-3;
% 
% M0 = diag([1 1]) * 2;
% 
% % virtual control
% alp = xd2 - K1 * e1;
% fil_alp_grad = (alp - opt.fil_alp) / 0.1;
% opt.fil_alp = opt.fil_alp + fil_alp_grad * ctrl_dt;
% % opt.fil_alp = alp; fil_alp = [0;0];
% BSC_e2 = x2-opt.fil_alp;
% % u = -e1-K2*BSC_e2-K3*z -M0*( - fil_alp_grad);
% u = -e1-K2*BSC_e2-K3*z -M0*( -u_NN- fil_alp_grad);
% 
% % auxiliary system update
% del_u(1,1) = min(max(u(1), -opt.cstr.uMax1), opt.cstr.uMax1) - u(1);
% del_u(2,1) = min(max(u(2), -opt.cstr.uMax2), opt.cstr.uMax2) - u(2);
% % del_u = -del_u;
% 
% % NONLINEAR AUXILIARY SYSTEM UPDATE
% mu_z = 1e2;
% mu_du = 1e-8;
% 
% % if norm(z) <= mu_z && norm(del_u) <= mu_du
% %     cpx_ = 0;
% % else
% %     den_z = max(norm(z)^2, mu_z^2);
% %     cpx_ = ((abs(BSC_e2'*del_u) + 1/2*(del_u'*B2*del_u)) / (den_z));   
% 
% %     Ad_z = expm((-B1-cpx_*eye(2))*ctrl_dt);
% %     z = Ad_z * z + (-B1-cpx_*eye(2)) \ (Ad_z - eye(2)) * B2 * del_u;
% % end
% cpx_ = [0;0];
% for z_idx = 1:2
%     if abs(z(z_idx)) <= mu_z && abs(del_u(z_idx)) <= mu_du
%         cpx_(z_idx) = 0;
%     else
%         den_z = max(abs(z(z_idx)), mu_z);
%         cpx_(z_idx) = ((abs(BSC_e2(z_idx)*del_u(z_idx)) + 1/2*(del_u(z_idx)*B2(z_idx,z_idx)*del_u(z_idx))) / (den_z));   
%         Ad_z = expm((-B1(z_idx,z_idx)-cpx_(z_idx))*ctrl_dt);
%         z(z_idx) = Ad_z * z(z_idx) + (-B1(z_idx,z_idx)-cpx_(z_idx)) \ (Ad_z - 1) * B2(z_idx,z_idx) * del_u(z_idx);
%     end
% end
% 
% 
% 
% % NN backward propagation, update NN weights
% [nn, opt] = nnBackward(nn, opt, BSC_e2, u_NN);
