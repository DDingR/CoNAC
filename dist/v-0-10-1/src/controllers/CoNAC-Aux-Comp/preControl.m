
x_in = [x1;x2;r];
% x_in = [x1;x2];

[nn, u_NN] = nnForward(nn, opt, x_in);
u = u_NN;

% G = -eye(2) *5e-1; % control input matrix
% u = u - G*(z);
% % u = u - G*(z - r);



