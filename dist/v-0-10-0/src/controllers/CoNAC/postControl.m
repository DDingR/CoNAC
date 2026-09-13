% [nn, opt] = nnBackward(nn, opt, r, u_NN);       % NN backward propagation, update NN weights and Lagrange multipliers
[nn, opt] = nnBackward(nn, opt, fil_x-fil_xd, u_NN);       % NN backward propagation, update NN weights and Lagrange multipliers
