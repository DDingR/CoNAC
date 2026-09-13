z = zeros(length(u), 1); % auxiliary state initialization

%% NEURAL NETOWORK INITIAILIZATION
nn.th = (rand(opt.th_size,1)-1/2)*2*opt.init_range;     % NN weight initialization
nn.tape = zeros(opt.tp_size, 1);                        % gradient tape for storing intermediate values during forward propagation, used for backpropagation