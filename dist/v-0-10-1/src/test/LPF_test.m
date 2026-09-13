clear

%%
T = 10;
dt  = 1e-4;
t = 0:dt:T;

t_len = length(t);

%%
LPF_tau = 1/(2*pi * 1);
LPF_Ts = dt;
% LPF_Ts = 1;


%%
% sig = @(x) sin(x);
sig = @(x) heaviside(x-2);
u = sig(0);
u_LPF = 0;

u_list = zeros(1, t_len);
uLPF_list = zeros(1, t_len);
u_list(1) = u;
uLPF_list(1) = u_LPF;

%%
for idx = 2:1:t_len

    u = sig(t(idx));

    u_LPF = (LPF_tau * u_LPF + LPF_Ts * u) / (LPF_tau + LPF_Ts);


    u_list(idx) = u;
    uLPF_list(idx) = u_LPF;
end

%%
figure(1); clf
plot(t, u_list, 'r', 'LineWidth', 2); hold on
plot(t, uLPF_list, 'b', 'LineWidth', 2); hold off