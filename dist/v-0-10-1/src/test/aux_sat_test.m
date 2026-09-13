% ********************************************************
%    * Optimization based Neuro Adaptive Controller *
%
%   Crafted By - Ryu Myeongseok
%   Version 0.1
%   Date : 2024.08.-  
%
% ********************************************************
 
%% INITIALIZATION
% ********************************************************
clear  
clc

addpath("../dynamics_model/")
%% SIMULATION SETTING
% ********************************************************
dt = 1e-4;
% T = 0.01;
T = 10;
% T = 15;
% T = 200;
t = 0:dt:T;
rpt_dt = 1;

[grad_x, IC] = model1_load();
[r_f, rd_f, rdd_f] = ref1_load(); % sin func

%% INITIAL CONDITION
% ********************************************************
x = IC.x;
u = IC.u;

%% RECORDER
% ********************************************************
num_x = length(x); num_u = length(u);
num_t = length(t);

x_hist = zeros(num_x, num_t); x_hist(:, 1) = x;
r_hist = zeros(num_x/2, num_t); r_hist(:, 1) = r_f(x,0);
u_hist = zeros(num_u, num_t); u_hist(:, 1) = u;
uSat_hist = zeros(num_u, num_t); uSat_hist(:, 1) = u;
aux_hist = zeros(2, num_t); 

k1 = 10; 
k2 = 10;

delta = 1e2;
zeta = rand(2,1)/norm(rand(2,1))* sqrt(delta);
pre_zeta = zeta;
aux_hist(:, 1) = zeta;

for t_idx = 2:1:num_t
    
    r1 = r_f(x, t(t_idx));
    r1d = rd_f(x, t(t_idx));
    r1dd = rdd_f(x, t(t_idx));
  
    % e = x-[r;rd];
    e1 = r1 - x(1:2);
    r2 = k1*e1 + r1d;
    e2 = r2 - x(3:4);
    e = [e1;e2];
    r2d = r1dd + k1*(x(3:4) + r1d);
       
    M = eye(2);
    C = eye(2);
    G = zeros(2,1);
    % k2 = 1e1;
    
    K_zeta = diag([1 1]) * 1e1 * 3;
    Kp = diag([1 1]) * 1e0  * 5;
%     k2 = 1e1 * 10;
%     Kp = k2;

    u1 = + M * k2 * e2 + M * e1 - C * x(3:4) + G*0 - M * r2d;
%     u1 = 50 * sin(t(t_idx)) * [1;0.1];
    u2 = + M * Kp * zeta;
    u = u1 + u2;
%     u = u1;

    u_bar = 52;
    u_sat = min(max(u, -u_bar),u_bar);

    u_del = u - u_sat;
    u_del = -u_del;
    u_del = u_del * 20;

    if norm(zeta) >= sqrt(delta)
%     if 1


        grad_zeta = (-K_zeta * zeta ...
                - ( abs(e2'*u_del)+0.5*(u_del'*u_del) )/norm(zeta)^2*zeta + u_del);

        assert(~isnan(norm(grad_zeta)));

        pre_zeta = zeta;
        zeta = zeta + grad_zeta * dt;

        
%         if norm(zeta) < sqrt(delta)
%             zeta = zeta/norm(zeta) * sqrt(delta);
%         end

    else
        
        grad_zeta = (-K_zeta * zeta + u_del);

        assert(~isnan(norm(grad_zeta)));

        pre_zeta = zeta;
        zeta = zeta + grad_zeta * dt;

        
%         if norm(zeta) < sqrt(delta)
%             zeta = zeta/norm(zeta) * sqrt(delta);
%         end
    end
    
    % p = 100;
    % u_sat = u / norm(u) * u_bar * ( norm(u)/u_bar / (1+(norm(u)/u_bar)^p)^(1/p) );

    % gradient calculation
    grad = grad_x(x, u_sat, t(t_idx));

    % error check
    assert(~isnan(norm(u)));
    assert(~isnan(norm(grad)));

    % step forward
    x = x + grad * dt;
    
    % record
    x_hist(:, t_idx) = x;
    r_hist(:, t_idx) = r1;
    u_hist(:, t_idx) = u;
    uSat_hist(:, t_idx) = u_sat;
    aux_hist(:, t_idx) = zeta;
    
end


figure(1);clf;
tiledlayout(2,1);
foc_ti = 3.57;
foc_tf = 3.58;
foc_ti = 1;
foc_tf = T;
foc_t = foc_ti/dt:foc_tf/dt;
foc_t = uint32(foc_t);

nexttile;
plot(t(foc_t), x_hist(1,foc_t), 'r', 'LineWidth', 2); hold on;
plot(t(foc_t), r_hist(1,foc_t), 'b', 'LineWidth', 2); hold off;
legend('x1', 'r1');
grid on

nexttile
plot(t(foc_t), x_hist(2,foc_t), 'r', 'LineWidth', 2); hold on;
plot(t(foc_t), r_hist(2,foc_t), 'b', 'LineWidth', 2); hold off;
legend('x2', 'r2');
grid on

figure(2);clf;
tiledlayout(2,1);
nexttile
plot(t(foc_t), u_hist(1,foc_t), 'r', 'LineWidth', 2); hold on;
plot(t(foc_t), uSat_hist(1,foc_t), 'b', 'LineWidth', 2); hold off;
legend('u', 'uSat');
% ylim([-70, 70]);
grid on

nexttile
plot(t(foc_t), u_hist(2,foc_t), 'r', 'LineWidth', 2); hold on;
plot(t(foc_t), uSat_hist(2,foc_t), 'b', 'LineWidth', 2); hold off; 
legend('u', 'uSat');
% ylim([-70, 70]);
grid on

figure(3);clf;
tiledlayout(2,1)
nexttile
plot(t(foc_t), aux_hist(1,foc_t), 'r', 'LineWidth', 2); hold on;
grid on
nexttile
plot(t(foc_t), aux_hist(2,foc_t), 'r', 'LineWidth', 2); hold on;
grid on

figure(4);clf
plot(aux_hist(1,foc_t), aux_hist(2,foc_t))
pbaspect([1 1 1])
grid on
