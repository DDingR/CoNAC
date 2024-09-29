clear 

t = 0:0.01:60;

e = 1;
M = 50;
m = -50;

alp = (M+m)/2;
mu = (M-m)/2;
k = 100;

tanh_t = alp + mu * ( ((t-alp)/mu) ./ ((1+((t-alp)/mu).^k).^(1/k)) );
% tanh_t = (1-a) * M* tanh( (t-a*M) ./ ((1-a)*M) );
% tanh_t =  M* tanh( t./M );


clf
plot(t, tanh_t); hold on
plot(t, ones(size(t))*M); hold on
plot(t, t); hold on

%%
clear

% syms x [2 1] matrix
syms x1 x2

M = 50; p = 100;

dfdx_11(x1,x2) = diff( ...
    x1 / (1+(sqrt(x1^2+x2^2)/M)^p).^(1/p) , x1);
    % dfdx = diff(x1*x2, x1)
dfdx_12(x1,x2) = diff( ...
    x1 / (1+(sqrt(x1^2+x2^2)/M)^p).^(1/p) , x2);
    % dfdx = diff(x1*x2, x1)
dfdx_21(x1,x2) = diff( ...
    x2 / (1+(sqrt(x1^2+x2^2)/M)^p).^(1/p) , x1);
    % dfdx = diff(x1*x2, x1)
dfdx_22(x1,x2) = diff( ...
    x2 / (1+(sqrt(x1^2+x2^2)/M)^p).^(1/p) , x2);
    % dfdx = diff(x1*x2, x1)

matlabFunction(dfdx_11, "File", "dfdx_11", "Vars", {x1,x2});
matlabFunction(dfdx_12, "File", "dfdx_12", "Vars", {x1,x2});
matlabFunction(dfdx_21, "File", "dfdx_21", "Vars", {x1,x2});
matlabFunction(dfdx_22, "File", "dfdx_22", "Vars", {x1,x2});
%%

% x = [1;1];
x1 = 49;
x2 = 1;


vpa(subs(dfdx(x1,x2)))
subs(x2);
pretty(dfdx)
vpa(dfdx)
