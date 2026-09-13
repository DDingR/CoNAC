
%% PREPARE
font_size = 8;
line_width = 2;
lgd_size = 12;
% fig_height = 210; fig_width = 450;
fig_height = 210*4; fig_width = 450*2;

%% DATA EXTRACTION
cstr = opt.cstr;

x1 = x1_hist(1,:);
x2 = x1_hist(2,:);
xd1 = x2_hist(1,:);
xd2 = x2_hist(2,:);

r1 = xd1_hist(1,:);
r2 = xd1_hist(2,:);
rd1 = xd2_hist(1,:);
rd2 = xd2_hist(2,:);

u1 = u_hist(1,:);
u2 = u_hist(2,:);
u1_sat = uSat_hist(1,:);
u2_sat = uSat_hist(2,:);

th = th_hist;

if CONTROL_NUM == 1
    % CoNAC
    L = lbd_hist;
elseif CONTROL_NUM == 2
    % Aux.
    z1 = zeta_hist(1,:);
    z2 = zeta_hist(2,:);
end

%% 
plotter_long_1
% if CONTROL_NUM == 1
%     plotter_long_1
% elseif CONTROL_NUM == 2
%     plotter2
% end