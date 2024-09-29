
%% PREPARE
font_size = 16;
line_width = 2;
lgd_size = 12;

%% DATA EXTRACTION
cstr = nnOpt.cstr;

x1 = x_hist(1,:);
x2 = x_hist(2,:);
xd1 = x_hist(3,:);
xd2 = x_hist(4,:);

r1 = r_hist(1,:);
r2 = r_hist(2,:);
rd1 = rd_hist(1,:);
rd2 = rd_hist(2,:);

u1 = u_hist(1,:);
u2 = u_hist(2,:);
u1_sat = uSat_hist(1,:);
u2_sat = uSat_hist(2,:);

th = V_hist;
L = L_hist;

dot_L = dot_L_hist;

%% 
plotter1