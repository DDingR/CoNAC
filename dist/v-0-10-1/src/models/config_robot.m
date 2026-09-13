function p = config_robot(dof)
%CONFIG_ROBOT Robot model parameters aligned with robot3 firmware config.h.

if nargin < 1 || isempty(dof)
    dof = 2;
end

if ~isscalar(dof) || ~ismember(dof, [1 2 3])
    error('config_robot:InvalidDof', 'dof must be one of 1, 2, or 3.');
end

% 03_code/Exp/firmware/robot3/src/common/config.h
link1_1d = struct( ...
    'l', 0.2, ...
    'lc', 8.658e-2, ...
    'm', 1.330, ...
    'Ixx', 1.248e5 * 1.0e-7, ...
    'Im', 9.0 * 9.0 * 1002.0e-7, ...
    'b', 0.0766972627, ...
    'fc', 0.7101160378, ...
    'Fs', 0.9651779394, ...
    'vs', 0.2245092542, ...
    'k', 50.0);

% Link 2 includes a 0.280 kg point mass at 0.20 m from joint 2.
links_2d = [ ...
    struct( ...
        'l', 0.2, 'lc', 0.13888, 'm', 2.465, ...
        'Ixx', 0.06911, 'Im', 9.0 * 9.0 * 1002.0e-7, ...
        'b', 0.0766972627, 'fc', 0.7101160378, ...
        'Fs', 0.9651779394, 'vs', 0.2245092542, 'k', 50.0), ...
    struct( ...
        'l', 0.2, 'lc', 0.13888, 'm', 2.465, ...
        'Ixx', 0.06911, 'Im', 9.0 * 9.0 * 1002.0e-7, ...
        'b', 0.0766972627, 'fc', 0.7101160378, ...
        'Fs', 0.9651779394, 'vs', 0.2245092542, 'k', 50.0)];

% Link 2 includes the same additional point mass as the 2-DOF assembly.
% Joint 2 and Joint 3 use the Joint 1 friction result provisionally.
links_3d = [ ...
    struct( ...
        'l', 0.2, 'lc', 0.13888, 'm', 2.465, ...
        'Ixx', 0.06911, 'Im', 9.0 * 9.0 * 1002.0e-7, ...
        'b', 0.0766972627, 'fc', 0.7101160378, ...
        'Fs', 0.9651779394, 'vs', 0.2245092542, 'k', 50.0), ...
    struct( ...
        'l', 0.2, 'lc', 0.13888, 'm', 2.465, ...
        'Ixx', 0.06911, 'Im', 9.0 * 9.0 * 1002.0e-7, ...
        'b', 0.0766972627, 'fc', 0.7101160378, ...
        'Fs', 0.9651779394, 'vs', 0.2245092542, 'k', 50.0), ...
    struct( ...
        'l', 0.2, 'lc', 6.501e-2, 'm', 0.870264, ...
        'Ixx', 48949.343 * 1.0e-7, 'Im', 9.0 * 9.0 * 1002.0e-7, ...
        'b', 0.0766972627, 'fc', 0.7101160378, ...
        'Fs', 0.9651779394, 'vs', 0.2245092542, 'k', 50.0)];

p = struct();
p.dof = dof;
p.g = 9.81;
p.dt = 0.002;
p.kt_nom = 1.3137;
p.kt_nom_vec = p.kt_nom * ones(1, max(dof, 1));

switch dof
    case 1
        p.l1 = link1_1d.l;
        p.lc1 = link1_1d.lc;
        p.m1 = link1_1d.m;
        p.I1 = link1_1d.Ixx;
        p.I1m = link1_1d.Im;
        p.b1 = link1_1d.b;
        p.fc1 = link1_1d.fc;
        p.Fs1 = link1_1d.Fs;
        p.vs1 = link1_1d.vs;
        p.k1 = link1_1d.k;

    case 2
        p.l1 = links_2d(1).l;
        p.l2 = links_2d(2).l;
        p.lc1 = links_2d(1).lc;
        p.lc2 = links_2d(2).lc;
        p.m1 = links_2d(1).m;
        p.m2 = links_2d(2).m;
        p.I1 = links_2d(1).Ixx;
        p.I2 = links_2d(2).Ixx;
        p.I1m = links_2d(1).Im;
        p.I2m = links_2d(2).Im;
        p.b1 = links_2d(1).b;
        p.fc1 = links_2d(1).fc;
        p.Fs1 = links_2d(1).Fs;
        p.vs1 = links_2d(1).vs;
        p.k1 = links_2d(1).k;
        p.b2 = links_2d(2).b;
        p.fc2 = links_2d(2).fc;
        p.Fs2 = links_2d(2).Fs;
        p.vs2 = links_2d(2).vs;
        p.k2 = links_2d(2).k;

    case 3
        p.l1 = links_3d(1).l;
        p.l2 = links_3d(2).l;
        p.l3 = links_3d(3).l;
        p.lc1 = links_3d(1).lc;
        p.lc2 = links_3d(2).lc;
        p.lc3 = links_3d(3).lc;
        p.m1 = links_3d(1).m;
        p.m2 = links_3d(2).m;
        p.m3 = links_3d(3).m;
        p.I1 = links_3d(1).Ixx;
        p.I2 = links_3d(2).Ixx;
        p.I3 = links_3d(3).Ixx;
        p.I1m = links_3d(1).Im;
        p.I2m = links_3d(2).Im;
        p.I3m = links_3d(3).Im;
        p.b1 = links_3d(1).b;
        p.fc1 = links_3d(1).fc;
        p.Fs1 = links_3d(1).Fs;
        p.vs1 = links_3d(1).vs;
        p.k1 = links_3d(1).k;
        p.b2 = links_3d(2).b;
        p.fc2 = links_3d(2).fc;
        p.Fs2 = links_3d(2).Fs;
        p.vs2 = links_3d(2).vs;
        p.k2 = links_3d(2).k;
        p.b3 = links_3d(3).b;
        p.fc3 = links_3d(3).fc;
        p.Fs3 = links_3d(3).Fs;
        p.vs3 = links_3d(3).vs;
        p.k3 = links_3d(3).k;
end

% MATLAB-side controller defaults.
p.Kp = diag(25 * ones(1, dof));
p.Kd = diag(10 * ones(1, dof));

