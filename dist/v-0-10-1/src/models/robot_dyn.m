function [M, C, G, F] = robot_dyn(q, q_dot, p, dof)
%ROBOT_DYN Planar serial manipulator dynamics for 1R/2R/3R models.
%   The rigid-body terms follow the symbolic models exported in
%   03_code/Exp/mathematica/planar_{1R,2R,3R}_dynamics.m.
%   Viscous + Stribeck friction terms and optional motor inertias are added
%   on top of the rigid-body model using fields from the parameter struct p.

q = q(:);
q_dot = q_dot(:);

if nargin < 4 || isempty(dof)
    dof = numel(q);
end

if numel(q) < dof || numel(q_dot) < dof
    error('robot_dyn:DimensionMismatch', ...
        'q and q_dot must contain at least %d elements.', dof);
end

g = 9.81;

switch dof
    case 1
        q1 = q(1);

        M = p.I1 + p.lc1^2 * p.m1 + get_field_local(p, 'I1m', 0.0);
        C = 0;
        G = g * p.lc1 * p.m1 * cos(q1);

    case 2
        q1 = q(1);
        q2 = q(2);
        q_dot1 = q_dot(1);
        q_dot2 = q_dot(2);

        M11 = p.I1 + p.I2 + p.lc1^2 * p.m1 + p.l1^2 * p.m2 + p.lc2^2 * p.m2 + ...
            2 * p.l1 * p.lc2 * p.m2 * cos(q2) + get_field_local(p, 'I1m', 0.0);
        M12 = p.I2 + p.lc2 * (p.l1 * cos(q2) + p.lc2) * p.m2;
        M21 = M12;
        M22 = p.I2 + p.lc2^2 * p.m2 + get_field_local(p, 'I2m', 0.0);

        C11 = -p.l1 * p.lc2 * p.m2 * q_dot2 * sin(q2);
        C12 = -p.l1 * p.lc2 * p.m2 * (q_dot1 + q_dot2) * sin(q2);
        C21 = p.l1 * p.lc2 * p.m2 * q_dot1 * sin(q2);
        C22 = 0;

        G1 = g * (p.lc1 * p.m1 * cos(q1) + p.l1 * p.m2 * cos(q1) + p.lc2 * p.m2 * cos(q1 + q2));
        G2 = g * p.lc2 * p.m2 * cos(q1 + q2);

        M = [M11 M12; M21 M22];
        C = [C11 C12; C21 C22];
        G = [G1; G2];

    case 3
        q1 = q(1);
        q2 = q(2);
        q3 = q(3);
        q_dot1 = q_dot(1);
        q_dot2 = q_dot(2);
        q_dot3 = q_dot(3);

        c2 = cos(q2);
        c3 = cos(q3);
        c23 = cos(q2 + q3);
        s2 = sin(q2);
        s3 = sin(q3);
        s23 = sin(q2 + q3);

        h2 = p.l1 * ((p.lc2 * p.m2 + p.l2 * p.m3) * s2 + p.lc3 * p.m3 * s23);
        h3 = p.lc3 * p.m3 * (p.l2 * s3 + p.l1 * s23);
        h23 = p.l2 * p.lc3 * p.m3 * s3;

        M11 = p.I1 + p.I2 + p.I3 + p.lc1^2 * p.m1 + p.l1^2 * p.m2 + p.lc2^2 * p.m2 + ...
            p.l1^2 * p.m3 + p.l2^2 * p.m3 + p.lc3^2 * p.m3 + ...
            2 * p.l1 * (p.lc2 * p.m2 + p.l2 * p.m3) * c2 + ...
            2 * p.l2 * p.lc3 * p.m3 * c3 + ...
            2 * p.l1 * p.lc3 * p.m3 * c23 + ...
            get_field_local(p, 'I1m', 0.0);
        M12 = p.I2 + p.I3 + p.lc2^2 * p.m2 + p.l2^2 * p.m3 + p.lc3^2 * p.m3 + ...
            p.l1 * (p.lc2 * p.m2 + p.l2 * p.m3) * c2 + ...
            2 * p.l2 * p.lc3 * p.m3 * c3 + ...
            p.l1 * p.lc3 * p.m3 * c23;
        M13 = p.I3 + p.lc3^2 * p.m3 + p.l2 * p.lc3 * p.m3 * c3 + p.l1 * p.lc3 * p.m3 * c23;
        M21 = M12;
        M22 = p.I2 + p.I3 + p.lc2^2 * p.m2 + p.l2^2 * p.m3 + p.lc3^2 * p.m3 + ...
            2 * p.l2 * p.lc3 * p.m3 * c3 + get_field_local(p, 'I2m', 0.0);
        M23 = p.I3 + p.lc3^2 * p.m3 + p.l2 * p.lc3 * p.m3 * c3;
        M31 = M13;
        M32 = M23;
        M33 = p.I3 + p.lc3^2 * p.m3 + get_field_local(p, 'I3m', 0.0);

        C11 = -h2 * q_dot2 - h3 * q_dot3;
        C12 = -h2 * (q_dot1 + q_dot2) - h3 * q_dot3;
        C13 = -h3 * (q_dot1 + q_dot2 + q_dot3);
        C21 = h2 * q_dot1 - h23 * q_dot3;
        C22 = -h23 * q_dot3;
        C23 = -h23 * (q_dot1 + q_dot2 + q_dot3);
        C31 = p.lc3 * p.m3 * ((p.l2 * s3 + p.l1 * s23) * q_dot1 + p.l2 * s3 * q_dot2);
        C32 = h23 * (q_dot1 + q_dot2);
        C33 = 0;

        G1 = g * ((p.lc1 * p.m1 + p.l1 * (p.m2 + p.m3)) * cos(q1) + ...
            (p.lc2 * p.m2 + p.l2 * p.m3) * cos(q1 + q2) + ...
            p.lc3 * p.m3 * cos(q1 + q2 + q3));
        G2 = g * ((p.lc2 * p.m2 + p.l2 * p.m3) * cos(q1 + q2) + ...
            p.lc3 * p.m3 * cos(q1 + q2 + q3));
        G3 = g * p.lc3 * p.m3 * cos(q1 + q2 + q3);

        M = [M11 M12 M13; M21 M22 M23; M31 M32 M33];
        C = [C11 C12 C13; C21 C22 C23; C31 C32 C33];
        G = [G1; G2; G3];

    otherwise
        error('Unsupported number of degrees of freedom: %d', dof);
end

F = friction_model(q_dot, p, dof);
end

function F = friction_model(q_dot, p, dof)
switch dof
    case 1
        b = p.b1;
        fc = p.fc1;
        Fs = p.Fs1;
        vs = p.vs1;
        k = p.k1;

    case 2
        b = [p.b1; p.b2];
        fc = [p.fc1; p.fc2];
        Fs = [p.Fs1; p.Fs2];
        vs = [p.vs1; p.vs2];
        k = [p.k1; p.k2];

    case 3
        b = [p.b1; p.b2; p.b3];
        fc = [p.fc1; p.fc2; p.fc3];
        Fs = [p.Fs1; p.Fs2; p.Fs3];
        vs = [p.vs1; p.vs2; p.vs3];
        k = [p.k1; p.k2; p.k3];
end

q_dot = q_dot(1:dof);
stribeck = fc + (Fs - fc) .* exp(-(q_dot ./ vs).^2);
F = stribeck .* tanh(k .* q_dot) + b .* q_dot;
end

function value = get_field_local(s, name, default_value)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default_value;
end
end
