function u = sat(u, opt, CONTROL_NUM)
    % u = u; return;
    u_ball = opt.cstr.u_ball;
    uMax1 = opt.cstr.uMax1;
    uMax2 = opt.cstr.uMax2;

    if (u(2) > uMax2 || u(2) < -uMax2)
        uMax2 = uMax2*sign(u(2));
        u(2) = uMax2;
    end

    if norm(u) > u_ball
        u = u/norm(u)*u_ball;
    end

end


    % if (u(1) > uMax1 || u(1) < -uMax1)
    %     uMax1 = uMax1*sign(u(1));
    %     u(1) = uMax1;
    % end


    % if u(1) > uMax1 || u(1) < -uMax1
    %     uMax1 = uMax1*sign(u(1));
    %     u(2) = uMax1*u(2)/u(1);
    %     u(1) = uMax1;
    % end