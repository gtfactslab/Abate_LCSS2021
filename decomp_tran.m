


% decomposition function for closed-loop dynamics (WITH linear transformation)

function out = decomp_tran(x, xh, w, wh, J, Jinv, umax, kp, kd)
    out = zeros(7,1);
    out(1:4, 1) = d_bits_4(x, xh, Jinv);
    out(5:7, 1) = d_bits_1(x, xh, Jinv) + sigma(d_bits_2(x, xh, Jinv) + d_bits_3(x, xh, J) - kd*x(5:7)) + w;



    function out = sigma(in)
        out = umax*tanh((1/umax)*in);
    end

    % -cross(Jinv*om, om)
    function out = d_bits_1(x, xh, Jinv)
    om1 = x(5);
    om2 = x(6);
    om3 = x(7);
    oh1 = xh(5);
    oh2 = xh(6);
    oh3 = xh(7);

    Jinv1_1 = Jinv(1, 1);
    Jinv2_1 = Jinv(2, 1);
    Jinv3_1 = Jinv(3, 1);
    Jinv2_2 = Jinv(2, 2);
    Jinv3_2 = Jinv(2, 3);
    Jinv3_3 = Jinv(3, 3);

    a = zeros(5,1);
    b = zeros(5,1);
    c = zeros(5,1);
    out = zeros(3,1);


    if prod(x <= xh)
        % first part: om1_dot
        % om2*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3) - om3*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3)
        % or
        % Jinv3_1*om1*om2 + Jinv3_2*om2^2 + (Jinv3_3- Jinv2_2)*om3*om2 - Jinv2_1*om1*om3 - Jinv3_2*om3^2



        a(1) = min([ Jinv3_1*om1*om2, Jinv3_1*om1*oh2]);

        if om2 <= 0 && 0 <= oh2
            a(2) = min([Jinv3_2*om2^2, Jinv3_2*oh2^2, 0]);
        else
            a(2) = min([Jinv3_2*om2^2, Jinv3_2*oh2^2]);
        end

        a(3) = min([ (Jinv3_3-Jinv2_2)*om3*om2, (Jinv3_3-Jinv2_2)*oh3*om2, (Jinv3_3-Jinv2_2)*om3*oh2, (Jinv3_3-Jinv2_2)*oh3*oh2]);

        a(4) = min([ -Jinv2_1*om1*om3, -Jinv2_1*om1*oh3]);

        if om3 <= 0 && 0 <= oh3
            a(5) = min([- Jinv3_2*om3^2, - Jinv3_2*oh3^2, 0]);
        else
            a(5) = min([- Jinv3_2*om3^2, - Jinv3_2*oh3^2]);
        end

        out(1, 1) = sum(a);


        % second part: om2_dot
        % om3*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3) - om1*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3)
        % or
        % (Jinv1_1 - Jinv3_3)*om1*om3 + Jinv2_1*om2*om3 + Jinv3_1*om3^2 - Jinv3_1*om1^2 - Jinv3_2*om1*om2
        b(1) = min([(Jinv1_1 - Jinv3_3)*om1*om3, (Jinv1_1 - Jinv3_3)*oh1*om3, (Jinv1_1 - Jinv3_3)*om1*oh3, (Jinv1_1 - Jinv3_3)*oh1*oh3]);
        b(2) = min([Jinv2_1*om2*om3, Jinv2_1*om2*oh3]);

        if om3 <= 0 && 0 <= oh3
            b(3) = min([Jinv3_1*om3^2, Jinv3_1*oh3^2, 0]);
        else
            b(3) = min([Jinv3_1*om3^2, Jinv3_1*oh3^2]);
        end
        if om1 <= 0 && 0 <= oh1
            b(4) = min([- Jinv3_1*om1^2, - Jinv3_1*oh1^2, 0]);
        else
            b(4) = min([- Jinv3_1*om1^2, - Jinv3_1*oh1^2]);
        end

        b(5) = min([- Jinv3_2*om1*om2, - Jinv3_2*oh1*om2]);

        out(2, 1) = sum(b);


        % third part: om3_dot
        % om1*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3) - om2*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3)
        % or
        % Jinv2_1*om1^2 + (Jinv2_2-Jinv1_1)*om1*om2 + Jinv3_2*om1*om3 - Jinv2_1*om2^2 - Jinv3_1*om2*om3

        if om1 <= 0 && 0 <= oh1
            c(1) = min([Jinv2_1*om1^2, Jinv2_1*oh1^2, 0]);
        else
            c(1) = min([Jinv2_1*om1^2, Jinv2_1*oh1^2]);
        end
        c(2) = min([(Jinv2_2-Jinv1_1)*om1*om2, (Jinv2_2-Jinv1_1)*oh1*om2, (Jinv2_2-Jinv1_1)*om1*oh2, (Jinv2_2-Jinv1_1)*oh1*oh2]);
        c(3) = min([Jinv3_2*om1*om3, Jinv3_2*oh1*om3]);
        if om2 <= 0 && 0 <= oh2
            c(4) = min([- Jinv2_1*om2^2, - Jinv2_1*oh2^2, 0]);
        else
            c(4) = min([- Jinv2_1*om2^2, - Jinv2_1*oh2^2]);
        end
        c(5) = min([- Jinv3_1*om2*om3, - Jinv3_1*oh2*om3]);
        out(3, 1) = sum(c);



    elseif prod(xh <= x)
        % first part: om1_dot
        % om2*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3) - om3*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3)
        % or
        % Jinv3_1*om1*om2 + Jinv3_2*om2^2 + (Jinv3_3- Jinv2_2)*om3*om2 - Jinv2_1*om1*om3 - Jinv3_2*om3^2

        a(1) = max([ Jinv3_1*om1*om2, Jinv3_1*om1*oh2]);

        if oh2 <= 0 && 0 <= om2
            a(2) = max([Jinv3_2*om2^2, Jinv3_2*oh2^2, 0]);
        else
            a(2) = max([Jinv3_2*om2^2, Jinv3_2*oh2^2]);
        end

        a(3) = max([ (Jinv3_3- Jinv2_2)*om3*om2, (Jinv3_3- Jinv2_2)*oh3*om2, (Jinv3_3- Jinv2_2)*om3*oh2, (Jinv3_3- Jinv2_2)*oh3*oh2]);

        a(4) = max([ -Jinv2_1*om1*om3, -Jinv2_1*om1*oh3]);

        if oh3 <= 0 && 0 <= om3
            a(5) = max([- Jinv3_2*om3^2, - Jinv3_2*oh3^2, 0]);
        else
            a(5) = max([- Jinv3_2*om3^2, - Jinv3_2*oh3^2]);
        end

        out(1, 1) = sum(a);


        % second part: om2_dot
        % om3*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3) - om1*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3)
        % or
        % (Jinv1_1 - Jinv3_3)*om1*om3 + Jinv2_1*om2*om3 + Jinv3_1*om3^2 - Jinv3_1*om1^2 - Jinv3_2*om1*om2
        b(1) = max([(Jinv1_1 - Jinv3_3)*om1*om3, (Jinv1_1 - Jinv3_3)*oh1*om3, (Jinv1_1 - Jinv3_3)*om1*oh3, (Jinv1_1 - Jinv3_3)*oh1*oh3]);
        b(2) = max([Jinv2_1*om2*om3, Jinv2_1*om2*oh3]);

        if oh3 <= 0 && 0 <= om3
            b(3) = max([Jinv3_1*om3^2, Jinv3_1*oh3^2, 0]);
        else
            b(3) = max([Jinv3_1*om3^2, Jinv3_1*oh3^2]);
        end
        if oh1 <= 0 && 0 <= om1
            b(4) = max([- Jinv3_1*om1^2, - Jinv3_1*oh1^2, 0]);
        else
            b(4) = max([- Jinv3_1*om1^2, - Jinv3_1*oh1^2]);
        end

        b(5) = max([- Jinv3_2*om1*om2, - Jinv3_2*oh1*om2]);

        out(2, 1) = sum(b);


        % third part: om3_dot
        % om1*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3) - om2*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3)
        % or
        % Jinv2_1*om1^2 + (Jinv2_2-Jinv1_1)*om1*om2 + Jinv3_2*om1*om3 - Jinv2_1*om2^2 - Jinv3_1*om2*om3

        if oh1 <= 0 && 0 <= om1
            c(1) = max([Jinv2_1*om1^2, Jinv2_1*oh1^2, 0]);
        else
            c(1) = max([Jinv2_1*om1^2, Jinv2_1*oh1^2]);
        end
        c(2) = max([(Jinv2_2-Jinv1_1)*om1*om2, (Jinv2_2-Jinv1_1)*oh1*om2, (Jinv2_2-Jinv1_1)*om1*oh2, (Jinv2_2-Jinv1_1)*oh1*oh2]);
        c(3) = max([Jinv3_2*om1*om3, Jinv3_2*oh1*om3]);
        if oh2 <= 0 && 0 <= om2
            c(4) = max([- Jinv2_1*om2^2, - Jinv2_1*oh2^2, 0]);
        else
            c(4) = max([- Jinv2_1*om2^2, - Jinv2_1*oh2^2]);
        end
        c(5) = max([- Jinv3_1*om2*om3, - Jinv3_1*oh2*om3]);

        out(3, 1) = sum(c);


    end

    end

    % cross(Jinv*om, om)
    function out = d_bits_2(x, xh, Jinv)
    om1 = x(5);
    om2 = x(6);
    om3 = x(7);
    oh1 = xh(5);
    oh2 = xh(6);
    oh3 = xh(7);

    Jinv1_1 = Jinv(1, 1);
    Jinv2_1 = Jinv(2, 1);
    Jinv3_1 = Jinv(3, 1);
    Jinv2_2 = Jinv(2, 2);
    Jinv3_2 = Jinv(2, 3);
    Jinv3_3 = Jinv(3, 3);

    a = zeros(5,1);
    b = zeros(5,1);
    c = zeros(5,1);
    out = zeros(3,1);


    if prod(x <= xh)
        % first part: om1_dot
        % om2*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3) - om3*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3)
        % or
        % Jinv3_1*om1*om2 + Jinv3_2*om2^2 + (Jinv3_3- Jinv2_2)*om3*om2 - Jinv2_1*om1*om3 - Jinv3_2*om3^2

        a(1) = -max([ Jinv3_1*om1*om2, Jinv3_1*om1*oh2]);

        if om2 <= 0 && 0 <= oh2
            a(2) = -max([Jinv3_2*om2^2, Jinv3_2*oh2^2, 0]);
        else
            a(2) = -max([Jinv3_2*om2^2, Jinv3_2*oh2^2]);
        end

        a(3) = -max([ (Jinv3_3-Jinv2_2)*om3*om2, (Jinv3_3-Jinv2_2)*oh3*om2, (Jinv3_3-Jinv2_2)*om3*oh2, (Jinv3_3-Jinv2_2)*oh3*oh2]);

        a(4) = -max([ -Jinv2_1*om1*om3, -Jinv2_1*om1*oh3]);

        if om3 <= 0 && 0 <= oh3
            a(5) = -max([- Jinv3_2*om3^2, - Jinv3_2*oh3^2, 0]);
        else
            a(5) = -max([- Jinv3_2*om3^2, - Jinv3_2*oh3^2]);
        end

        out(1, 1) = sum(a);


        % second part: om2_dot
        % om3*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3) - om1*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3)
        % or
        % (Jinv1_1 - Jinv3_3)*om1*om3 + Jinv2_1*om2*om3 + Jinv3_1*om3^2 - Jinv3_1*om1^2 - Jinv3_2*om1*om2
        b(1) = -max([(Jinv1_1 - Jinv3_3)*om1*om3, (Jinv1_1 - Jinv3_3)*oh1*om3, (Jinv1_1 - Jinv3_3)*om1*oh3, (Jinv1_1 - Jinv3_3)*oh1*oh3]);
        b(2) = -max([Jinv2_1*om2*om3, Jinv2_1*om2*oh3]);

        if om3 <= 0 && 0 <= oh3
            b(3) = -max([Jinv3_1*om3^2, Jinv3_1*oh3^2, 0]);
        else
            b(3) = -max([Jinv3_1*om3^2, Jinv3_1*oh3^2]);
        end
        if om1 <= 0 && 0 <= oh1
            b(4) = -max([- Jinv3_1*om1^2, - Jinv3_1*oh1^2, 0]);
        else
            b(4) = -max([- Jinv3_1*om1^2, - Jinv3_1*oh1^2]);
        end

        b(5) = -max([- Jinv3_2*om1*om2, - Jinv3_2*oh1*om2]);

        out(2, 1) = sum(b);


        % third part: om3_dot
        % om1*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3) - om2*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3)
        % or
        % Jinv2_1*om1^2 + (Jinv2_2-Jinv1_1)*om1*om2 + Jinv3_2*om1*om3 - Jinv2_1*om2^2 - Jinv3_1*om2*om3

        if om1 <= 0 && 0 <= oh1
            c(1) = -max([Jinv2_1*om1^2, Jinv2_1*oh1^2, 0]);
        else
            c(1) = -max([Jinv2_1*om1^2, Jinv2_1*oh1^2]);
        end
        c(2) = -max([(Jinv2_2-Jinv1_1)*om1*om2, (Jinv2_2-Jinv1_1)*oh1*om2, (Jinv2_2-Jinv1_1)*om1*oh2, (Jinv2_2-Jinv1_1)*oh1*oh2]);
        c(3) = -max([Jinv3_2*om1*om3, Jinv3_2*oh1*om3]);
        if om2 <= 0 && 0 <= oh2
            c(4) = -max([- Jinv2_1*om2^2, - Jinv2_1*oh2^2, 0]);
        else
            c(4) = -max([- Jinv2_1*om2^2, - Jinv2_1*oh2^2]);
        end
        c(5) = -max([- Jinv3_1*om2*om3, - Jinv3_1*oh2*om3]);
        out(3, 1) = sum(c);



    elseif prod(xh <= x)
        % first part: om1_dot
        % om2*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3) - om3*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3)
        % or
        % Jinv3_1*om1*om2 + Jinv3_2*om2^2 + (Jinv3_3- Jinv2_2)*om3*om2 - Jinv2_1*om1*om3 - Jinv3_2*om3^2

        a(1) = -min([ Jinv3_1*om1*om2, Jinv3_1*om1*oh2]);

        if oh2 <= 0 && 0 <= om2
            a(2) = -min([Jinv3_2*om2^2, Jinv3_2*oh2^2, 0]);
        else
            a(2) = -min([Jinv3_2*om2^2, Jinv3_2*oh2^2]);
        end

        a(3) = -min([ (Jinv3_3- Jinv2_2)*om3*om2, (Jinv3_3- Jinv2_2)*oh3*om2, (Jinv3_3- Jinv2_2)*om3*oh2, (Jinv3_3- Jinv2_2)*oh3*oh2]);

        a(4) = -min([ -Jinv2_1*om1*om3, -Jinv2_1*om1*oh3]);

        if oh3 <= 0 && 0 <= om3
            a(5) = -min([- Jinv3_2*om3^2, - Jinv3_2*oh3^2, 0]);
        else
            a(5) = -min([- Jinv3_2*om3^2, - Jinv3_2*oh3^2]);
        end

        out(1, 1) = sum(a);


        % second part: om2_dot
        % om3*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3) - om1*(Jinv3_1*om1 + Jinv3_2*om2 + Jinv3_3*om3)
        % or
        % (Jinv1_1 - Jinv3_3)*om1*om3 + Jinv2_1*om2*om3 + Jinv3_1*om3^2 - Jinv3_1*om1^2 - Jinv3_2*om1*om2
        b(1) = -min([(Jinv1_1 - Jinv3_3)*om1*om3, (Jinv1_1 - Jinv3_3)*oh1*om3, (Jinv1_1 - Jinv3_3)*om1*oh3, (Jinv1_1 - Jinv3_3)*oh1*oh3]);
        b(2) = -min([Jinv2_1*om2*om3, Jinv2_1*om2*oh3]);

        if oh3 <= 0 && 0 <= om3
            b(3) = -min([Jinv3_1*om3^2, Jinv3_1*oh3^2, 0]);
        else
            b(3) = -min([Jinv3_1*om3^2, Jinv3_1*oh3^2]);
        end
        if oh1 <= 0 && 0 <= om1
            b(4) = -min([- Jinv3_1*om1^2, - Jinv3_1*oh1^2, 0]);
        else
            b(4) = -min([- Jinv3_1*om1^2, - Jinv3_1*oh1^2]);
        end

        b(5) = -min([- Jinv3_2*om1*om2, - Jinv3_2*oh1*om2]);

        out(2, 1) = sum(b);


        % third part: om3_dot
        % om1*(Jinv2_1*om1 + Jinv2_2*om2 + Jinv3_2*om3) - om2*(Jinv1_1*om1 + Jinv2_1*om2 + Jinv3_1*om3)
        % or
        % Jinv2_1*om1^2 + (Jinv2_2-Jinv1_1)*om1*om2 + Jinv3_2*om1*om3 - Jinv2_1*om2^2 - Jinv3_1*om2*om3

        if oh1 <= 0 && 0 <= om1
            c(1) = -min([Jinv2_1*om1^2, Jinv2_1*oh1^2, 0]);
        else
            c(1) = -min([Jinv2_1*om1^2, Jinv2_1*oh1^2]);
        end
        c(2) = -min([(Jinv2_2-Jinv1_1)*om1*om2, (Jinv2_2-Jinv1_1)*oh1*om2, (Jinv2_2-Jinv1_1)*om1*oh2, (Jinv2_2-Jinv1_1)*oh1*oh2]);
        c(3) = -min([Jinv3_2*om1*om3, Jinv3_2*oh1*om3]);
        if oh2 <= 0 && 0 <= om2
            c(4) = -min([- Jinv2_1*om2^2, - Jinv2_1*oh2^2, 0]);
        else
            c(4) = -min([- Jinv2_1*om2^2, - Jinv2_1*oh2^2]);
        end
        c(5) = -min([- Jinv3_1*om2*om3, - Jinv3_1*oh2*om3]);
        out(3, 1) = sum(c);


    end

    end

    % - kp*J*eta
    function out = d_bits_3(x, xh, J)

    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    qh0 = xh(1);
    qh1 = xh(2);
    qh2 = xh(3);
    qh3 = xh(4);

    J1_1 = J(1, 1);
    J2_1 = J(2, 1);
    J3_1 = J(3, 1);
    J2_2 = J(2, 2);
    J3_2 = J(2, 3);
    J3_3 = J(3, 3);


    a = zeros(4,1);
    b = zeros(4,1);
    c = zeros(4,1);
    out = zeros(3,1);



    % - J1_1*kp*2*q0*q1 - J1_1*kp*2*q2*q3 - J2_1*kp*2*q0*q2 + J2_1*kp*2*q1*q3
    % - J2_1*kp*2*q0*q1 - J2_1*kp*2*q2*q3 - J2_2*kp*2*q0*q2 + J2_2*kp*2*q1*q3
    % - J3_1*kp*2*q0*q1 - J3_1*kp*2*q2*q3 - J3_2*kp*2*q0*q2 + J3_2*kp*2*q1*q3

    if prod(x <= xh)
        % first part: om1_dot
        % - J1_1*kp*2*q0*q1 - J1_1*kp*2*q2*q3 - J2_1*kp*2*q0*q2 + J2_1*kp*2*q1*q3

        a(1) = min([- J1_1*kp*2*q0*q1, - J1_1*kp*2*qh0*q1, - J1_1*kp*2*q0*qh1, - J1_1*kp*2*qh0*qh1]);
        a(2) = min([- J1_1*kp*2*q2*q3, - J1_1*kp*2*qh2*q3, - J1_1*kp*2*q2*qh3, - J1_1*kp*2*qh2*qh3]);
        a(3) = min([- J2_1*kp*2*q0*q2, - J2_1*kp*2*qh0*q2, - J2_1*kp*2*q0*qh2, - J2_1*kp*2*qh0*qh2]);
        a(4) = min([J2_1*kp*2*q1*q3, J2_1*kp*2*qh1*q3, J2_1*kp*2*q1*qh3, J2_1*kp*2*qh1*qh3]);
        out(1, 1) = sum(a);


        % second part: om2_dot
        % - J2_1*kp*2*q0*q1 - J2_1*kp*2*q2*q3 - J2_2*kp*2*q0*q2 + J2_2*kp*2*q1*q3
        b(1) = min([- J2_1*kp*2*q0*q1, - J2_1*kp*2*qh0*q1, - J2_1*kp*2*q0*qh1, - J2_1*kp*2*qh0*qh1]);
        b(2) = min([- J2_1*kp*2*q2*q3, - J2_1*kp*2*qh2*q3, - J2_1*kp*2*q2*qh3, - J2_1*kp*2*qh2*qh3]);
        b(3) = min([- J2_2*kp*2*q0*q2, - J2_2*kp*2*qh0*q2, - J2_2*kp*2*q0*qh2, - J2_2*kp*2*qh0*qh2]);
        b(4) = min([J2_2*kp*2*q1*q3, J2_2*kp*2*qh1*q3, J2_2*kp*2*q1*qh3, J2_2*kp*2*qh1*qh3]);

        out(2, 1) = sum(b);


        % third part: om3_dot
        % - J3_1*kp*2*q0*q1 - J3_1*kp*2*q2*q3 - J3_2*kp*2*q0*q2 + J3_2*kp*2*q1*q3
        c(1) = min([- J3_1*kp*2*q0*q1, - J3_1*kp*2*qh0*q1, - J3_1*kp*2*q0*qh1, - J3_1*kp*2*qh0*qh1]);
        c(2) = min([- J3_1*kp*2*q2*q3, - J3_1*kp*2*qh2*q3, - J3_1*kp*2*q2*qh3, - J3_1*kp*2*qh2*qh3]);
        c(3) = min([- J3_2*kp*2*q0*q2, - J3_2*kp*2*qh0*q2, - J3_2*kp*2*q0*qh2, - J3_2*kp*2*qh0*qh2]);
        c(4) = min([J3_2*kp*2*q1*q3, J3_2*kp*2*qh1*q3, J3_2*kp*2*q1*qh3, J3_2*kp*2*qh1*qh3]);

        out(3, 1) = sum(c);



    elseif prod(xh <= x)
        % first part: om1_dot
        % - J1_1*kp*2*q0*q1 - J1_1*kp*2*q2*q3 - J2_1*kp*2*q0*q2 + J2_1*kp*2*q1*q3

        a(1) = max([- J1_1*kp*2*q0*q1, - J1_1*kp*2*qh0*q1, - J1_1*kp*2*q0*qh1, - J1_1*kp*2*qh0*qh1]);
        a(2) = max([- J1_1*kp*2*q2*q3, - J1_1*kp*2*qh2*q3, - J1_1*kp*2*q2*qh3, - J1_1*kp*2*qh2*qh3]);
        a(3) = max([- J2_1*kp*2*q0*q2, - J2_1*kp*2*qh0*q2, - J2_1*kp*2*q0*qh2, - J2_1*kp*2*qh0*qh2]);
        a(4) = max([J2_1*kp*2*q1*q3, J2_1*kp*2*qh1*q3, J2_1*kp*2*q1*qh3, J2_1*kp*2*qh1*qh3]);
        out(1, 1) = sum(a);


        % second part: om2_dot
        % - J2_1*kp*2*q0*q1 - J2_1*kp*2*q2*q3 - J2_2*kp*2*q0*q2 + J2_2*kp*2*q1*q3
        b(1) = max([- J2_1*kp*2*q0*q1, - J2_1*kp*2*qh0*q1, - J2_1*kp*2*q0*qh1, - J2_1*kp*2*qh0*qh1]);
        b(2) = max([- J2_1*kp*2*q2*q3, - J2_1*kp*2*qh2*q3, - J2_1*kp*2*q2*qh3, - J2_1*kp*2*qh2*qh3]);
        b(3) = max([- J2_2*kp*2*q0*q2, - J2_2*kp*2*qh0*q2, - J2_2*kp*2*q0*qh2, - J2_2*kp*2*qh0*qh2]);
        b(4) = max([J2_2*kp*2*q1*q3, J2_2*kp*2*qh1*q3, J2_2*kp*2*q1*qh3, J2_2*kp*2*qh1*qh3]);

        out(2, 1) = sum(b);


        % third part: om3_dot
        % - J3_1*kp*2*q0*q1 - J3_1*kp*2*q2*q3 - J3_2*kp*2*q0*q2 + J3_2*kp*2*q1*q3
        c(1) = max([- J3_1*kp*2*q0*q1, - J3_1*kp*2*qh0*q1, - J3_1*kp*2*q0*qh1, - J3_1*kp*2*qh0*qh1]);
        c(2) = max([- J3_1*kp*2*q2*q3, - J3_1*kp*2*qh2*q3, - J3_1*kp*2*q2*qh3, - J3_1*kp*2*qh2*qh3]);
        c(3) = max([- J3_2*kp*2*q0*q2, - J3_2*kp*2*qh0*q2, - J3_2*kp*2*q0*qh2, - J3_2*kp*2*qh0*qh2]);
        c(4) = max([J3_2*kp*2*q1*q3, J3_2*kp*2*qh1*q3, J3_2*kp*2*q1*qh3, J3_2*kp*2*qh1*qh3]);

        out(3, 1) = sum(c);


    end

    end


    %
    function out = d_bits_4(x, xh, Jinv)
    %- Jinv1_1*om1*q1 - Jinv2_1*om1*q2
    % - Jinv2_1*om2*q1
    % - Jinv2_2*om2*q2
    % - Jinv3_1*om1*q3 - Jinv3_1*om3*q1 - Jinv3_2*om2*q3 - Jinv3_2*om3*q2 - Jinv3_3*om3*q3
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    om1 = x(5);
    om2 = x(6);
    om3 = x(7);

    qh0 = xh(1);
    qh1 = xh(2);
    qh2 = xh(3);
    qh3 = xh(4);
    oh1 = xh(5);
    oh2 = xh(6);
    oh3 = xh(7);

    Jinv1_1 = Jinv(1, 1);
    Jinv2_1 = Jinv(2, 1);
    Jinv3_1 = Jinv(3, 1);
    Jinv2_2 = Jinv(2, 2);
    Jinv3_2 = Jinv(2, 3);
    Jinv3_3 = Jinv(3, 3);

    a = zeros(9,1);
    b = zeros(9,1);
    c = zeros(9,1);
    d = zeros(9,1);
    out = zeros(4,1);


    if prod(x <= xh)
        % first part: q0_dot
        %   - Jinv1_1*om1*q1 - Jinv2_1*om1*q2 - Jinv2_1*om2*q1 - Jinv2_2*om2*q2
        %   - Jinv3_1*om1*q3 - Jinv3_1*om3*q1 - Jinv3_2*om2*q3 - Jinv3_2*om3*q2 - Jinv3_3*om3*q3

        a(1) = min([ - Jinv1_1*om1*q1, - Jinv1_1*oh1*q1, - Jinv1_1*om1*qh1, - Jinv1_1*oh1*qh1]);
        a(2) = min([ - Jinv2_1*om1*q2, - Jinv2_1*oh1*q2, - Jinv2_1*om1*qh2, - Jinv2_1*oh1*qh2]);
        a(3) = min([ - Jinv2_1*om2*q1, - Jinv2_1*oh2*q1, - Jinv2_1*om2*qh1, - Jinv2_1*oh2*qh1]);
        a(4) = min([ - Jinv2_2*om2*q2, - Jinv2_2*oh2*q2, - Jinv2_2*om2*qh2, - Jinv2_2*oh2*qh2]);
        a(5) = min([ - Jinv3_1*om1*q3, - Jinv3_1*oh1*q3, - Jinv3_1*om1*qh3, - Jinv3_1*oh1*qh3]);
        a(6) = min([ - Jinv3_1*om3*q1, - Jinv3_1*oh3*q1, - Jinv3_1*om3*qh1, - Jinv3_1*oh3*qh1]);
        a(7) = min([ - Jinv3_2*om2*q3, - Jinv3_2*oh2*q3, - Jinv3_2*om2*qh3, - Jinv3_2*oh2*qh3]);
        a(8) = min([ - Jinv3_2*om3*q2, - Jinv3_2*oh3*q2, - Jinv3_2*om3*qh2, - Jinv3_2*oh3*qh2]);
        a(9) = min([ - Jinv3_3*om3*q3, - Jinv3_3*oh3*q3, - Jinv3_3*om3*qh3, - Jinv3_3*oh3*qh3]);

        out(1, 1) = 1/2*sum(a);


        % second part: q1_dot
        %     Jinv1_1*om1*q0 + Jinv2_1*om2*q0 - Jinv2_1*om1*q3 - Jinv2_2*om2*q3
        %   + Jinv3_1*om1*q2 + Jinv3_1*om3*q0 + Jinv3_2*om2*q2 - Jinv3_2*om3*q3 + Jinv3_3*om3*q2
        b(1) = min([  Jinv1_1*om1*q0,   Jinv1_1*oh1*q0,   Jinv1_1*om1*qh0,   Jinv1_1*oh1*qh0]);
        b(2) = min([  Jinv2_1*om2*q0,   Jinv2_1*oh2*q0,   Jinv2_1*om2*qh0,   Jinv2_1*oh2*qh0]);
        b(3) = min([- Jinv2_1*om1*q3, - Jinv2_1*oh1*q3, - Jinv2_1*om1*qh3, - Jinv2_1*oh1*qh3]);
        b(4) = min([- Jinv2_2*om2*q3, - Jinv2_2*oh2*q3, - Jinv2_2*om2*qh3, - Jinv2_2*oh2*qh3]);
        b(5) = min([  Jinv3_1*om1*q2,   Jinv3_1*oh1*q2,   Jinv3_1*om1*qh2,   Jinv3_1*oh1*qh2]);
        b(6) = min([  Jinv3_1*om3*q0,   Jinv3_1*oh3*q0,   Jinv3_1*om3*qh0,   Jinv3_1*oh3*qh0]);
        b(7) = min([  Jinv3_2*om2*q2,   Jinv3_2*oh2*q2,   Jinv3_2*om2*qh2,   Jinv3_2*oh2*qh2]);
        b(8) = min([- Jinv3_2*om3*q3, - Jinv3_2*oh3*q3, - Jinv3_2*om3*qh3, - Jinv3_2*oh3*qh3]);
        b(9) = min([  Jinv3_3*om3*q2,   Jinv3_3*oh3*q2,   Jinv3_3*om3*qh2,   Jinv3_3*oh3*qh2]);

        out(2, 1) = 1/2*sum(b);


        % third part: q2_dot
        %     Jinv1_1*om1*q3 + Jinv2_1*om1*q0 + Jinv2_2*om2*q0 + Jinv2_1*om2*q3
        %   - Jinv3_1*om1*q1 - Jinv3_2*om2*q1 + Jinv3_2*om3*q0 + Jinv3_1*om3*q3 - Jinv3_3*om3*q1
        c(1) = min([  Jinv1_1*om1*q3,   Jinv1_1*oh1*q3,   Jinv1_1*om1*qh3,   Jinv1_1*oh1*qh3]);
        c(2) = min([  Jinv2_1*om1*q0,   Jinv2_1*oh1*q0,   Jinv2_1*om1*qh0,   Jinv2_1*oh1*qh0]);
        c(3) = min([  Jinv2_2*om2*q0,   Jinv2_2*oh2*q0,   Jinv2_2*om2*qh0,   Jinv2_2*oh2*qh0]);
        c(4) = min([  Jinv2_1*om2*q3,   Jinv2_1*oh2*q3,   Jinv2_1*om2*qh3,   Jinv2_1*oh2*qh3]);
        c(5) = min([- Jinv3_1*om1*q1, - Jinv3_1*oh1*q1, - Jinv3_1*om1*qh1, - Jinv3_1*oh1*qh1]);
        c(6) = min([- Jinv3_2*om2*q1, - Jinv3_2*oh2*q1, - Jinv3_2*om2*qh1, - Jinv3_2*oh2*qh1]);
        c(7) = min([  Jinv3_2*om3*q0,   Jinv3_2*oh3*q0,   Jinv3_2*om3*qh0,   Jinv3_2*oh3*qh0]);
        c(8) = min([  Jinv3_1*om3*q3,   Jinv3_1*oh3*q3,   Jinv3_1*om3*qh3,   Jinv3_1*oh3*qh3]);
        c(9) = min([- Jinv3_3*om3*q1, - Jinv3_3*oh3*q1, - Jinv3_3*om3*qh1, - Jinv3_3*oh3*qh1]);

        out(3, 1) = 1/2*sum(c);


        % fourth part: q3_dot
        %    Jinv2_1*om1*q1 - Jinv1_1*om1*q2 - Jinv2_1*om2*q2 + Jinv2_2*om2*q1
        %  + Jinv3_1*om1*q0 + Jinv3_2*om2*q0 - Jinv3_1*om3*q2 + Jinv3_2*om3*q1 + Jinv3_3*om3*q0
        d(1) = min([  Jinv2_1*om1*q1,   Jinv2_1*oh1*q1,   Jinv2_1*om1*qh1,   Jinv2_1*oh1*qh1]);
        d(2) = min([- Jinv1_1*om1*q2, - Jinv1_1*oh1*q2, - Jinv1_1*om1*qh2, - Jinv1_1*oh1*qh2]);
        d(3) = min([- Jinv2_1*om2*q2, - Jinv2_1*oh2*q2, - Jinv2_1*om2*qh2, - Jinv2_1*oh2*qh2]);
        d(4) = min([  Jinv2_2*om2*q1,   Jinv2_2*oh2*q1,   Jinv2_2*om2*qh1,   Jinv2_2*oh2*qh1]);
        d(5) = min([  Jinv3_1*om1*q0,   Jinv3_1*oh1*q0,   Jinv3_1*om1*qh0,   Jinv3_1*oh1*qh0]);
        d(6) = min([  Jinv3_2*om2*q0,   Jinv3_2*oh2*q0,   Jinv3_2*om2*qh0,   Jinv3_2*oh2*qh0]);
        d(7) = min([- Jinv3_1*om3*q2, - Jinv3_1*oh3*q2, - Jinv3_1*om3*qh2, - Jinv3_1*oh3*qh2]);
        d(8) = min([  Jinv3_2*om3*q1,   Jinv3_2*oh3*q1,   Jinv3_2*om3*qh1,   Jinv3_2*oh3*qh1]);
        d(9) = min([  Jinv3_3*om3*q0,   Jinv3_3*oh3*q0,   Jinv3_3*om3*qh0,   Jinv3_3*oh3*qh0]);

        out(4, 1) = 1/2*sum(d);

    elseif prod(xh <= x)
        % first part: q0_dot
        %   - Jinv1_1*om1*q1 - Jinv2_1*om1*q2 - Jinv2_1*om2*q1 - Jinv2_2*om2*q2
        %   - Jinv3_1*om1*q3 - Jinv3_1*om3*q1 - Jinv3_2*om2*q3 - Jinv3_2*om3*q2 - Jinv3_3*om3*q3

        a(1) = max([ - Jinv1_1*om1*q1, - Jinv1_1*oh1*q1, - Jinv1_1*om1*qh1, - Jinv1_1*oh1*qh1]);
        a(2) = max([ - Jinv2_1*om1*q2, - Jinv2_1*oh1*q2, - Jinv2_1*om1*qh2, - Jinv2_1*oh1*qh2]);
        a(3) = max([ - Jinv2_1*om2*q1, - Jinv2_1*oh2*q1, - Jinv2_1*om2*qh1, - Jinv2_1*oh2*qh1]);
        a(4) = max([ - Jinv2_2*om2*q2, - Jinv2_2*oh2*q2, - Jinv2_2*om2*qh2, - Jinv2_2*oh2*qh2]);
        a(5) = max([ - Jinv3_1*om1*q3, - Jinv3_1*oh1*q3, - Jinv3_1*om1*qh3, - Jinv3_1*oh1*qh3]);
        a(6) = max([ - Jinv3_1*om3*q1, - Jinv3_1*oh3*q1, - Jinv3_1*om3*qh1, - Jinv3_1*oh3*qh1]);
        a(7) = max([ - Jinv3_2*om2*q3, - Jinv3_2*oh2*q3, - Jinv3_2*om2*qh3, - Jinv3_2*oh2*qh3]);
        a(8) = max([ - Jinv3_2*om3*q2, - Jinv3_2*oh3*q2, - Jinv3_2*om3*qh2, - Jinv3_2*oh3*qh2]);
        a(9) = max([ - Jinv3_3*om3*q3, - Jinv3_3*oh3*q3, - Jinv3_3*om3*qh3, - Jinv3_3*oh3*qh3]);

        out(1, 1) = 1/2*sum(a);


        % second part: q1_dot
        %     Jinv1_1*om1*q0 + Jinv2_1*om2*q0 - Jinv2_1*om1*q3 - Jinv2_2*om2*q3
        %   + Jinv3_1*om1*q2 + Jinv3_1*om3*q0 + Jinv3_2*om2*q2 - Jinv3_2*om3*q3 + Jinv3_3*om3*q2
        b(1) = max([  Jinv1_1*om1*q0,   Jinv1_1*oh1*q0,   Jinv1_1*om1*qh0,   Jinv1_1*oh1*qh0]);
        b(2) = max([  Jinv2_1*om2*q0,   Jinv2_1*oh2*q0,   Jinv2_1*om2*qh0,   Jinv2_1*oh2*qh0]);
        b(3) = max([- Jinv2_1*om1*q3, - Jinv2_1*oh1*q3, - Jinv2_1*om1*qh3, - Jinv2_1*oh1*qh3]);
        b(4) = max([- Jinv2_2*om2*q3, - Jinv2_2*oh2*q3, - Jinv2_2*om2*qh3, - Jinv2_2*oh2*qh3]);
        b(5) = max([  Jinv3_1*om1*q2,   Jinv3_1*oh1*q2,   Jinv3_1*om1*qh2,   Jinv3_1*oh1*qh2]);
        b(6) = max([  Jinv3_1*om3*q0,   Jinv3_1*oh3*q0,   Jinv3_1*om3*qh0,   Jinv3_1*oh3*qh0]);
        b(7) = max([  Jinv3_2*om2*q2,   Jinv3_2*oh2*q2,   Jinv3_2*om2*qh2,   Jinv3_2*oh2*qh2]);
        b(8) = max([- Jinv3_2*om3*q3, - Jinv3_2*oh3*q3, - Jinv3_2*om3*qh3, - Jinv3_2*oh3*qh3]);
        b(9) = max([  Jinv3_3*om3*q2,   Jinv3_3*oh3*q2,   Jinv3_3*om3*qh2,   Jinv3_3*oh3*qh2]);

        out(2, 1) = 1/2*sum(b);


        % third part: q2_dot
        %     Jinv1_1*om1*q3 + Jinv2_1*om1*q0 + Jinv2_2*om2*q0 + Jinv2_1*om2*q3
        %   - Jinv3_1*om1*q1 - Jinv3_2*om2*q1 + Jinv3_2*om3*q0 + Jinv3_1*om3*q3 - Jinv3_3*om3*q1
        c(1) = max([  Jinv1_1*om1*q3,   Jinv1_1*oh1*q3,   Jinv1_1*om1*qh3,   Jinv1_1*oh1*qh3]);
        c(2) = max([  Jinv2_1*om1*q0,   Jinv2_1*oh1*q0,   Jinv2_1*om1*qh0,   Jinv2_1*oh1*qh0]);
        c(3) = max([  Jinv2_2*om2*q0,   Jinv2_2*oh2*q0,   Jinv2_2*om2*qh0,   Jinv2_2*oh2*qh0]);
        c(4) = max([  Jinv2_1*om2*q3,   Jinv2_1*oh2*q3,   Jinv2_1*om2*qh3,   Jinv2_1*oh2*qh3]);
        c(5) = max([- Jinv3_1*om1*q1, - Jinv3_1*oh1*q1, - Jinv3_1*om1*qh1, - Jinv3_1*oh1*qh1]);
        c(6) = max([- Jinv3_2*om2*q1, - Jinv3_2*oh2*q1, - Jinv3_2*om2*qh1, - Jinv3_2*oh2*qh1]);
        c(7) = max([  Jinv3_2*om3*q0,   Jinv3_2*oh3*q0,   Jinv3_2*om3*qh0,   Jinv3_2*oh3*qh0]);
        c(8) = max([  Jinv3_1*om3*q3,   Jinv3_1*oh3*q3,   Jinv3_1*om3*qh3,   Jinv3_1*oh3*qh3]);
        c(9) = max([- Jinv3_3*om3*q1, - Jinv3_3*oh3*q1, - Jinv3_3*om3*qh1, - Jinv3_3*oh3*qh1]);

        out(3, 1) = 1/2*sum(c);


        % fourth part: q3_dot
        %    Jinv2_1*om1*q1 - Jinv1_1*om1*q2 - Jinv2_1*om2*q2 + Jinv2_2*om2*q1
        %  + Jinv3_1*om1*q0 + Jinv3_2*om2*q0 - Jinv3_1*om3*q2 + Jinv3_2*om3*q1 + Jinv3_3*om3*q0
        d(1) = max([  Jinv2_1*om1*q1,   Jinv2_1*oh1*q1,   Jinv2_1*om1*qh1,   Jinv2_1*oh1*qh1]);
        d(2) = max([- Jinv1_1*om1*q2, - Jinv1_1*oh1*q2, - Jinv1_1*om1*qh2, - Jinv1_1*oh1*qh2]);
        d(3) = max([- Jinv2_1*om2*q2, - Jinv2_1*oh2*q2, - Jinv2_1*om2*qh2, - Jinv2_1*oh2*qh2]);
        d(4) = max([  Jinv2_2*om2*q1,   Jinv2_2*oh2*q1,   Jinv2_2*om2*qh1,   Jinv2_2*oh2*qh1]);
        d(5) = max([  Jinv3_1*om1*q0,   Jinv3_1*oh1*q0,   Jinv3_1*om1*qh0,   Jinv3_1*oh1*qh0]);
        d(6) = max([  Jinv3_2*om2*q0,   Jinv3_2*oh2*q0,   Jinv3_2*om2*qh0,   Jinv3_2*oh2*qh0]);
        d(7) = max([- Jinv3_1*om3*q2, - Jinv3_1*oh3*q2, - Jinv3_1*om3*qh2, - Jinv3_1*oh3*qh2]);
        d(8) = max([  Jinv3_2*om3*q1,   Jinv3_2*oh3*q1,   Jinv3_2*om3*qh1,   Jinv3_2*oh3*qh1]);
        d(9) = max([  Jinv3_3*om3*q0,   Jinv3_3*oh3*q0,   Jinv3_3*om3*qh0,   Jinv3_3*oh3*qh0]);

        out(4, 1) = 1/2*sum(d);



    end




    end


end