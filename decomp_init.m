

% decomposition function for closed-loop dynamics (without linear transformation)


function out = decomp_init(x, xh, w, wh, J, Jinv, umax, kp, kd)
    J1  = J(1, 1);
    J2  = J(2, 2);
    J3  = J(3, 3);
    J12 = J(1, 2);
    J13 = J(1, 3);
    J23 = J(2, 3);
    Ji1  = Jinv(1, 1);
    Ji2  = Jinv(2, 2);
    Ji3  = Jinv(3, 3);
    Ji12 = Jinv(1, 2);
    Ji13 = Jinv(1, 3);
    Ji23 = Jinv(2, 3);

    out = [quat_decomp(x, w, xh, wh); ...
           dec_omg1(x, xh, w, wh, J, Jinv) + thing1(x, xh, J, Jinv); ...
           dec_omg2(x, xh, w, wh, J, Jinv) + thing2(x, xh, J, Jinv);...
           dec_omg3(x, xh, w, wh, J, Jinv) + thing3(x, xh, J, Jinv)];
       
    % Saturation function
    function out = sigma(in)
        out = umax*tanh((1/umax)*in);
    end

    % decomposes -inv(J)*cross(omg, J*omg), where omg = x(5:7)
    function out = dec_omg1(x, xh, w, wh, J, Jinv)

        a11 = J13*Ji12 - J12*Ji13;
        a22 = J12*Ji13 - J23*Ji1;
        a33 = J23*Ji1 - J13*Ji12;
        a12 = J23*Ji12 - J13*Ji1 + Ji13*(J1 - J2);
        a13 = J12*Ji1 - J23*Ji13 - Ji12*(J1 - J3);
        a23 = J13*Ji13 - J12*Ji12 + Ji1*(J2 - J3);

        if prod(x <= xh)
            b1 = min([a12*x(5)*x(6), a12*x(5)*xh(6)]);
            b2 = min([a13*x(5)*x(7), a13*x(5)*xh(7)]);
            b3 = min([a23*x(6)*x(7), a23*xh(6)*x(7), a23*x(6)*xh(7), a23*xh(6)*xh(7)]);
            b4 = a11*x(5)^2;

            if x(6) <= 0 && 0 <= xh(6)
                b5 = min([a22*x(6)^2, a22*xh(6)^2, 0]);
            else
                b5 = min([a22*x(6)^2, a22*xh(6)^2]);
            end
            if x(7) <= 0 && 0 <= xh(7)
                b6 = min([a33*x(7)^2, a33*xh(7)^2, 0]);
            else
                b6 = min([a33*x(7)^2, a33*xh(7)^2]);
            end

            out = b1+b2+b3+b4+b5+b6;

        elseif prod(xh <= x)
            b1 = max([a12*x(5)*x(6), a12*x(5)*xh(6)]);
            b2 = max([a13*x(5)*x(7), a13*x(5)*xh(7)]);
            b3 = max([a23*x(6)*x(7), a23*xh(6)*x(7), a23*x(6)*xh(7), a23*xh(6)*xh(7)]);
            b4 = a11*x(5)^2;

            if x(6) <= 0 && 0 <= xh(6)
                b5 = max([a22*x(6)^2, a22*xh(6)^2, 0]);
            else
                b5 = max([a22*x(6)^2, a22*xh(6)^2]);
            end
            if x(7) <= 0 && 0 <= xh(7)
                b6 = max([a33*x(7)^2, a33*xh(7)^2, 0]);
            else
                b6 = max([a33*x(7)^2, a33*xh(7)^2]);
            end

            out = b1+b2+b3+b4+b5+b6;

        else
            out = nan; % bad
        end

        if Jinv(1, 1) >= 0
            c1 = Jinv(1, 1)*w(1);
        else
            c1 = Jinv(1, 1)*wh(1);
        end
        if Jinv(1, 2) >= 0
            c2 = Jinv(1, 2)*w(2);
        else
            c2 = Jinv(1, 2)*wh(2);
        end
        if Jinv(1, 3) >= 0
            c3 = Jinv(1, 3)*w(3);
        else
            c3 = Jinv(1, 3)*wh(3);
        end

        out = out + c1 + c2 + c3;
    end

    % decomposes -inv(J)*cross(omg, J*omg), where omg = x(5:7)
    function out = dec_omg2(x, xh, w, wh, J, Jinv)


        a11 = J13*Ji2 - J12*Ji23;
        a22 = J12*Ji23 - J23*Ji12;
        a33 = - J13*Ji2+ J23*Ji12;
        a13 = J3*Ji2- J1*Ji2+ J12*Ji12- J23*Ji23;
        a23 = J2*Ji12- J12*Ji2- J3*Ji12+ J13*Ji23;
        a12 = J1*Ji23 - J2*Ji23- J13*Ji12+ J23*Ji2;

        if prod(x <= xh)
            b1 = min([a12*x(5)*x(6), a12*xh(5)*x(6)]);
            b2 = min([a13*x(5)*x(7), a13*xh(5)*x(7), a13*x(5)*xh(7), a13*xh(5)*xh(7)]);
            b3 = min([a23*x(6)*x(7), a23*x(6)*xh(7)]);

            b4 = a22*x(6)^2;

            if x(5) <= 0 && 0 <= xh(5)
                b5 = min([a11*x(5)^2, a11*xh(5)^2, 0]);
            else
                b5 = min([a11*x(5)^2, a11*xh(5)^2]);
            end
            if x(7) <= 0 && 0 <= xh(7)
                b6 = min([a33*x(7)^2, a33*xh(7)^2, 0]);
            else
                b6 = min([a33*x(7)^2, a33*xh(7)^2]);
            end

            out = (b1+b2+b3+b4+b5+b6);

        elseif prod(xh <= x)
            b1 = max([a12*x(5)*x(6), a12*xh(5)*x(6)]);
            b2 = max([a13*x(5)*x(7), a13*xh(5)*x(7), a13*x(5)*xh(7), a13*xh(5)*xh(7)]);
            b3 = max([a23*x(6)*x(7), a23*x(6)*xh(7)]);

            b4 = a22*x(6)^2;

            if x(5) <= 0 && 0 <= xh(5)
                b5 = max([a11*x(5)^2, a11*xh(5)^2, 0]);
            else
                b5 = max([a11*x(5)^2, a11*xh(5)^2]);
            end
            if x(7) <= 0 && 0 <= xh(7)
                b6 = max([a33*x(7)^2, a33*xh(7)^2, 0]);
            else
                b6 = max([a33*x(7)^2, a33*xh(7)^2]);
            end

            out = (b1+b2+b3+b4+b5+b6);
        else
            out = nan;
        end

        if Jinv(2, 1) >= 0
            c1 = Jinv(2, 1)*w(1);
        else
            c1 = Jinv(2, 1)*wh(1);
        end
        if Jinv(2, 2) >= 0
            c2 = Jinv(2, 2)*w(2);
        else
            c2 = Jinv(2, 2)*wh(2);
        end
        if Jinv(2, 3) >= 0
            c3 = Jinv(2, 3)*w(3);
        else
            c3 = Jinv(2, 3)*wh(3);
        end

        out = out + c1 + c2 + c3;
    end

    % decomposes -inv(J)*cross(omg, J*omg), where omg = x(5:7)
    function out = dec_omg3(x, xh, w, wh, J, Jinv)

        a11 = J13*Ji23- J12*Ji3;
        a22 = J12*Ji3- J23*Ji13;
        a33 = J23*Ji13 - J13*Ji23;
        a13 = - J1*Ji23 + J12*Ji13+ J3*Ji23- J23*Ji3;
        a23 = J2*Ji13- J3*Ji13 + J13*Ji3- J12*Ji23;
        a12 = J1*Ji3- J2*Ji3- J13*Ji13+ J23*Ji23;

        if prod(x <= xh)
            b1 = min([a12*x(5)*x(6), a12*xh(5)*x(6), a12*x(5)*xh(6), a12*xh(5)*xh(6)]);
            b2 = min([a13*x(5)*x(7), a13*xh(5)*x(7)]);
            b3 = min([a23*x(6)*x(7), a23*xh(6)*x(7)]);

            b4 = a33*x(7)^2;

            if x(5) <= 0 && 0 <= xh(5)
                b5 = min([a11*x(5)^2, a11*xh(5)^2, 0]);
            else
                b5 = min([a11*x(5)^2, a11*xh(5)^2]);
            end
            if x(6) <= 0 && 0 <= xh(6)
                b6 = min([a22*x(6)^2, a22*xh(6)^2, 0]);
            else
                b6 = min([a22*x(6)^2, a22*xh(6)^2]);
            end

            out = (b1+b2+b3+b4+b5+b6);

        elseif prod(xh <= x)
            b1 = max([a12*x(5)*x(6), a12*xh(5)*x(6), a12*x(5)*xh(6), a12*xh(5)*xh(6)]);
            b2 = max([a13*x(5)*x(7), a13*xh(5)*x(7)]);
            b3 = max([a23*x(6)*x(7), a23*xh(6)*x(7)]);

            b4 = a33*x(7)^2;

            if x(5) <= 0 && 0 <= xh(5)
                b5 = max([a11*x(5)^2, a11*xh(5)^2, 0]);
            else
                b5 = max([a11*x(5)^2, a11*xh(5)^2]);
            end
            if x(6) <= 0 && 0 <= xh(6)
                b6 = max([a22*x(6)^2, a22*xh(6)^2, 0]);
            else
                b6 = max([a22*x(6)^2, a22*xh(6)^2]);
            end

            out = (b1+b2+b3+b4+b5+b6);
        else
            out = nan;
        end

        if Jinv(3, 1) >= 0
            c1 = Jinv(3, 1)*w(1);
        else
            c1 = Jinv(3, 1)*wh(1);
        end
        if Jinv(3, 2) >= 0
            c2 = Jinv(3, 2)*w(2);
        else
            c2 = Jinv(3, 2)*wh(2);
        end
        if Jinv(3, 3) >= 0
            c3 = Jinv(3, 3)*w(3);
        else
            c3 = Jinv(3, 3)*wh(3);
        end

        out = out + c1 + c2 + c3;
    end



    function out = thing1(x, xh, J, Jinv)

        q0 = x(1);
        q1 = x(2);
        q2 = x(3);
        q3 = x(4);

        qh0 = xh(1);
        qh1 = xh(2);
        qh2 = xh(3);
        qh3 = xh(4);

        om1 = x(5);
        om2 = x(6);
        om3 = x(7);

        oh2 = xh(6);
        oh3 = xh(7);

        if prod(x <= xh)        
            b = -sign(Ji12);

            a1 = min([b*(J3-J1)*om1*om3, b*(J3-J1)*om1*oh3]);
            a2 = min([-b*J12*om2*om3, -b*J12*oh2*om3, -b*J12*om2*oh3, -b*J12*oh2*oh3]);
            if om3 <= 0  && 0 <= oh3
                a3 = min([-b*J13*om3^2, -b*J13*oh3^2, 0]);
            else
                a3 = min([-b*J13*om3^2, -b*J13*oh3^2]);
            end
            a4 =  b*J13*om1^2;
            a5 = min([ b*J23*om1*om2,  b*J23*om1*oh2]);

            
            a6 = b*J12*kd*om1;
            a7 = min([b*J2*kd*om2, b*J2*kd*oh2]);
            a8 = min([b*J23*kd*om3, b*J23*kd*oh3]);


            a9 = min([b*J12*kp*2*q0*q1, b*J12*kp*2*qh0*q1, b*J12*kp*2*q0*qh1,  b*J12*kp*2*qh0*qh1]);
            a10 = min([b*J12*kp*2*q2*q3, b*J12*kp*2*qh2*q3,  b*J12*kp*2*q2*qh3,  b*J12*kp*2*qh2*qh3]);
            a11 = min([J2*b*kp*2*q0*q2, J2*b*kp*2*qh0*q2,  J2*b*kp*2*q0*qh2,  J2*b*kp*2*qh0*qh2]);
            a12 = min([-J2*b*kp*2*q1*q3, -J2*b*kp*2*qh1*q3, -J2*b*kp*2*q1*qh3, -J2*b*kp*2*qh1*qh3]);

            ZZ1 = umax*sign(Ji12)*Ji12*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works





            %%%part 2

            b = -sign(Ji1);

            a1 = min([b*J12*om1*om3, b*J12*om1*oh3]);
            a2 = min([b*(J2- J3)*om2*om3, b*(J2- J3)*oh2*om3, b*(J2 - J3)*om2*oh3, b*(J2- J3)*oh2*oh3]);
            if om3 <= 0  && 0 <= oh3
                a3 = min([ b*J23*om3^2, b*J23*oh3^2, 0]);
            else
                a3 = min([ b*J23*om3^2, b*J23*oh3^2]);
            end
            a4 = min([-J13*b*om1*om2, -J13*b*om1*oh2]);
            if om2 <= 0  && 0 <= oh2
                a5 = min([-J23*b*om2^2, -J23*b*oh2^2, 0]);
            else
                a5 = min([-J23*b*om2^2, -J23*b*oh2^2]);
            end

            b = -b;
            a6 = -b*J12*kd*om1;
            a7 = min([-b*J2*kd*om2, -b*J2*kd*oh2]);
            a8 = min([- b*J23*kd*om3, - b*J23*kd*oh3]);

            a9 = min([- b*J1*kp*2*q0*q1, - b*J1*kp*2*qh0*q1, - b*J1*kp*2*q0*qh1, - b*J1*kp*2*qh0*qh1]);
            a10 = min([- b*J1*kp*2*q2*q3, - b*J1*kp*2*qh2*q3, - b*J1*kp*2*q2*qh3, - b*J1*kp*2*qh2*qh3]);
            a11 = min([- J12*b*kp*2*q0*q2, - J12*b*kp*2*qh0*q2, - J12*b*kp*2*q0*qh2, - J12*b*kp*2*qh0*qh2]);
            a12 = min([J12*b*kp*2*q1*q3, J12*b*kp*2*qh1*q3, J12*b*kp*2*q1*qh3, J12*b*kp*2*qh1*qh3]);

            ZZ2 = umax*sign(Ji1)*Ji1*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works


            %%% Part 3
            b = -sign(Ji13);
            a1 = min([b*(J1-J2)*om1*om2, b*(J1-J2)*om1*oh2]);
            if om2 <= 0  && 0 <= oh2
                a2 = min([b*J12*om2^2, b*J12*oh2^2, 0]);
            else
                a2 = min([b*J12*om2^2, b*J12*oh2^2]);
            end
            a3 = min([b*J13*om2*om3, b*J13*oh2*om3, b*J13*om2*oh3, b*J13*oh2*oh3]);
            a4 = -b*J12*om1^2;
            a5 = min([-b*J23*om1*om3, -b*J23*om1*oh3]);

            b = -b;
            a6 = -b*J12*kd*om1;
            a7 = min([-b*J2*kd*om2, -b*J2*kd*oh2]);
            a8 = min([- b*J23*kd*om3, - b*J23*kd*oh3]);

            a9 = min([- b*J13*kp*2*q0*q1, - b*J13*kp*2*qh0*q1, - b*J13*kp*2*q0*qh1, - b*J13*kp*2*qh0*qh1]);
            a10 = min([- b*J13*kp*2*q2*q3, - b*J13*kp*2*qh2*q3, - b*J13*kp*2*q2*qh3, - b*J13*kp*2*qh2*qh3]);
            a11 = min([- J23*b*kp*2*q0*q2, - J23*b*kp*2*qh0*q2, - J23*b*kp*2*q0*qh2, - J23*b*kp*2*qh0*qh2]);
            a12 = min([J23*b*kp*2*q1*q3, J23*b*kp*2*qh1*q3, J23*b*kp*2*q1*qh3, J23*b*kp*2*qh1*qh3]);

            ZZ3 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works

            out = ZZ1 + ZZ2 + ZZ3;

        elseif prod(xh <= x)
            b = -sign(Ji12);

            a1 = max([b*(J3-J1)*om1*om3, b*(J3-J1)*om1*oh3]);
            a2 = max([-b*J12*om2*om3, -b*J12*oh2*om3, -b*J12*om2*oh3, -b*J12*oh2*oh3]);
            if oh3 <= 0  && 0 <= om3
                a3 = max([-b*J13*om3^2, -b*J13*oh3^2, 0]);
            else
                a3 = max([-b*J13*om3^2, -b*J13*oh3^2]);
            end
            a4 =  b*J13*om1^2;
            a5 = max([ b*J23*om1*om2,  b*J23*om1*oh2]);

            b = -b;
            a6 = -b*J12*kd*om1;
            a7 = max([-b*J2*kd*om2, -b*J2*kd*oh2]);
            a8 = max([- b*J23*kd*om3, - b*J23*kd*oh3]);

            a9 = max([- b*J12*kp*2*q0*q1, - b*J12*kp*2*qh0*q1, - b*J12*kp*2*q0*qh1, - b*J12*kp*2*qh0*qh1]);
            a10 = max([- b*J12*kp*2*q2*q3, - b*J12*kp*2*qh2*q3, - b*J12*kp*2*q2*qh3, - b*J12*kp*2*qh2*qh3]);
            a11 = max([- J2*b*kp*2*q0*q2, - J2*b*kp*2*qh0*q2, - J2*b*kp*2*q0*qh2, - J2*b*kp*2*qh0*qh2]);
            a12 = max([J2*b*kp*2*q1*q3, J2*b*kp*2*qh1*q3, J2*b*kp*2*q1*qh3, J2*b*kp*2*qh1*qh3]);

            ZZ1 = umax*sign(Ji12)*Ji12*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works


            %%%
            b = -sign(Ji1);

            a1 = max([b*J12*om1*om3, b*J12*om1*oh3]);
            a2 = max([b*(J2- J3)*om2*om3, b*(J2- J3)*oh2*om3, b*(J2 - J3)*om2*oh3, b*(J2- J3)*oh2*oh3]);
            if oh3 <= 0  && 0 <= om3
                a3 = max([ b*J23*om3^2, b*J23*oh3^2, 0]);
            else
                a3 = max([ b*J23*om3^2, b*J23*oh3^2]);
            end
            a4 = max([-J13*b*om1*om2, -J13*b*om1*oh2]);
            if oh2 <= 0  && 0 <= om2
                a5 = max([-J23*b*om2^2, -J23*b*oh2^2, 0]);
            else
                a5 = max([-J23*b*om2^2, -J23*b*oh2^2]);
            end

            b = -b;
            a6 = -b*J12*kd*om1;
            a7 = max([-b*J2*kd*om2, -b*J2*kd*oh2]);
            a8 = max([- b*J23*kd*om3, - b*J23*kd*oh3]);

            a9 = max([- b*J1*kp*2*q0*q1, - b*J1*kp*2*qh0*q1, - b*J1*kp*2*q0*qh1, - b*J1*kp*2*qh0*qh1]);
            a10 = max([- b*J1*kp*2*q2*q3, - b*J1*kp*2*qh2*q3, - b*J1*kp*2*q2*qh3, - b*J1*kp*2*qh2*qh3]);
            a11 = max([- J12*b*kp*2*q0*q2, - J12*b*kp*2*qh0*q2, - J12*b*kp*2*q0*qh2, - J12*b*kp*2*qh0*qh2]);
            a12 = max([J12*b*kp*2*q1*q3, J12*b*kp*2*qh1*q3, J12*b*kp*2*q1*qh3, J12*b*kp*2*qh1*qh3]);

            ZZ2 = umax*sign(Ji1)*Ji1*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works



            %%%% Part
            b = -sign(Ji13);
            a1 = max([b*(J1-J2)*om1*om2, b*(J1-J2)*om1*oh2]);
            if oh2 <= 0  && 0 <= om2
                a2 = max([b*J12*om2^2, b*J12*oh2^2, 0]);
            else
                a2 = max([b*J12*om2^2, b*J12*oh2^2]);
            end
            a3 = max([b*J13*om2*om3, b*J13*oh2*om3, b*J13*om2*oh3, b*J13*oh2*oh3]);
            a4 = -b*J12*om1^2;
            a5 = max([-b*J23*om1*om3, -b*J23*om1*oh3]);

            b = -b;
            a6 = -b*J12*kd*om1;
            a7 = max([-b*J2*kd*om2, -b*J2*kd*oh2]);
            a8 = max([- b*J23*kd*om3, - b*J23*kd*oh3]);

            a9 = max([- b*J13*kp*2*q0*q1, - b*J13*kp*2*qh0*q1, - b*J13*kp*2*q0*qh1, - b*J13*kp*2*qh0*qh1]);
            a10 = max([- b*J13*kp*2*q2*q3, - b*J13*kp*2*qh2*q3, - b*J13*kp*2*q2*qh3, - b*J13*kp*2*qh2*qh3]);
            a11 = max([- J23*b*kp*2*q0*q2, - J23*b*kp*2*qh0*q2, - J23*b*kp*2*q0*qh2, - J23*b*kp*2*qh0*qh2]);
            a12 = max([J23*b*kp*2*q1*q3, J23*b*kp*2*qh1*q3, J23*b*kp*2*q1*qh3, J23*b*kp*2*qh1*qh3]);

            ZZ3 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works


            out = ZZ1 + ZZ2 + ZZ3;
        else
            out = nan;
        end

    end

    function out = thing2(x, xh, J, Jinv)

        q0 = x(1);
        q1 = x(2);
        q2 = x(3);
        q3 = x(4);

        qh0 = xh(1);
        qh1 = xh(2);
        qh2 = xh(3);
        qh3 = xh(4);

        om1 = x(5);
        om2 = x(6);
        om3 = x(7);

        oh1 = xh(5);
        oh2 = xh(6);
        oh3 = xh(7);

        if prod(x <= xh)

            b = -sign(Ji2);

            a1 = min([b*(J3-J1)*om1*om3, b*(J3-J1)*oh1*om3, b*(J3-J1)*om1*oh3, b*(J3-J1)*oh1*oh3]);
            a2 = min([- b*J12*om2*om3, - b*J12*om2*oh3]);
            if om3 <= 0  && 0 <= oh3
                a3 = min([- b*J13*om3^2, - b*J13*oh3^2, 0]);
            else
                a3 = min([- b*J13*om3^2, - b*J13*oh3^2]);
            end
            if om1 <= 0  && 0 <= oh1
                a4 = min([b*J13*om1^2, b*J13*oh1^2, 0]);
            else
                a4 = min([b*J13*om1^2, b*J13*oh1^2]);
            end
            a5 = min([b*J23*om1*om2, b*J23*oh1*om2]);

            b = -b;
            a6 = min([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a7 = -b*J2*kd*oh2;
            a8 = min([-b*J23*kd*om3, -b*J23*kd*oh3]);

            a9 = min([- b*J12*kp*2*q0*q1, - b*J12*kp*2*qh0*q1, - b*J12*kp*2*q0*qh1, - b*J12*kp*2*qh0*qh1]);
            a10 = min([- b*J12*kp*2*q2*q3, - b*J12*kp*2*qh2*q3, - b*J12*kp*2*q2*qh3, - b*J12*kp*2*qh2*qh3]);
            a11 = min([- J2*b*kp*2*q0*q2, - J2*b*kp*2*qh0*q2, - J2*b*kp*2*q0*qh2, - J2*b*kp*2*qh0*qh2]);
            a12 = min([J2*b*kp*2*q1*q3, J2*b*kp*2*qh1*q3, J2*b*kp*2*q1*qh3, J2*b*kp*2*qh1*qh3]);

            ZZ1 = umax*sign(Ji2)*Ji2*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works

            b = -sign(Ji12);

            a1 = min([b*J12*om1*om3, b*J12*oh1*om3, b*J12*om1*oh3, b*J12*oh1*oh3]);
            a2 = min([(b*J2- b*J3)*om2*om3, (b*J2- b*J3)*om2*oh3]);
            if om3 <= 0  && 0 <= oh3
                a3 = min([b*J23*om3^2, b*J23*oh3^2, 0]);
            else
                a3 = min([b*J23*om3^2, b*J23*oh3^2]);
            end
            a4 = min([-b*J13*om1*om2, -b*J13*oh1*om2]);
            a5 = - b*J23*om2^2;

            b = -b;
            a6 = min([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a7 = -b*J2*kd*oh2;
            a8 = min([-b*J23*kd*om3, -b*J23*kd*oh3]);

            a9 = min([- b*J1*kp*2*q0*q1, - b*J1*kp*2*qh0*q1, - b*J1*kp*2*q0*qh1, - b*J1*kp*2*qh0*qh1]);
            a10 = min([- b*J1*kp*2*q2*q3, - b*J1*kp*2*qh2*q3, - b*J1*kp*2*q2*qh3, - b*J1*kp*2*qh2*qh3]);
            a11 = min([- J12*b*kp*2*q0*q2, - J12*b*kp*2*qh0*q2, - J12*b*kp*2*q0*qh2, - J12*b*kp*2*qh0*qh2]);
            a12 = min([J12*b*kp*2*q1*q3, J12*b*kp*2*qh1*q3, J12*b*kp*2*q1*qh3, J12*b*kp*2*qh1*qh3]);

            ZZ2 = umax*sign(Ji12)*Ji12*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works

            b = -sign(Ji23);

            a1 = min([b*(J1-J2)*om1*om2, b*(J1-J2)*oh1*om2]);
            a2 = b*J12*om2^2;
            a3 = min([b*J13*om2*om3, b*J13*om2*oh3]);
            if om1 <= 0  && 0 <= oh1
                a4 = min([-b*J12*om1^2, -b*J12*oh1^2, 0]);
            else
                a4 = min([-b*J12*om1^2, -b*J12*oh1^2]);
            end
            a5 = min([- b*J23*om1*om3, - b*J23*oh1*om3, - b*J23*om1*oh3, - b*J23*oh1*oh3]);

            b = -b;
            a6 = min([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a7 = -b*J2*kd*oh2;
            a8 = min([-b*J23*kd*om3, -b*J23*kd*oh3]);

            a9 = min([- b*J13*kp*2*q0*q1, - b*J13*kp*2*qh0*q1, - b*J13*kp*2*q0*qh1, - b*J13*kp*2*qh0*qh1]);
            a10 = min([- b*J13*kp*2*q2*q3, - b*J13*kp*2*qh2*q3, - b*J13*kp*2*q2*qh3, - b*J13*kp*2*qh2*qh3]);
            a11 = min([- J23*b*kp*2*q0*q2, - J23*b*kp*2*qh0*q2, - J23*b*kp*2*q0*qh2, - J23*b*kp*2*qh0*qh2]);
            a12 = min([J23*b*kp*2*q1*q3, J23*b*kp*2*qh1*q3, J23*b*kp*2*q1*qh3, J23*b*kp*2*qh1*qh3]);

            ZZ3 = umax*sign(Ji23)*Ji23*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works

            out = ZZ1 + ZZ2 + ZZ3;

        elseif prod(xh <= x)
             b = -sign(Ji2);

            a1 = max([b*(J3-J1)*om1*om3, b*(J3-J1)*oh1*om3, b*(J3-J1)*om1*oh3, b*(J3-J1)*oh1*oh3]);
            a2 = max([- b*J12*om2*om3, - b*J12*om2*oh3]);
            if oh3 <= 0  && 0 <= om3
                a3 = max([- b*J13*om3^2, - b*J13*oh3^2, 0]);
            else
                a3 = max([- b*J13*om3^2, - b*J13*oh3^2]);
            end
            if oh1 <= 0  && 0 <= om1
                a4 = max([b*J13*om1^2, b*J13*oh1^2, 0]);
            else
                a4 = max([b*J13*om1^2, b*J13*oh1^2]);
            end
            a5 = max([b*J23*om1*om2, b*J23*oh1*om2]);

            b = -b;
            a6 = max([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a7 = -b*J2*kd*oh2;
            a8 = max([-b*J23*kd*om3, -b*J23*kd*oh3]);

            a9 = max([- b*J12*kp*2*q0*q1, - b*J12*kp*2*qh0*q1, - b*J12*kp*2*q0*qh1, - b*J12*kp*2*qh0*qh1]);
            a10 = max([- b*J12*kp*2*q2*q3, - b*J12*kp*2*qh2*q3, - b*J12*kp*2*q2*qh3, - b*J12*kp*2*qh2*qh3]);
            a11 = max([- J2*b*kp*2*q0*q2, - J2*b*kp*2*qh0*q2, - J2*b*kp*2*q0*qh2, - J2*b*kp*2*qh0*qh2]);
            a12 = max([J2*b*kp*2*q1*q3, J2*b*kp*2*qh1*q3, J2*b*kp*2*q1*qh3, J2*b*kp*2*qh1*qh3]);

            ZZ1 = umax*sign(Ji2)*Ji2*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works


            b = -sign(Ji12);

            a1 = max([b*J12*om1*om3, b*J12*oh1*om3, b*J12*om1*oh3, b*J12*oh1*oh3]);
            a2 = max([(b*J2- b*J3)*om2*om3, (b*J2- b*J3)*om2*oh3]);
            if oh3 <= 0  && 0 <= om3
                a3 = max([b*J23*om3^2, b*J23*oh3^2, 0]);
            else
                a3 = max([b*J23*om3^2, b*J23*oh3^2]);
            end
            a4 = max([-b*J13*om1*om2, -b*J13*oh1*om2]);
            a5 = - b*J23*om2^2;

            b = -b;
            a6 = max([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a7 = -b*J2*kd*oh2;
            a8 = max([-b*J23*kd*om3, -b*J23*kd*oh3]);

            a9 = max([- b*J1*kp*2*q0*q1, - b*J1*kp*2*qh0*q1, - b*J1*kp*2*q0*qh1, - b*J1*kp*2*qh0*qh1]);
            a10 = max([- b*J1*kp*2*q2*q3, - b*J1*kp*2*qh2*q3, - b*J1*kp*2*q2*qh3, - b*J1*kp*2*qh2*qh3]);
            a11 = max([- J12*b*kp*2*q0*q2, - J12*b*kp*2*qh0*q2, - J12*b*kp*2*q0*qh2, - J12*b*kp*2*qh0*qh2]);
            a12 = max([J12*b*kp*2*q1*q3, J12*b*kp*2*qh1*q3, J12*b*kp*2*q1*qh3, J12*b*kp*2*qh1*qh3]);

            ZZ2 = umax*sign(Ji12)*Ji12*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works



            b = -sign(Ji23);

            a1 = max([b*(J1-J2)*om1*om2, b*(J1-J2)*oh1*om2]);
            a2 = b*J12*om2^2;
            a3 = max([b*J13*om2*om3, b*J13*om2*oh3]);
            if oh1 <= 0  && 0 <= om1
                a4 = max([-b*J12*om1^2, -b*J12*oh1^2, 0]);
            else
                a4 = max([-b*J12*om1^2, -b*J12*oh1^2]);
            end
            a5 = max([- b*J23*om1*om3, - b*J23*oh1*om3, - b*J23*om1*oh3, - b*J23*oh1*oh3]);


            b = -b;
            a6 = max([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a7 = -b*J2*kd*oh2;
            a8 = max([-b*J23*kd*om3, -b*J23*kd*oh3]);

            a9 = max([- b*J13*kp*2*q0*q1, - b*J13*kp*2*qh0*q1, - b*J13*kp*2*q0*qh1, - b*J13*kp*2*qh0*qh1]);
            a10 = max([- b*J13*kp*2*q2*q3, - b*J13*kp*2*qh2*q3, - b*J13*kp*2*q2*qh3, - b*J13*kp*2*qh2*qh3]);
            a11 = max([- J23*b*kp*2*q0*q2, - J23*b*kp*2*qh0*q2, - J23*b*kp*2*q0*qh2, - J23*b*kp*2*qh0*qh2]);
            a12 = max([J23*b*kp*2*q1*q3, J23*b*kp*2*qh1*q3, J23*b*kp*2*q1*qh3, J23*b*kp*2*qh1*qh3]);

            ZZ3 = umax*sign(Ji23)*Ji23*tanh((1/umax)*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)); % works

            out = ZZ1 + ZZ2 + ZZ3;

        else
            out = nan;
        end


    end

    function out = thing3(x, xh, J, Jinv)


        q0 = x(1);
        q1 = x(2);
        q2 = x(3);
        q3 = x(4);

        qh0 = xh(1);
        qh1 = xh(2);
        qh2 = xh(3);
        qh3 = xh(4);

        om1 = x(5);
        om2 = x(6);
        om3 = x(7);

        oh1 = xh(5);
        oh2 = xh(6);
        oh3 = xh(7);

        a = zeros(12, 1);
        if prod(x <= xh)
            b = -sign(Ji23);

  
            a(1) = min([(b*J3-b*J1)*om1*om3, (b*J3-b*J1)*oh1*om3]);
            a(2) = min([- b*J12*om2*om3, - b*J12*oh2*om3]);
            a(3) = - b*J13*om3^2;
            if om1 <= 0  && 0 <= oh1
                a(4) = min([b*J13*om1^2, b*J13*oh1^2, 0]);
            else
                a(4) = min([b*J13*om1^2, b*J13*oh1^2]);
            end
            a(5) = min([ b*J23*om1*om2,  b*J23*oh1*om2,  b*J23*om1*oh2,  b*J23*oh1*oh2]);

            b = -b;
        % -J1*kd*om1 - J12*kd*om2 - J13*kd*om3 - J1*kp*(2*q0*q1 + 2*q2*q3) - J12*kp*(2*q0*q2 - 2*q1*q3)
        % -J12*kd*om1 - J2*kd*om2 - J23*kd*om3 - J12*kp*(2*q0*q1 + 2*q2*q3) - J2*kp*(2*q0*q2 - 2*q1*q3)
        % -J13*kd*om1 - J23*kd*om2 - J3*kd*om3 - J13*kp*(2*q0*q1 + 2*q2*q3) - J23*kp*(2*q0*q2 - 2*q1*q3)

        % -J12*kd*om1 - J2*kd*om2 - J23*kd*om3 - J12*kp*(2*q0*q1 + 2*q2*q3) - J2*kp*(2*q0*q2 - 2*q1*q3)
            a(6) = min([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a(7) = min([-b*J2*kd*om2, -b*J2*kd*oh2]);
            a(8) = - b*J23*kd*om3;

            a(9) = min([- b*J12*kp*2*q0*q1, - b*J12*kp*2*qh0*q1, - b*J12*kp*2*q0*qh1, - b*J12*kp*2*qh0*qh1]);
            a(10) = min([- b*J12*kp*2*q2*q3, - b*J12*kp*2*qh2*q3, - b*J12*kp*2*q2*qh3, - b*J12*kp*2*qh2*qh3]);
            a(11) = min([- J2*b*kp*2*q0*q2, - J2*b*kp*2*qh0*q2, - J2*b*kp*2*q0*qh2, - J2*b*kp*2*qh0*qh2]);
            a(12) = min([J2*b*kp*2*q1*q3, J2*b*kp*2*qh1*q3, J2*b*kp*2*q1*qh3, J2*b*kp*2*qh1*qh3]);

            ZZ1 = umax*sign(Ji23)*Ji23*tanh((1/umax)*(sum(a))); % works






    % - Jin31*tanh(kd*om1 + om3*(J21*om1 + J22*om2 + J23*om3) - om2*(J31*om1 + J32*om2 + J33*om3) + kp*(2*q0*q1 + 2*q2*q3)) 

            b = -sign(Ji13);

    % b*kd*om1 + b*J12*om1*om3 + b*(J2 - J3)*om2*om3 + b*J23*om3^2 - b*J13*om1*om2 - b*J23*om2^2 + kp*(b*2*q0*q1 + b*2*q2*q3) 


            a(1) = min([b*J12*om1*om3, b*J12*oh1*om3]);

            a(2) = min([b*(J2 - J3)*om2*om3, b*(J2 - J3)*oh2*om3]);

            a(3) = b*J23*om3^2;

            a(4) = min([- b*J13*om1*om2, - b*J13*oh1*om2, - b*J13*om1*oh2, - b*J13*oh1*oh2]);

            if om2 <= 0  && 0 <= oh2
                a(5) = min([- b*J23*om2^2, - b*J23*oh2^2, 0]);
            else
                a(5) = min([- b*J23*om2^2, - b*J23*oh2^2]);
            end

            %ZZ2 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(a(1)+a(2)+a(3)+a(4)+a(5))); % works

        b = -b;
        % -J1*kd*om1 - J12*kd*om2 - J13*kd*om3 - J1*kp*(2*q0*q1 + 2*q2*q3) - J12*kp*(2*q0*q2 - 2*q1*q3)
            a(6) = min([-b*J1*kd*om1, -b*J1*kd*oh1]);
            a(7) = min([-b*J12*kd*om2, -b*J12*kd*oh2]);
            a(8) = - b*J13*kd*om3;

            a(9) = min([- b*J1*kp*2*q0*q1, - b*J1*kp*2*qh0*q1, - b*J1*kp*2*q0*qh1, - b*J1*kp*2*qh0*qh1]);
            a(10) = min([- b*J1*kp*2*q2*q3, - b*J1*kp*2*qh2*q3, - b*J1*kp*2*q2*qh3, - b*J1*kp*2*qh2*qh3]);
            a(11) = min([- J12*b*kp*2*q0*q2, - J12*b*kp*2*qh0*q2, - J12*b*kp*2*q0*qh2, - J12*b*kp*2*qh0*qh2]);
            a(12) = min([J12*b*kp*2*q1*q3, J12*b*kp*2*qh1*q3, J12*b*kp*2*q1*qh3, J12*b*kp*2*qh1*qh3]);

            ZZ2 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(sum(a))); % works









    % - Jin3*tanh(kd*om3 + om2*(J1*om1 + J12*om2 + J13*om3) - om1*(J12*om1 + J2*om2 + J23*om3))

            b = -sign(Ji3);

    % b*kd*om3 + (b*J1- b*J2)*om1*om2 + b*J12*om2^2 + b*J13*om2*om3 -b*J12*om1^2 - b*J23*om1*om3

            a(1) = min([(b*J1- b*J2)*om1*om2, (b*J1- b*J2)*oh1*om2, (b*J1- b*J2)*om1*oh2, (b*J1- b*J2)*oh1*oh2]);
            if om2 <= 0  && 0 <= oh2
                a(2) = min([b*J12*om2^2, b*J12*oh2^2, 0]);
            else
                a(2) = min([b*J12*om2^2, b*J12*oh2^2]);
            end
            a(3) = min([b*J13*om2*om3, b*J13*oh2*om3]);
            if om1 <= 0  && 0 <= oh1
                a(4) = min([-b*J12*om1^2, -b*J12*oh1^2, 0]);
            else
                a(4) = min([-b*J12*om1^2, -b*J12*oh1^2]);
            end
            a(5) = min([- b*J23*om1*om3, - b*J23*oh1*om3]);

            b = -b;
            % -J13*kd*om1 - J23*kd*om2 - J3*kd*om3 - J13*kp*(2*q0*q1 + 2*q2*q3) - J23*kp*(2*q0*q2 - 2*q1*q3)
            a(6) = min([-b*J13*kd*om1, -b*J13*kd*oh1]);
            a(7) = min([-b*J23*kd*om2, -b*J23*kd*oh2]);
            a(8) = - b*J3*kd*om3;

            a(9) = min([- b*J13*kp*2*q0*q1, - b*J13*kp*2*qh0*q1, - b*J13*kp*2*q0*qh1, - b*J13*kp*2*qh0*qh1]);
            a(10) = min([- b*J13*kp*2*q2*q3, - b*J13*kp*2*qh2*q3, - b*J13*kp*2*q2*qh3, - b*J13*kp*2*qh2*qh3]);
            a(11) = min([- J23*b*kp*2*q0*q2, - J23*b*kp*2*qh0*q2, - J23*b*kp*2*q0*qh2, - J23*b*kp*2*qh0*qh2]);
            a(12) = min([J23*b*kp*2*q1*q3, J23*b*kp*2*qh1*q3, J23*b*kp*2*q1*qh3, J23*b*kp*2*qh1*qh3]);

            %ZZ2 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+a(7)+a(8)+a(9)+a(10)+a(11)+a(12))); % works

            ZZ3 = umax*sign(Ji3)*Ji3*tanh((1/umax)*(sum(a))); % works

            out = ZZ1 + ZZ2 + ZZ3;


        elseif prod(xh <= x)
            b = -sign(Ji23);

    % b*kd*om2 (b*J3-b*J1)*om1*om3 - b*J12*om2*om3 - b*J13*om3^2 + b*J13*om1^2 + b*J23*om1*om2 + kp*(b*2*q0*q2 - b*2*q1*q3))

            a(1) = max([(b*J3-b*J1)*om1*om3, (b*J3-b*J1)*oh1*om3]);
            a(2) = max([- b*J12*om2*om3, - b*J12*oh2*om3]);
            a(3) = - b*J13*om3^2;
            if oh1 <= 0  && 0 <= om1
                a(4) = max([b*J13*om1^2, b*J13*oh1^2, 0]);
            else
                a(4) = max([b*J13*om1^2, b*J13*oh1^2]);
            end
            a(5) = max([ b*J23*om1*om2,  b*J23*oh1*om2,  b*J23*om1*oh2,  b*J23*oh1*oh2]);

            b = -b;
        % -J1*kd*om1 - J12*kd*om2 - J13*kd*om3 - J1*kp*(2*q0*q1 + 2*q2*q3) - J12*kp*(2*q0*q2 - 2*q1*q3)
        % -J12*kd*om1 - J2*kd*om2 - J23*kd*om3 - J12*kp*(2*q0*q1 + 2*q2*q3) - J2*kp*(2*q0*q2 - 2*q1*q3)
        % -J13*kd*om1 - J23*kd*om2 - J3*kd*om3 - J13*kp*(2*q0*q1 + 2*q2*q3) - J23*kp*(2*q0*q2 - 2*q1*q3)

        % -J12*kd*om1 - J2*kd*om2 - J23*kd*om3 - J12*kp*(2*q0*q1 + 2*q2*q3) - J2*kp*(2*q0*q2 - 2*q1*q3)
            a(6) = max([-b*J12*kd*om1, -b*J12*kd*oh1]);
            a(7) = max([-b*J2*kd*om2, -b*J2*kd*oh2]);
            a(8) = - b*J23*kd*om3;

            a(9) = max([- b*J12*kp*2*q0*q1, - b*J12*kp*2*qh0*q1, - b*J12*kp*2*q0*qh1, - b*J12*kp*2*qh0*qh1]);
            a(10) = max([- b*J12*kp*2*q2*q3, - b*J12*kp*2*qh2*q3, - b*J12*kp*2*q2*qh3, - b*J12*kp*2*qh2*qh3]);
            a(11) = max([- J2*b*kp*2*q0*q2, - J2*b*kp*2*qh0*q2, - J2*b*kp*2*q0*qh2, - J2*b*kp*2*qh0*qh2]);
            a(12) = max([J2*b*kp*2*q1*q3, J2*b*kp*2*qh1*q3, J2*b*kp*2*q1*qh3, J2*b*kp*2*qh1*qh3]);

            ZZ1 = umax*sign(Ji23)*Ji23*tanh((1/umax)*(sum(a))); % works






    % - Jin31*tanh(kd*om1 + om3*(J21*om1 + J22*om2 + J23*om3) - om2*(J31*om1 + J32*om2 + J33*om3) + kp*(2*q0*q1 + 2*q2*q3)) 

            b = -sign(Ji13);

    % b*kd*om1 + b*J12*om1*om3 + b*(J2 - J3)*om2*om3 + b*J23*om3^2 - b*J13*om1*om2 - b*J23*om2^2 + kp*(b*2*q0*q1 + b*2*q2*q3) 


            a(1) = max([b*J12*om1*om3, b*J12*oh1*om3]);

            a(2) = max([b*(J2 - J3)*om2*om3, b*(J2 - J3)*oh2*om3]);

            a(3) = b*J23*om3^2;

            a(4) = max([- b*J13*om1*om2, - b*J13*oh1*om2, - b*J13*om1*oh2, - b*J13*oh1*oh2]);

            if oh2 <= 0  && 0 <= om2
                a(5) = max([- b*J23*om2^2, - b*J23*oh2^2, 0]);
            else
                a(5) = max([- b*J23*om2^2, - b*J23*oh2^2]);
            end

            %ZZ2 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(a(1)+a(2)+a(3)+a(4)+a(5))); % works

        b = -b;
        % -J1*kd*om1 - J12*kd*om2 - J13*kd*om3 - J1*kp*(2*q0*q1 + 2*q2*q3) - J12*kp*(2*q0*q2 - 2*q1*q3)
            a(6) = max([-b*J1*kd*om1, -b*J1*kd*oh1]);
            a(7) = max([-b*J12*kd*om2, -b*J12*kd*oh2]);
            a(8) = - b*J13*kd*om3;

            a(9) = max([- b*J1*kp*2*q0*q1, - b*J1*kp*2*qh0*q1, - b*J1*kp*2*q0*qh1, - b*J1*kp*2*qh0*qh1]);
            a(10) = max([- b*J1*kp*2*q2*q3, - b*J1*kp*2*qh2*q3, - b*J1*kp*2*q2*qh3, - b*J1*kp*2*qh2*qh3]);
            a(11) = max([- J12*b*kp*2*q0*q2, - J12*b*kp*2*qh0*q2, - J12*b*kp*2*q0*qh2, - J12*b*kp*2*qh0*qh2]);
            a(12) = max([J12*b*kp*2*q1*q3, J12*b*kp*2*qh1*q3, J12*b*kp*2*q1*qh3, J12*b*kp*2*qh1*qh3]);

            ZZ2 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(sum(a))); % works









    % - Jin3*tanh(kd*om3 + om2*(J1*om1 + J12*om2 + J13*om3) - om1*(J12*om1 + J2*om2 + J23*om3))

            b = -sign(Ji3);

    % b*kd*om3 + (b*J1- b*J2)*om1*om2 + b*J12*om2^2 + b*J13*om2*om3 -b*J12*om1^2 - b*J23*om1*om3

            a(1) = max([(b*J1- b*J2)*om1*om2, (b*J1- b*J2)*oh1*om2, (b*J1- b*J2)*om1*oh2, (b*J1- b*J2)*oh1*oh2]);
            if oh2 <= 0  && 0 <= om2
                a(2) = max([b*J12*om2^2, b*J12*oh2^2, 0]);
            else
                a(2) = max([b*J12*om2^2, b*J12*oh2^2]);
            end
            a(3) = max([b*J13*om2*om3, b*J13*oh2*om3]);
            if oh1 <= 0  && 0 <= om1
                a(4) = max([-b*J12*om1^2, -b*J12*oh1^2, 0]);
            else
                a(4) = max([-b*J12*om1^2, -b*J12*oh1^2]);
            end
            a(5) = max([- b*J23*om1*om3, - b*J23*oh1*om3]);

            b = -b;
            % -J13*kd*om1 - J23*kd*om2 - J3*kd*om3 - J13*kp*(2*q0*q1 + 2*q2*q3) - J23*kp*(2*q0*q2 - 2*q1*q3)
            a(6) = max([-b*J13*kd*om1, -b*J13*kd*oh1]);
            a(7) = max([-b*J23*kd*om2, -b*J23*kd*oh2]);
            a(8) = - b*J3*kd*om3;

            a(9) = max([- b*J13*kp*2*q0*q1, - b*J13*kp*2*qh0*q1, - b*J13*kp*2*q0*qh1, - b*J13*kp*2*qh0*qh1]);
            a(10) = max([- b*J13*kp*2*q2*q3, - b*J13*kp*2*qh2*q3, - b*J13*kp*2*q2*qh3, - b*J13*kp*2*qh2*qh3]);
            a(11) = max([- J23*b*kp*2*q0*q2, - J23*b*kp*2*qh0*q2, - J23*b*kp*2*q0*qh2, - J23*b*kp*2*qh0*qh2]);
            a(12) = max([J23*b*kp*2*q1*q3, J23*b*kp*2*qh1*q3, J23*b*kp*2*q1*qh3, J23*b*kp*2*qh1*qh3]);

            %ZZ2 = umax*sign(Ji13)*Ji13*tanh((1/umax)*(a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+a(7)+a(8)+a(9)+a(10)+a(11)+a(12))); % works

            ZZ3 = umax*sign(Ji3)*Ji3*tanh((1/umax)*(sum(a))); % works

            out = ZZ1 + ZZ2 + ZZ3;

        else
            out = nan;
        end
    end

end