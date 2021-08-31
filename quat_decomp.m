

function out = quat_decomp(x, w, xhat, what)

out = [d1(x,w,xhat,what); ...
       d2(x,w,xhat,what); ...
       d3(x,w,xhat,what); ...
       d4(x,w,xhat,what)];

    function out = d1(x,w,xhat,what)
        if x(end)<=xhat(end)
            t1 = psi(x(2),x(5),xhat(2),xhat(5));
            t2 = psi(x(3),x(6),xhat(3),xhat(6));
            t3 = psi(x(4),x(7),xhat(4),xhat(7));
            out = 0.5*(-t1-t2-t3);
        else
            t1 = phi(x(2),x(5),xhat(2),xhat(5));
            t2 = phi(x(3),x(6),xhat(3),xhat(6));
            t3 = phi(x(4),x(7),xhat(4),xhat(7));
            out = 0.5*(-t1-t2-t3);
        end
    end

    function out = d2(x,w,xhat,what)
        if x(end)<=xhat(end)
            t1 =  phi(x(1),x(5),xhat(1),xhat(5));
            t2 = -psi(x(4),x(6),xhat(4),xhat(6));
            t3 =  phi(x(3),x(7),xhat(3),xhat(7));
            out = 0.5*(t1+t2+t3);
        else
            t1 =  psi(x(1),x(5),xhat(1),xhat(5));
            t2 = -phi(x(4),x(6),xhat(4),xhat(6));
            t3 =  psi(x(3),x(7),xhat(3),xhat(7));
            out = 0.5*(t1+t2+t3);
        end
    end

    function out = d3(x,w,xhat,what)
        if x(end)<=xhat(end)
            t1 = phi(x(4),x(5),xhat(4),xhat(5));
            t2 = phi(x(1),x(6),xhat(1),xhat(6));
            t3 = -psi(x(2),x(7),xhat(2),xhat(7));
            out = 0.5*(t1+t2+t3);
        else
            t1 = psi(x(4),x(5),xhat(4),xhat(5));
            t2 = psi(x(1),x(6),xhat(1),xhat(6));
            t3 = -phi(x(2),x(7),xhat(2),xhat(7));
            out = 0.5*(t1+t2+t3);
        end
    end

    function out = d4(x,w,xhat,what)
        if x(end)<=xhat(end)
            t1 = -psi(x(3),x(5),xhat(3),xhat(5));
            t2 =  phi(x(2),x(6),xhat(2),xhat(6));
            t3 =  phi(x(1),x(7),xhat(1),xhat(7));
            out = 0.5*(t1+t2+t3);
        else
            t1 = -phi(x(3),x(5),xhat(3),xhat(5));
            t2 =  psi(x(2),x(6),xhat(2),xhat(6));
            t3 =  psi(x(1),x(7),xhat(1),xhat(7));
            out = 0.5*(t1+t2+t3);
        end
    end

    function out = phi(xi, xj, xhati, xhatj)
        out = min([xi*xj, xhati*xj, xi*xhatj, xhati*xhatj]);
    end

    function out = psi(xi, xj, xhati, xhatj)
        out = max([xi*xj, xhati*xj, xi*xhatj, xhati*xhatj]);
    end

end
