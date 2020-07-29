function [xi, eta] = prolateCoordsToParamCoords(nurbs, x,y,z, Upsilon, chi, A_2, x_0)


switch nurbs.type
    case '3Dvolume'
        X = [x,y,z];
        Xt = A_2*(X-x_0).';
        xt = Xt(1);
        yt = Xt(2);
        zt = Xt(3);
        [~, theta, phi] = evaluateProlateCoords(xt,yt,zt,Upsilon);

        xt = sqrt(chi^2-Upsilon^2)*sin(theta)*cos(phi);
        yt = sqrt(chi^2-Upsilon^2)*sin(theta)*sin(phi);
        zt = chi*cos(theta);
        Xt = [xt,yt,zt];
        X = A_2\Xt.' + x_0.';
        
        f = @(parms) norm(evaluateNURBS(nurbs, [parms, 1]) - X);

        options = optimset('fminsearch');
        options.TolX = 1.e-13;
        options.TolFun = 1.e-13;
%         options.Display = 'iter';
%         options.TolX = 1.e-12;
        [parms, fval] = fminsearchbnd(f,[0, 0],[0, 0], [1, 1], options);
        if fval > 1e-10
            parms = fminsearchbnd(f,[1, 0],[0, 0], [1, 1], options);
        end
        if fval > 1e-10
            parms = fminsearchbnd(f,[0.5, 0],[0, 0], [1, 1], options);
        end
        if fval > 1e-10
            parms = fminsearchbnd(f,[0.25, 0],[0, 0], [1, 1], options);
        end
        if fval > 1e-10
            parms = fminsearchbnd(f,[0.75, 0],[0, 0], [1, 1], options);
        end
        xi = parms(1);
        eta = parms(2);
    case '2Dsurface'
        [~, theta, phi] = evaluateProlateCoords(x,y,z,Upsilon);

        x = sqrt(chi^2-Upsilon^2)*sin(theta)*cos(phi);
        y = sqrt(chi^2-Upsilon^2)*sin(theta)*sin(phi);
        z = chi*cos(theta);
        if x > 0
            xi = 0;
            parms = physical2paramMapping([x,y,z], nurbs, [xi, 0.5]);
        else
            xi = 0.5;
            parms = physical2paramMapping([x,y,z], nurbs, [xi, 0.5]);
        end

        eta = parms(1);
end