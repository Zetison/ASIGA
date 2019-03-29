function parms = physical2paramMapping(v, nurbs, guess)

if strcmp(nurbs.type, '3Dvolume')

    if nargin == 2
        xi = 0.5;
        eta = 0.5;
        zeta = 0.5;
    else 
        xi = guess(1);
        eta = guess(2);
        zeta = guess(3);
    end


    f = @(parms) norm(evaluateNURBS(nurbs, parms) - v');

    % options = optimset('Display','iter');

    parms = fminsearchbnd(f,[xi, eta, zeta],[0, 0, 0], [1, 1, 1]);
elseif strcmp(nurbs.type, '3Dsurface')

    if nargin == 2
        xi = 0.5;
        eta = 0.5;
    else 
        xi = guess(1);
        eta = guess(2);
    end


    f = @(eta) norm(evaluateNURBS(nurbs, [xi, eta]) - v');

    options = optimset('TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',10000, 'MaxIter', 10000);
%     options = optimset('Display','iter', 'TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',10000, 'MaxIter', 10000);

    parms = fminsearchbnd(f,eta,0,1,options);
end