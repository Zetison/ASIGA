
x_0 = [0, 0, 0]; % The origin of the model
switch method
    case {'IE','IENSG'}
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        A_2 = [1 0 0;
              0 1 0;
              0 0 1];
        varCol.x_0 = x_0;
        varCol.A_2 = A_2;
end

if varCol.boundaryMethod
    principalLengthXiDir = R_o;
    principalLengthEtaDir = R_o;

    L_gamma = R_o;
    eta2 = (R_o+L)/(L+2*R_o);
    eta1 = R_o/(L+2*R_o);
    switch model
        case 'M5A'
            fluid = getModel5Data(R_o, eta1, eta2, L, l, 'Xaxis');
        case 'M5B'
            fluid = getModel5Data(R_o, eta1, eta2, L, l, 'Zaxis');
    end

    fluid = elevateNURBSdegree(fluid,[1 1]*degreeElev);
    if M ~= 0
        fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, 1) [linspace2(eta1/2,eta2/2,round(L/R_o)) linspace2((eta1+1)/2,(eta2+1)/2,round(L/R_o))]});
    end
    noNewXiKnots = 2^M-1; % 8*i_mesh
    noNewEtaKnots = noNewXiKnots;
    fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, noNewXiKnots) ...
                                      insertUniform2(fluid.knots{2}, noNewEtaKnots)});
end