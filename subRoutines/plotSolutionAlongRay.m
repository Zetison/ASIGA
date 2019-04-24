if true % plot only inside artificial boundary
%     xi = 0;
%     eta = 0.5;
    npts = 10000;
    zeta_arr = linspace(0,1,npts);

    Xi = varCol.nurbs.knots{1};
    Eta = varCol.nurbs.knots{2};
    Zeta = varCol.nurbs.knots{3};
    p_xi = varCol.nurbs.degree(1);
    p_eta = varCol.nurbs.degree(2);
    p_zeta = varCol.nurbs.degree(3);

    weights = varCol.weights;
    controlPts = varCol.controlPts;
    noDofs = varCol.noDofs;
    nurbs = varCol.nurbs;
    d = varCol.dimension;
    runInParallell = varCol.runInParallell;
    omega = varCol.omega;
    v = zeros(npts,3);
    p_h_arr = zeros(npts,1);
    gp_h_arr = zeros(npts,3);
    alpha_f = alpha_f_arr;
    beta_f = beta_f_arr;
    d_vec2 = [cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)].';
    parm_pt = pointInversion(nurbs,R_o*d_vec2,10*eps);
    xi = parm_pt(1);
    eta = parm_pt(2);
%     for i = 1:length(zeta_arr)
    parfor (i = 1:length(zeta_arr), runInParallell)
        zeta = zeta_arr(i);
        [p, X, dpdx, dpdy, dpdz] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U_fluid_o);
        v(i,:) = X;
        p_h_arr(i) = p;
        gp_h_arr(i,:) = [dpdx, dpdy, dpdz];
    end

    data = e3Dss(v, e3Dss_options);
    d_vec = k_vec/k(1);
    p_inc = @(v) P_inc*exp(1i*k*dot3(v,d_vec));
    gp_inc = @(v) -1i*p_inc(v)*k(1)*dot(d_vec,d_vec2);

    p_tot = data(1).p+p_inc(v);
    p_h_tot = p_h_arr+p_inc(v);
    vel = -[data(1).dpdx, data(1).dpdy, data(1).dpdz]*d_vec2+gp_inc(v);
    vel_h = -gp_h_arr*d_vec2+gp_inc(v);
    if true
        plot(norm2(v), abs(p_tot), '--', 'color','blue')
        hold on
        plot(norm2(v), abs(p_h_tot), 'color','blue')
        plot(norm2(v), abs(vel), '--', 'color','green')
        plot(norm2(v), abs(vel_h), 'color','green')
        legend('|p_{tot}|', '|p_{tot,h}|', '|v|', '|v_h|')
    else
        plot(norm2(v), real(p_tot), '--', 'color','blue')
        hold on
        plot(norm2(v), real(p_h_tot), 'color','blue')
        plot(norm2(v), imag(p_tot), '--', 'color','green')
        plot(norm2(v), imag(p_h_tot), 'color','green')

        plot(norm2(v), real(vel), '--', 'color','red')
        plot(norm2(v), real(vel_h), 'color','red')
        plot(norm2(v), imag(vel), '--', 'color','magenta')
        plot(norm2(v), imag(vel_h), 'color','magenta')
        legend('real(p)', 'real(p_h)', 'imag(p)', 'imag(p_h)', 'real(v)', 'real(v_h)', 'imag(v)', 'imag(v_h)')
    end
    

    subFolderName = [folderName '/' saveName '/Error_plots'];
    if ~exist(subFolderName, 'dir')
        mkdir(subFolderName);
    end
    savefig([subFolderName '/' saveName '_' coreMethod method '_mesh' num2str(M) '_degree' num2str(degreeElev+2) ...
                            '_formulation' formulation])
    figure(42)
    error_p = abs(p_tot-p_h_tot)/max(abs(p_h_tot));
    error_vel = abs(vel-vel_h)/max(abs(vel));
    semilogy(norm2(v), error_p, norm2(v), error_vel)
    legend('relative error in p_h', 'relative error in v_h')
    hold on
    uniqueZeta = unique(Zeta);
    for i = 1:length(uniqueZeta)
        zeta = uniqueZeta(i);
        v = evaluateNURBS(nurbs,[xi,eta,zeta]);
        r = norm(v);
        plot(r*[1, 1], ylim, '--', 'color','black')
        if i > 1
            plot((r+r_old)/2*[1, 1], ylim, ':', 'color','black')
        end
        r_old = r;            
    end
    fullSaveName = [subFolderName '/' saveName '_' coreMethod method '_mesh' num2str(M) '_degree' num2str(degreeElev+2) ...
                            '_formulation' formulation '_Error'];
    savefig(fullSaveName)

    printResultsToFile(fullSaveName, h_max(:,f_Nr), Errors(:,f_Nr), varCol, 0, 1); 
else
    v = P_far;
    plot(norm2(v),real(p_h))
    data = e3Dss(v, e3Dss_options);
    d_vec = k_vec/k(1);
    p_inc = @(v) P_inc*exp(1i*k*dot3(v,d_vec));
    gp_inc = @(v) -1i*p_inc(v)*k(1);

    plot(norm2(v), real(p_h+p_inc(v)), norm2(v), real(data(1).p+p_inc(v)), norm2(v), real(-[data(1).dpdx, data(1).dpdy, data(1).dpdz]*d_vec+gp_inc(v)))
    legend('p_h', 'p', 'v')
end