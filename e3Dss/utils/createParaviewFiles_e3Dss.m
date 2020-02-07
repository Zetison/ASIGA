function createParaviewFiles_e3Dss(extraPts, vtfFileName, options)

ESBC = options.ESBC;
SHBC = options.SHBC;
SSBC = options.SSBC;
R_i = options.R_i;
R_o = options.R_o;
R_a = options.R_a;
% P_inc = options.P_inc;
d_vec = options.d_vec;

M = length(R_o);
plotTimeOscillation = options.plotTimeOscillation;
plotInTimeDomain = options.plotInTimeDomain;
if ~plotInTimeDomain
    omega = options.omega;
end

rho_f = options.rho_f;
c_f = options.c_f;
noDomains = 2*M+1-ESBC-2*SHBC-SSBC;
nodes = cell(noDomains,1);
visElements = cell(noDomains,1);
if isfield(options, 'computeForSolidDomain')
    computeForSolidDomain = options.computeForSolidDomain;
else
    computeForSolidDomain = 0;
end
if computeForSolidDomain
     options.calc_sigma_xx = 1;
     options.calc_sigma_yy = 1; 
     options.calc_sigma_zz = 1;
     options.calc_sigma_yz = 1;
     options.calc_sigma_xz = 1;
     options.calc_sigma_xy = 1; 
     
     options.calc_u_x = 1; 
     options.calc_u_y = 1; 
     options.calc_u_z = 1; 
     
     options.calc_du_xdx = 1; 
     options.calc_du_xdy = 1; 
     options.calc_du_xdz = 1; 
     options.calc_du_ydx = 1; 
     options.calc_du_ydy = 1; 
     options.calc_du_ydz = 1; 
     options.calc_du_zdx = 1; 
     options.calc_du_zdy = 1; 
     options.calc_du_zdz = 1; 
end
% meshRectangle
tic
for m = 1:M+1-ESBC-SHBC-SSBC
    if m ~= M+1 && ~(SHBC && m == M) && computeForSolidDomain
        if ESBC && m == M
%             solid = getESBCData(R_o(m),'Xaxis');
%             t_solid = R_o(m);
%             R_solid = R_o(m);
            [visElements{2*m}, nodes{2*m}] = mesh2DDisk(R_o(m),round(1.18*extraPts*2*R_o(m)));
            nodes{2*m} = [nodes{2*m}, zeros(size(nodes{2*m},1),1)];
        else
            [x,y,z] = sphere(round(R_o(m)*pi*extraPts));
            surfObj = surf2patch(x,y,z,'triangles');
            visElements{2*m} = [surfObj.faces, surfObj.faces+size(surfObj.vertices,1)];
            nodes{2*m} = [R_i(m)*surfObj.vertices; R_o(m)*surfObj.vertices];
        end
    end
    
    if m == 1
        [visElements{2*m-1}, nodes{2*m-1}] = meshRectangleWcircHole([-1.3*R_a, -R_a],[1.3*R_a, R_a],R_o(1),round(extraPts*1.3*R_a));
        nodes{2*m-1} = [nodes{2*m-1}, zeros(size(nodes{2*m-1},1),1)];
    elseif m < M+1
        [visElements{2*m-1}, nodes{2*m-1}] = mesh2DDonut(R_i(m-1),R_o(m),round(1.18*extraPts*2*R_i(m-1)));
        nodes{2*m-1} = [nodes{2*m-1}, zeros(size(nodes{2*m-1},1),1)];
    else
        [visElements{2*m-1}, nodes{2*m-1}] = mesh2DDisk(R_i(m-1),round(1.18*extraPts*2*R_i(m-1)));
        nodes{2*m-1} = [nodes{2*m-1}, zeros(size(nodes{2*m-1},1),1)];
    end
end
if plotInTimeDomain
    [visTemp, temp_nodes] = meshRectangle([-1.3*R_a, 0],[1.3*R_a, 0.5*R_a],round(extraPts*1.3*R_a));
    nodes{1} = [nodes{1}; temp_nodes(:,1), R_a*ones(size(temp_nodes,1),1), temp_nodes(:,2)];
    visElements{1} = [visElements{1}; visTemp+max(max(visElements{1}))];
end
% keyboard
% nodes(2:2:end) = zeros(1,0);
npts = 0;
for i = 1:length(nodes)
    npts = npts + size(nodes{i},1);
end
npts
% keyboard
% nodes = nodes(1);
disp(['Time spent building mesh ' num2str(toc) ' seconds.'])
% keyboard
% nodes = nodes(1:end-1);
% nodes = nodes(1);
if plotTimeOscillation
    N = 30;
    N_fine = 30;
    type = 1;
    startIdx = 1;
    options.P_inc = 1;
elseif plotInTimeDomain
    f_c = options.f_c;
    N = options.N; % 200
    N_fine = 2*N;
    switch options.applyLoad
        case 'planeWave'
            startIdx = 1900; % 900
        case 'pointCharge'
            startIdx = 2000;
        case 'radialPulsation'
%             N = 2*N;
            startIdx = 2*1900;
    end
    T = options.T; %N/M;
    B = N/T; % bandwidth
    
    f_R = B/2;
    df = 1/T;
    f = linspace(0,f_R-df,N/2);
    omega = 2*pi*f;
    options.omega = omega(2:end);
    type = 1;
    totnpts = npts*N*4
    omega_c = 2*pi*f_c;
    options.P_inc = @(omega) P_inc_(omega,omega_c,type);
%     keyboard
else
    N = 1;
    N_fine = 1;
    type = 1;
    startIdx = 1;
    options.P_inc = 1;
end
tic
data = e3Dss(nodes, options);
disp(['Time spent on computing solution: ' num2str(toc) ' seconds.'])
m = 1;
[pathstr,filename] = fileparts(vtfFileName);

for j = 1:length(nodes)
    if ~(SSBC && j == length(nodes))
        tic
        clear VTKdata
        if mod(j,2) == 0 && computeForSolidDomain
            toggleJacobianMatrix = false(1,9);
            if ESBC && m == M
                VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', plotTimeOscillation, ...
                            'plotSphericalRadialDisplacement',1, 'plotdu_xdx',1, 'plotdu_xdy',1, 'plotdu_xdz',1, 'plotdu_ydx',1, 'plotdu_ydy',1, ...
                            'plotdu_ydz',1, 'plotdu_zdx',1, 'plotdu_zdy',1, 'plotdu_zdz',1, ...
                            'plotDisplacementVectors',0,'plotSphericalStress_rr',1,'plotStressXX',1,'plotStressZZ',1,'plotVonMisesStress',1); 
                toggleJacobianMatrix = [VTKoptions.plotdu_xdx VTKoptions.plotdu_xdy VTKoptions.plotdu_xdz ...
                                        VTKoptions.plotdu_ydx VTKoptions.plotdu_ydy VTKoptions.plotdu_ydz ...
                                        VTKoptions.plotdu_zdx VTKoptions.plotdu_zdy VTKoptions.plotdu_zdz];
            else
                VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_WEDGE', 'plotTimeOscillation', plotTimeOscillation, ...
                            'plotSphericalRadialDisplacement',0, 'plotDisplacementVectors',0,'plotSphericalStress_rr',1,'plotVonMisesStress',0); 
            end
                    
            if VTKoptions.plotDisplacementVectors || VTKoptions.plotSphericalRadialDisplacement
                displacement = zeros(size(nodes{j},1),3,length(omega));
                for i = 2:length(omega)
                    displacement(:,:,i) = [data(m).u_x(:,i-1) data(m).u_y(:,i-1) data(m).u_z(:,i-1)];
                end
                if ~plotInTimeDomain
                    VTKdata.displacement = [data(m).u_x data(m).u_y data(m).u_z];
                    
                    temp = VTKdata.displacement;
                    VTKdata.displacement = zeros([size(VTKdata.displacement), N]);
                    for i = 1:N
                        t = (i-1)/N*2*pi/omega;
                        VTKdata.displacement(:,:,i) = real(temp*exp(-1i*omega*t));
                    end
                end
            end
            if any(toggleJacobianMatrix)
                VTKdata.jacobian = zeros(size(nodes{j},1),6,length(omega));
                
                if ~plotInTimeDomain
                    VTKdata.jacobian = [data(m).du_xdx data(m).du_xdy data(m).du_xdz data(m).du_ydx data(m).du_ydy data(m).du_ydz data(m).du_zdx data(m).du_zdy data(m).du_zdz];
                    
                    temp = VTKdata.jacobian;
                    VTKdata.jacobian = zeros([size(VTKdata.jacobian), N]);
                    for i = 1:N
                        t = (i-1)/N*2*pi/omega;
                        VTKdata.jacobian(:,:,i) = real(temp*exp(-1i*omega*t));
                    end
                end
            end
            if VTKoptions.plotVonMisesStress || VTKoptions.plotSphericalStress_rr 
                VTKdata.stress = zeros(size(nodes{j},1),6,length(omega));
                
                if plotInTimeDomain
                    for i = 2:length(omega)
                        VTKdata.stress(:,:,i) = [data(m).sigma_xx(:,i-1) data(m).sigma_yy(:,i-1) data(m).sigma_zz(:,i-1) data(m).sigma_yz(:,i-1) data(m).sigma_xz(:,i-1) data(m).sigma_xy(:,i-1)];
                    end
                    VTKdata.stress = 2/T*real(fft(VTKdata.stress,N_fine,3));
                    temp = VTKdata.stress;
                    VTKdata.stress(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
                    VTKdata.stress(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
                else
                    VTKdata.stress = [data(m).sigma_xx data(m).sigma_yy data(m).sigma_zz data(m).sigma_yz data(m).sigma_xz data(m).sigma_xy];
                    
                    temp = VTKdata.stress;
                    VTKdata.stress = zeros([size(VTKdata.stress), N]);
                    for i = 1:N
                        t = (i-1)/N*2*pi/omega;
                        VTKdata.stress(:,:,i) = real(temp*exp(-1i*omega*t));
                    end
                end
            end
        elseif mod(j,2) ~= 0
            VTKoptions = struct('name',[pathstr '/fluid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', ...
                             'plotTimeOscillation', plotTimeOscillation, 'plotDisplacementVectors', 0, ...
                             'plotSphericalRadialDisplacement',0, 'plotTotField', 1, 'plotScalarField', 1, 'plotSPL', 1);
            if ~(plotInTimeDomain || plotTimeOscillation)
                VTKoptions.plotTotFieldAbs = 1;
            end
            if VTKoptions.plotDisplacementVectors 
                displacement = zeros(size(nodes{j},1),3,length(omega));
            end
            VTKdata.totField = zeros(size(nodes{j},1),1,length(omega));
            if m == 1
                if plotInTimeDomain
                    VTKdata.P_inc = zeros(size(nodes{j},1),1,length(omega));
                    for i = 2:length(omega)                     
                        switch options.applyLoad
                            case 'pointCharge'
                                r_s = options.r_s;
                                k = omega(i)/c_f(1);
                                r = @(y) norm2(repmat(-r_s*d_vec.',size(y,1),1)-y);
                                Phi_k = @(y) exp(1i*k*r(y))./(4*pi*r(y));     
                                p_inc = @(y) 4*pi*r_s*P_inc_(omega(i),omega_c,type)*Phi_k(y);
                            case 'planeWave'
                                k = omega(i)/c_f(1);
                                k_vec = options.d_vec*k;
                                p_inc = @(v) P_inc_(omega(i),omega_c,type).*exp(1i*dot3(v, k_vec));
                            case 'radialPulsation'
                                k = omega(i)/c_f(1);
                                p_inc = @(v) P_inc_(omega(i),omega_c,type)*R_o(1).*exp(-1i*k*(norm2(v)-R_o(1)))./norm2(v);
                        end
                        VTKdata.P_inc(:,:,i) = p_inc(nodes{1});
                        VTKdata.totField(:,:,i) = data(m).p(:,i-1) + VTKdata.P_inc(:,:,i);
                    end
                else
                    k = omega/c_f(1);
                    k_vec = options.d_vec*k;
                    p_inc = @(v) exp(1i*dot3(v, k_vec));   
                    VTKdata.P_inc = p_inc(nodes{1});
                    VTKdata.totField = data(m).p(:,1) + VTKdata.P_inc;
                    VTKdata.scalarField = data(m).p(:,1);
                end
                VTKdata = rmfield(VTKdata,'P_inc');
            else
                VTKoptions.plotScalarField = false;
                VTKoptions.plotP_inc = false;
                if plotInTimeDomain
                    for i = 2:length(omega)
    %                     gScalarField = [data(m).dpdx(:,i) data(m).dpdy(:,i) data(m).dpdz(:,i)];
    %                     displacement(:,:,i) = gScalarField/(rho_f(m)*omega(i)^2);
                        VTKdata.totField(:,:,i) = data(m).p(:,i-1);
                    end
                else
                    VTKdata.totField = data(m).p;
                end
            end
            if plotInTimeDomain
                VTKoptions.plotScalarField = false;
                VTKoptions.plotSPL = false;
                VTKdata.totField = 2/T*real(fft(VTKdata.totField,N_fine,3));
                temp = VTKdata.totField;
                VTKdata.totField(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
                VTKdata.totField(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
            elseif plotTimeOscillation
                temp = VTKdata.totField;
                VTKdata.totField = zeros([size(VTKdata.totField), N]);
                for i = 1:N
                    t = (i-1)/N*2*pi/omega;
                    VTKdata.totField(:,:,i) = real(temp*exp(-1i*omega*t));
                end
            else
                VTKdata.totFieldAbs = abs(VTKdata.totField);
                VTKdata.totField = real(VTKdata.totField);
            end
        end
        if mod(j,2) == 0 
            m = m + 1;
        end
%         VTKdata.displacement = displacement;
        VTKdata.nodes = nodes{j};
        VTKdata.visElements = visElements{j};
        VTKdata.omega = omega;
        if plotInTimeDomain
            VTKoptions.T = T;
        end
        VTKoptions.N = N_fine;
        disp(['Time spent on assembling data ' num2str(toc) ' seconds.'])
        tic
        if mod(j,2) || computeForSolidDomain
            makeVTKfile(VTKdata, VTKoptions);
        end
        disp(['Time spent on making VTK file ' num2str(toc) ' seconds.'])
    end
end
       