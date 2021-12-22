function createParaviewFiles_exact2(extraPts, vtfFileName, options)

SSBC = options.SSBC;
HWBC = options.HWBC;
R = options.R;
R_o = options.R_o;
R_a = options.R_a;
r_s = options.r_s;
% P_inc = options.P_inc;
d_vec = options.d_vec;

M = length(R_o);
plotTimeOscillation = 0;
addBackSideWall = 0;

rho_f = options.rho_f;
omega = options.omega;
c_f = options.c_f;
VFBC = options.VFBC;
noDomains = 2*M+1-SSBC-2*HWBC-VFBC;
nodes = cell(noDomains,1);
visElements = cell(noDomains,1);
computeForSolidDomain = 0;
tic
for m = 1:M+1-SSBC-HWBC-VFBC
    if m ~= M+1 && ~(HWBC && m == M) && computeForSolidDomain
        if SSBC && m == M
            solid = getSSBCData(R_o(m),'Xaxis');
            t_solid = R_o(m);
            R_solid = R_o(m);
        else
            solid = getSphericalShellData(R(m), R_o(m),'Xaxis');
            t_solid = R_o(m)-R(m);
            R_solid = R_o(m);
        end
        Xi = solid.knots{1};
        Eta = solid.knots{2};
        Zeta = solid.knots{3};
        
        [nodes{2*m}, ~, visElements{2*m}] = buildVisualization3dMesh_new3(Xi, Eta, Zeta, ...
          round(R_solid*pi/2*extraPts/R_a), round(R_solid*pi/2*extraPts/R_a), round(t_solid*extraPts/R_a), solid);
    end
    
    if m == 1
        [visElements{2*m-1}, nodes{2*m-1}] = meshRectangleWcircHole([-1.3*R_a, -R_a],[1.3*R_a, R_a],R_o(1),round(extraPts*1.3*R_a));
        nodes{2*m-1} = [nodes{2*m-1}, zeros(size(nodes{2*m-1},1),1)];
    elseif m < M+1
        [visElements{2*m-1}, nodes{2*m-1}] = mesh2DDonut(R(m-1),R_o(m),round(1.18*extraPts*2*R(m-1)));
        nodes{2*m-1} = [nodes{2*m-1}, zeros(size(nodes{2*m-1},1),1)];
    else
        [visElements{2*m-1}, nodes{2*m-1}] = mesh2DDisk(R_o(m),round(1.18*extraPts*2*R_o(m)));
        nodes{2*m-1} = [nodes{2*m-1}, zeros(size(nodes{2*m-1},1),1)];
    end
end
if addBackSideWall
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
if 0
    N = 1;
    Nq = 1;
    omega = 1;
    T = NaN;
    startIdx = 1;
else
    f_c = 1500;
    N = 2^9; % 200
    M = 2*N;
    T = 60/f_c; %N/M;
    B = N/T; % bandwidth
    
    f_R = B/2;
    df = 1/T;
    f = linspace(0,f_R-df,N/2);
    omega = 2*pi*f;
    options.omega = omega(2:end);
    if options.usePlaneWave
        startIdx = 900;
    elseif options.usePointChargeWave
        startIdx = 1000;
    end
end
type = 1;
totnpts = npts*N*4
omega_c = 2*pi*f_c;
options.P_inc = @(omega) P_inc_(omega,omega_c,type);
tic
data = e3Dss(nodes, options);
disp(['Time spent on computing solution: ' num2str(toc) ' seconds.'])
m = 1;
for j = 1:length(nodes)
    if ~(VFBC && j == length(nodes))
        tic
        clear VTKdata
        if mod(j,2) == 0 && computeForSolidDomain
            VTKoptions = struct('name',[vtfFileName 'solid' num2str(m)], 'celltype', 'VTK_HEXAHEDRON', 'plotTimeOscillation', plotTimeOscillation, ...
                        'plotSphericalRadialDisplacement',0, 'plotDisplacementVectors',0,'plotSphericalStress_rr',0); 
            if VTKoptions.plotDisplacementVectors 
                displacement = zeros(size(nodes{j},1),3,length(omega));
                for i = 2:length(omega)
                    displacement(:,:,i) = [data(m).u_x(:,i-1) data(m).u_y(:,i-1) data(m).u_z(:,i-1)];
                end
            end
            if VTKoptions.plotSphericalStress_rr 
                VTKdata.stress = zeros(size(nodes{j},1),6,length(omega));
                if options.calc_sigma_xx || options.calc_sigma_yy || options.calc_sigma_zz || options.calc_sigma_yz || options.calc_sigma_xz || options.calc_sigma_xy 

                    for i = 2:length(omega)
                        VTKdata.stress(:,:,i) = [data(m).sigma_xx(:,i-1) data(m).sigma_yy(:,i-1) data(m).sigma_zz(:,i-1) data(m).sigma_yz(:,i-1) data(m).sigma_xz(:,i-1) data(m).sigma_xy(:,i-1)];
                    end
                end
                VTKdata.stress = 2/T*real(fft(VTKdata.stress,M,3));
                temp = VTKdata.stress;
                VTKdata.stress(:,:,1:M-startIdx+1) = temp(:,:,startIdx:end);
                VTKdata.stress(:,:,M-startIdx+2:end) = temp(:,:,1:startIdx-1);
            end
        elseif mod(j,2) ~= 0
            VTKoptions = struct('name',[vtfFileName 'fluid' num2str(m)], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', plotTimeOscillation, 'plotDisplacementVectors', 0, ...
                             'plotSphericalRadialDisplacement',0, 'plotTotField', 1);
            if VTKoptions.plotDisplacementVectors 
                displacement = zeros(size(nodes{j},1),3,length(omega));
            end
            VTKdata.totField = zeros(size(nodes{j},1),1,length(omega));
            if m == 1
%                 VTKoptions.plotScalarField = 0;
%                 VTKoptions.plotP_inc = 0;
                VTKdata.P_inc = zeros(size(nodes{j},1),1,length(omega));
%                 VTKdata.scalarField = zeros(size(nodes{j},1),1,length(omega));
                for i = 2:length(omega)
                    if options.usePlaneWave
                        k = omega(i)/c_f(1);
                        k_vec = options.d_vec*k;
                        p_inc = @(v) P_inc_(omega(i),omega_c,type).*exp(1i*dot3(v, k_vec));
%                         gP_inc = @(v) 1i*p_inc(v)*k_vec.';
                    elseif options.usePointChargeWave
                        k = omega(i)/c_f(1);
                        r = @(y) norm2(repmat(-r_s*d_vec.',size(y,1),1)-y);
                        Phi_k = @(y) exp(1i*k*r(y))./(4*pi*r(y));     
                        p_inc = @(y) 4*pi*r_s*P_inc_(omega(i),omega_c,type)*Phi_k(y);
                    end
%                     gScalarField = [data(m).dpdx(:,i) data(m).dpdy(:,i) data(m).dpdz(:,i)];
%                     displacement(:,:,i) = (gScalarField+gP_inc(nodes{1}))/(rho_f(m)*omega(i)^2);
                    VTKdata.P_inc(:,:,i) = p_inc(nodes{1});
                    VTKdata.totField(:,:,i) = data(m).p(:,i-1) + VTKdata.P_inc(:,:,i);
%                     VTKdata.scalarField(:,:,i) = data(m).p(:,i);
                end
                VTKdata = rmfield(VTKdata,'P_inc');
%                 VTKdata.P_inc = 2/T*real(fft(VTKdata.P_inc,Nq,3));
%                 temp = VTKdata.P_inc;
%                 VTKdata.P_inc(:,:,1:Nq-startIdx+1) = temp(:,:,startIdx:end);
%                 VTKdata.P_inc(:,:,Nq-startIdx+2:end) = temp(:,:,1:startIdx-1);
%                 
%                 VTKdata.scalarField = 2/T*real(fft(VTKdata.scalarField,Nq,3));
%                 temp = VTKdata.scalarField;
%                 VTKdata.scalarField(:,:,1:Nq-startIdx+1) = temp(:,:,startIdx:end);
%                 VTKdata.scalarField(:,:,Nq-startIdx+2:end) = temp(:,:,1:startIdx-1);
            else
                VTKoptions.plotScalarField = false;
                VTKoptions.plotP_inc = false;
                for i = 2:length(omega)
%                     gScalarField = [data(m).dpdx(:,i) data(m).dpdy(:,i) data(m).dpdz(:,i)];
%                     displacement(:,:,i) = gScalarField/(rho_f(m)*omega(i)^2);
                    VTKdata.totField(:,:,i) = data(m).p(:,i-1);
                end
            end
            VTKdata.totField = 2/T*real(fft(VTKdata.totField,M,3));
            temp = VTKdata.totField;
            VTKdata.totField(:,:,1:M-startIdx+1) = temp(:,:,startIdx:end);
            VTKdata.totField(:,:,M-startIdx+2:end) = temp(:,:,1:startIdx-1);
%             VTKdata.gradient = zeros(size(nodes{j},1),3,length(omega));
%             for i = 1:length(omega)
%                 VTKdata.gradient(:,:,i) = gScalarField;
%             end
%             VTKdata.gradient = 2/T*real(fft(VTKdata.gradient,Nq,3));
%             temp = VTKdata.gradient;
%             VTKdata.gradient(:,:,1:Nq-startIdx+1) = temp(:,:,startIdx:end);
%             VTKdata.gradient(:,:,Nq-startIdx+2:end) = temp(:,:,1:startIdx-1);
        end
        if mod(j,2) == 0 
            m = m + 1;
        end
%         VTKdata.displacement = displacement;
        VTKdata.nodes = nodes{j};
        VTKdata.visElements = visElements{j};
        VTKdata.omega = omega;
        VTKoptions.T = T;
        VTKoptions.N = M;
        disp(['Time spent on FFT ' num2str(toc) ' seconds.'])
        tic
        if mod(j,2) || computeForSolidDomain
            makeVTKfile(VTKdata, VTKoptions);
        end
        disp(['Time spent on making VTK file ' num2str(toc) ' seconds.'])
    end
end
       