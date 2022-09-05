function studies = getTask_illustratePML()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 19 to 22 in Mi2021ilc
% Mi2021ilc is available at https://doi.org/10.1016/j.cma.2021.113925

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

saveStudies = false;

saveStudies = false;
misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'illustratePML';

c_x = 3;
c_y = 1;
c_z = 2;
% c_x = 1;
% c_y = 1;
% c_z = 1;
    
c_f = 1500;   % Speed of sound in outer fluid
k = 1;                 % Wave number for Simpson2014aib
misc.omega = c_f*k;    % Angular frequency
msh.parm = 1;
% c_x = c_z;
varCol{1} = struct('media', 'fluid', ...
                   'c_x',   c_x, ...
                   'c_y',   c_y, ...
                   'c_z',   c_z, ...
                   'c_f', 1500, ...
                   'rho', 1000);
msh.meshFile = 'createNURBSmesh_EL';
misc.preProcessOnly = true;
misc.compute_h_max  = true;        % Compute h_max and derived quantities like nepw (number of elements per wavelength)

misc.checkNURBSweightsCompatibility = 0;
prePlot.plotGeometryInfo = 0;       % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
prePlot.plotFullDomain   = 1;        % Plot volumetric domains
at = false(3,2);
at(3,2) = true;
prePlot.plotAt = at;
prePlot.view             = [0,90];
prePlot.view             = getView(3);
prePlot.plotSubsets      = {};
% prePlot.plotSubsets      = {'xy'};
% prePlot.plotSubsets      = {'xy','Gamma','Gamma_a','Gamma_b'};
prePlot.plotParmDir      = 0;
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
prePlot.plotControlPolygon = 0;       % Plot the control polygon for the NURBS mesh
prePlot.abortAfterPlotting = true;       % Abort simulation after pre plotting
prePlot.coarseLinearSampling = false;
% prePlot.colorFun = @(v) abs((v(:,1)/c_x).^2 + (v(:,2)/c_y).^2 + (v(:,3)/c_z).^2 - 1);
prePlot.colorFun = @(v,n,dXdzeta) log10(norm2(n-dXdzeta));
prePlot.colorParmDirs = {getColor(4),getColor(4),getColor(4)};
prePlot.useCamlight = false;
prePlot.quiverScale = 0.1;
prePlot.quiverLineWidth = 0.5;
prePlot.LineWidth = 1;
prePlot.QoI = @(v,n,dXdzeta) dXdzeta./norm2(dXdzeta);
prePlot.QoI_ref = @(v,n,dXdzeta) n;
warning('off','NURBS:weights')

msh.parm = 1;
msh.explodeNURBS = 0;   % Create patches from all C^0 interfaces
msh.alignWithAxis = 'Xaxis';
% msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
msh.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
msh.refineThetaOnly = 0;

misc.method = {'PML'};
misc.formulation = {'GSB'};

ffp.calculateFarFieldPattern = false;     % Calculate far field pattern
ffp.beta_s = 0;   
ffp.alpha_s = 0;   

misc.omega = 1;
misc.BC = 'SHBC';
msh.degree = 2;

pml.t = 1;         % thickness of PML
pml.dirichlet = true;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
misc.r_a = c_x*1.1; 
% misc.extraGP = [10,10,10];
% ffp.extraGP = [50,0,0];

loopParameters = {'msh.M','msh.degree','misc.coreMethod'};

pml.sigmaType = 3;   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n, sigmaType = 3: sigma(xi) = C/(1-xi)^n
pml.n = 1;
% misc.coreMethod = {'linear_FEM'};
misc.coreMethod = {'IGA'};
% misc.coreMethod = {'C0_IGA'};
switch misc.coreMethod{1}
    case 'linear_FEM'
        pml.X_bApprox = 'interp';
        msh.M = 4; % 4
        prePlot.resolution = [20,2,0];
    otherwise
        pml.X_bApprox = 'BA';
%         pml.X_bApprox = 'interp';
        msh.M = 3; % 4
        prePlot.resolution = [123,13,0];
%         prePlot.resolution = [2,2,0];
        prePlot.qstep = 4;
end
msh.M = 7; % 4
prePlot.resolution = [800,400,0];
prePlot.resolution = [100,50,0];
pml.refinement = @(M) 1;
varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 1];
prePlot.addCommands = @() addCommands(pml,misc,varCol);
% collectIntoTasks


misc.coreMethod = {'linear_FEM','IGA'};

% postPlot(1).xname        	= 'dofs';
postPlot(1).xname        	= 'surfDofs';
postPlot(1).yname        	= 'QoIError';    % Examples include: 'p_Re', 'p_Im', 'abs_p', 'TS', 'error_pAbs', 'error_p', 'surfaceError', 'energyError', 'L2Error', 'H1Error', 'H1sError'
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).addSlopes       = true;
postPlot(1).axisType      	= 'loglog';
postPlot(1).lineStyle    	= '*-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).noXLoopPrms   	= 1;

M_0 = 7;
msh.M = 1:M_0; % 4
misc.coreMethod = 'IGA';
pml.X_bApprox = 'BA';
msh.degree = 2:4;
collectIntoTasks

msh.M = 1:M_0; % 4
misc.coreMethod = 'C0_IGA';
collectIntoTasks

misc.coreMethod = 'hp_FEM';
collectIntoTasks

misc.coreMethod = 'h_FEM';
pml.X_bApprox = 'interp';
collectIntoTasks

msh.M = 1:M_0; % 4
misc.coreMethod = 'linear_FEM';
collectIntoTasks


function addCommands(pml,misc,varCol)
t = pml.t;
t_fluid = misc.r_a - varCol{1}.c_x;
c_x = varCol{1}.c_x + t_fluid;
c_y = varCol{1}.c_y + t_fluid;
c_z = varCol{1}.c_z + t_fluid;

X_a_1 = @(theta,phi) c_x.*sin(theta).*cos(phi);
X_a_2 = @(theta,phi) c_y.*sin(theta).*sin(phi);
X_a_3 = @(theta,phi) c_z.*cos(theta).*ones(size(phi));
q = @(theta,phi) sqrt(c_x^2*c_y^2*cos(theta).^2 + sin(theta).^2.*(c_x^2*c_z^2*sin(phi).^2 + c_y^2*c_z^2*cos(phi).^2));
X_b_1 = @(theta,phi) X_a_1(theta,phi) + (c_y.*c_z.*t./q(theta,phi)).*sin(theta).*cos(phi);
X_b_2 = @(theta,phi) X_a_2(theta,phi) + (c_x.*c_z.*t./q(theta,phi)).*sin(theta).*sin(phi);
X_b_3 = @(theta,phi) X_a_3(theta,phi) + (c_x.*c_y.*t./q(theta,phi)).*cos(theta);

if false
    theta = linspace(0,pi/2,10000);
    X_a = [c_x*cos(theta);c_y*sin(theta)];
    X_b = X_a - t*[(-1)*c_x.*c_y.^2.*cos(theta)./(sqrt(c_x.^2.*c_y.^2./(c_x.^2.*sin(theta).^2 - c_y.^2.*sin(theta).^2 + c_y.^2).^2).*(c_x.^2.*sin(theta).^2 + c_y.^2.*cos(theta).^2).^(3/2)); (-1)*c_x.^2.*c_y.*sin(theta)./(sqrt(c_x.^2.*c_y.^2./(c_x.^2.*sin(theta).^2 - c_y.^2.*sin(theta).^2 + c_y.^2).^2).*(c_x.^2.*sin(theta).^2 + c_y.^2.*cos(theta).^2).^(3/2))];
    hold on
    plot(X_a(1,:),X_a(2,:),'DisplayName','Exact X_a')
    plot(X_b(1,:),X_b(2,:),'DisplayName','Exact X_b')
    % 
    max(abs(X_b(1,:)-X_b_1(fliplr(theta),phi)))
    max(abs(X_b(2,:)-X_b_3(fliplr(theta),phi)))
end

if false
    phi = linspace(0,pi/2,10000);
    theta = pi/2;
    plot3(X_a_1(theta,phi),X_a_2(theta,phi),X_a_3(theta,phi),'DisplayName','Exact X_a','color','red','LineWidth',1)
    plot3(X_b_1(theta,phi),X_b_2(theta,phi),X_b_3(theta,phi),'DisplayName','Exact X_b','color','red','LineWidth',1) 
%     export_fig('~/Dropbox/Apps/Overleaf/PML/graphics/PML_IGA', '-pdf', '-transparent')
%     export_fig('~/Dropbox/Apps/Overleaf/PML/graphics/PML_linear_FEM', '-pdf', '-transparent')
else
    caxis([-5,0])
%     colorbar off

%     export_fig('~/Dropbox/Apps/Overleaf/PML/graphics/PML_IGA_3D', '-png', '-transparent', '-r300')
%     export_fig('~/Dropbox/Apps/Overleaf/PML/graphics/PML_linear_FEM_3D', '-png', '-transparent', '-r300')
end



if false
    % figure
    theta_0 = linspace(0,pi,100);
    phi_0 = linspace(0,2*pi,200);
    [theta,phi] = meshgrid(theta_0,phi_0);
    surf(X_a_1(theta,phi),X_a_2(theta,phi),X_a_3(theta,phi),'FaceColor', getColor(1), 'EdgeColor','none','LineStyle','none','FaceAlpha',1, 'DisplayName','Exact X_a')
    hold on
    surf(X_b_1(theta,phi),X_b_2(theta,phi),X_b_3(theta,phi),'FaceColor', getColor(1), 'EdgeColor','none','LineStyle','none','FaceAlpha',0.5, 'DisplayName','Exact X_b')
    axis equal
    % camlight
    % legend show
end
