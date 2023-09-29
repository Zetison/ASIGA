clear all
% close all

addpath(genpath('../export_fig'))
no2Dpoints = 1000;
% startMatlabPool



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot all shapes
% R = 1;
% models = {'Disk','Cone','WineGlass','Ellipsoid','Torus','Quadrilateral','Prism','Cube','Cylinder','CE','HalfSphere',...
%            'BeTSSiM1','BeTSSiM2','BeTSSiM3','BeTSSiSmoothM3','BeTSSiM4','BeTSSiM5','BeTSSiPH','MockShell','QuarterDisk','Barrel','BeTSSiM5','ScordelisLoRoof','Cube'};
% models = {'Mutter'};
% for model = models
%     for parm = 2
% %         close all
%         options.parm = parm;
% %         options.theta = 2*pi;r
% %         options.theta_eta = 0.9*2*pi;
% %         options.R = 2;
% %         options.t = 0.1;
%         nurbs = eval(['get' model{1} 'Data(options)']);
%         resolution = 32*[1,1,1];
%         M = 2;
%         refLength = 1;
%         geometry = getTopology(nurbs);
% %         nurbs = autoRefineNURBS(nurbs,geometry.topology.connection,refLength/2^(M-1));
% %         nurbs = elevateNURBSdegree(nurbs,ones(1,nurbs{1}.d_p));
%         figure
%         plotNURBS(nurbs,'resolution',resolution,'plotControlPolygon',0,...
%             'plotNormalVectors',0,'plotParmDir',0,'plotWeights',0,'coarseLinearSampling',false);
%         axis equal
%         grid off
%         axis off
%         colorbar
%         view(getView(0))
% %         view(getView(1))
% %         camlight
%         material dull
%         ax = gca;               % get the current axis
%         ax.Clipping = 'off';    % turn clipping off
% 
% %         write_g2(nurbs,['NURBSgeometries/g2files/' model{1} '_parm' num2str(parm)])
% %         axis on
% %         figureFullScreen(gcf)
% %         export_fig(['../../graphics/ASIGAmodels/' model{1} '_' num2str(parm)], '-png', '-transparent', '-r200')
%     end
% end
% % nurbs = translateNURBS(nurbs,[-L/2,0,0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test mutter
% close all
% figure
% nurbs = read_g2('../IGA-geometries/Mutter/Mutter.g2');
% axis equal
% view(1,40)
% camlight
% % plotNURBS(nurbs([6,15,13]),'plotControlPolygon',true)
% % plotNURBS(nurbs,'plotControlPolygon',true,'plotParmDir',1)
% 
% figure
% nurbs = getMutterData();
% axis equal
% view(1,40)
% camlight
% plotNURBS(nurbs,'plotControlPolygon',true,'plotParmDir',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create random cubic hole
% nurbs_vol = read_g2('../IGA-geometries/Chess/Capablanca.g2');
% clear all
% close all
% nurbs = getPrismData('L',[2,4,1],'x_0',[0,0,0]);
% nurbs = elevateNURBSdegree(nurbs,2);
% nurbs{1}.coeffs(1,:,:,:) = nurbs{1}.coeffs(1,:,:,:).*abs(1+0.1*nurbs{1}.coeffs(2,:,:,:));
% nurbs{1}.coeffs(2,:,:,:) = nurbs{1}.coeffs(2,:,:,:).*abs(1+0.1*nurbs{1}.coeffs(3,:,:,:));
% nurbs{1}.coeffs(3,:,:,:) = nurbs{1}.coeffs(3,:,:,:).*abs(1+0.1*nurbs{1}.coeffs(2,:,:,:));
% nurbs{1}.coeffs = nurbs{1}.coeffs + rand(size(nurbs{1}.coeffs(:,:,:,:)))*0.05;
% nurbs = insertKnotsInNURBS(nurbs,{{[ones(1,3)/3, 2*ones(1,3)/3],[ones(1,3)/3, 2*ones(1,3)/3], []}});
% nurbs = explodeNURBS(nurbs);
% nurbs = nurbs([1:4,6:9]);
% figure, plotNURBS(nurbs)
% axis equal
% view(getView())
% write_g2(nurbs,'NURBSgeometries/g2files/randomCubicHole.g2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test loftNURBS
% close all
% nurbs = getRectangleData();
% nurbs = ensure3DNURBS(nurbs);
% nurbs = insertKnotsInNURBS(nurbs,{0.5, [0.3, 0.7]});
% nurbs = explodeNURBS(nurbs);
% nurbs = elevateNURBSdegree(nurbs,1);
% nurbs2 = translateNURBS(nurbs,[0,-0.5,1]);
% nurbs3 = translateNURBS(nurbs,[0,0.5,2]);
% nurbs4 = translateNURBS(nurbs,[0.5,0,3]);
% nurbsCol = {nurbs,nurbs2,nurbs3,nurbs4};
% for i = 1:numel(nurbsCol)
%     for j = 1:numel(nurbsCol{i})
%         n = round(rand(1)*4);
%         m = round(rand(1)*4);
%         nurbsCol{i}(j) = insertKnotsInNURBS(nurbsCol{i}(j),{rand(1,n), rand(1,m)});
%         nurbsCol{i}{j}.coeffs = nurbsCol{i}{j}.coeffs + rand(size(nurbsCol{i}{j}.coeffs))*0.05;
%     end
%     plotNURBS(nurbsCol{i})
% end
% axis equal
% view(1,40)
% camlight
% nurbsVol = loftNURBS(nurbsCol,3,1);
% 
% plotNURBS(nurbsVol)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Capablanca
% nurbs_vol = read_g2('../IGA-geometries/Chess/Capablanca.g2');
% if false
%     nurbs_vol = nurbs_vol([121,280]);
%     coeffs1 = nurbs_vol{2}.coeffs(:,:,1,end);
%     coeffs2 = nurbs_vol{1}.coeffs(:,:,1,1);
%     knots = {[0,0,0,1,1,1]};
%     nurbs_circ1 = createNURBSobject(coeffs1,knots);
%     nurbs_circ2 = createNURBSobject(coeffs2,knots);
%     
%     xi = linspace(0,1,1000).';
%     C = evaluateNURBS(nurbs_circ1{1},xi);
%     NURBS_error = abs(norm2(C(:,1:2)+[3,3])-0.7528);
%     figure
%     semilogy(xi,NURBS_error)
%     hold on
%     C = evaluateNURBS(nurbs_circ2{1},xi);
%     NURBS_error = abs(norm2(C(:,1:2)+[3,3])-0.7528);
%     semilogy(xi,NURBS_error)
%     
%     close all
%     axis equal
%     plotNURBS([nurbs_circ1,nurbs_circ2],'plotControlPolygon',1,'plotParmDir',0,'resolution',[100,100,100]);
% end
% nurbs_vol = cleanNURBS(nurbs_vol,[],1e-6);
% % nurbs_vol = nurbs_vol([109:128,273:298]);
% % nurbs_vol = nurbs_vol([13,28]);
% nurbs_vol = nurbs_vol(setdiff(1:numel(nurbs_vol),[109:128,299:318]));
% geometry = getTopology(nurbs_vol);
% if false
%     nurbs = extractFreeSurface(nurbs_vol,'topologysets',geometry.topologysets);
%     save('../../results/ASIGA/Capablanca/nurbs.mat','nurbs')
%     save('../../results/ASIGA/Capablanca/nurbs_woRook.mat','nurbs')
% end
% refLength = 2;
% M = 5;
% % nurbs_vol = autoRefineNURBS(nurbs_vol,geometry.topology.connection,refLength/2^(M-1));
% if 1
%     close all
%     axis equal
%     color = getColor();
%     plotNURBS(nurbs_vol,'plotControlPolygon',0,'plotParmDir',0,'color',color);
% %     plotNURBS(nurbs_vol,'plotControlPolygon',1,'plotParmDir',0);
%     drawnow
%     view(getView());
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test surfaceToVolume
% close all
% % model = 'Capablanca_woRooks';
% % model = 'Capablanca';
% model = 'BCA';
% % model = 'S1_interior';
% % model = 'S1';
% % model = 'Mutter';
% switch model
%     case 'Mutter'
%         nurbs_vol = read_g2('../IGA-geometries/Mutter/Mutter.g2');
%         nurbs_vol = cleanNURBS(nurbs_vol,[],1e-2);
%         nurbs = extractFreeSurface(nurbs_vol);
% 
%         t = 1.1;
%         sharpAngle = 120*pi/180; % Threshold for a "sharp" angle
%     case 'Capablanca'
%         load('../../results/ASIGA/Capablanca/nurbs.mat')
% 
%         t = 1.1;
%         sharpAngle = 120*pi/180; % Threshold for a "sharp" angle
%         uiopen('../../results/ASIGA/Capablanca/CAD.fig',1)
%     case 'Capablanca_woRooks'
%         load('../../results/ASIGA/Capablanca/nurbs_woRook.mat')
%         t = 1.1;
%         sharpAngle = 120*pi/180; % Threshold for a "sharp" angle
%         uiopen('../../results/ASIGA/Capablanca/CADwoRook.fig',1)
%     case 'BCA'
%         load('BeTSSi_BCA_p2_unRefined.mat')
%         t = 1.1;
%         % sharpAngle = 100*pi/180; % Threshold for a "sharp" angle
%         sharpAngle = 120*pi/180; % Threshold for a "sharp" angle
%         % sharpAngle = 125*pi/180; % Threshold for a "sharp" angle
%         % sharpAngle = 135*pi/180; % Threshold for a "sharp" angle
%         % sharpAngle = 160*pi/180; % Threshold for a "sharp" angle
%         % sharpAngle = 161*pi/180; % Threshold for a "sharp" angle % In order to include connections over depth rudders
%         % sharpAngle = 173*pi/180; % Threshold for a "sharp" angle
%         uiopen('../../results/ASIGA/BCA/CAD.fig',1)
%     case 'S1_interior'
%         nurbs = getEllipsoidData('parm',2,'S2V_algorithm',{'A1_2'});
%         nurbs = flipNURBSparametrization(nurbs,1);
%         t = 0.4;
%         sharpAngle = 140*pi/180; % Threshold for a "sharp" angle
%     case 'S1'
%         nurbs = getEllipsoidData('parm',1);
%         t = 0.4;
%         sharpAngle = 140*pi/180; % Threshold for a "sharp" angle
% end
% if ~strcmp(model,'BCA') && ~strcmp(model,'Capablanca') && ~strcmp(model,'Capablanca_woRooks')
%     close all
%     axis equal
%     color = getColor();
%     plotNURBS(nurbs,'plotControlPolygon',0,'plotParmDir',0,'color',color);
%     drawnow
%     view(getView());
% end
% if 0
%     close all
%     axis equal
%     color = getColor();
%     plotNURBS(nurbs_vol,'plotControlPolygon',0,'plotParmDir',0,'color',color);
% %     plotNURBS(nurbs,'plotControlPolygon',0,'plotParmDir',0,'colorFun',@(x) log10(abs(norm2(x(:,[1,3]))-6)));
%     drawnow
%     view(getView());
% end
% if 1
%     geometry = getTopology(nurbs);
%     nurbs_vol = surfaceToVolume(nurbs,'t',t,'sharpAngle',sharpAngle);
% end

% close all
% load('BeTSSi_BCA_p2_unRefined.mat')
% % set(0, 'DefaultFigureRenderer', 'opengl');
% tic
% plotNURBS(nurbs,'plotControlPolygon', true);
% toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test convexifyNURBS
% close all
% % nurbs = getBeTSSiM4Data('t',1);
% % nurbs = extractFreeSurface(nurbs);
% nurbs = getEllipsoidData();
% axis equal
% plotNURBS(nurbs,'plotControlPolygon',true,'coarseLinearSampling',false)
% camlight
% view(120,40)
% nurbs2 = convexifyNURBS(nurbs);
% hold on
% plotNURBS(nurbs2,'plotControlPolygon',true,'plotObject',1,'coarseLinearSampling',false)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test collapseToEdge
% close all
% camlight
% nurbs = getPrismData();
% % nurbs = getHexEllipsoidData();
% for midx = 1:6
%     for midx_midx = 1:4
%         nurbs2 = collapseToEdge(nurbs{1},midx,midx_midx);
% 
%         nurbsSurf = subNURBS(nurbs);
%         h1 = plotNURBS(nurbsSurf,'plotControlPolygon',true);
%         h1.objectHandle(midx).FaceColor = 'red';
%         view(getView())
%         axis equal
%         pause(1)
%         hold off
%         h2 = plotNURBS(nurbs2,'plotControlPolygon',true);
%         delete(h1.objectHandle)
%         delete(h1.elementEdgesHandle)
%         delete(h1.controlPolygonHandle)
%         view(getView())
%         axis equal
%         pause(1)
%         hold off
%         delete(h2.objectHandle)
%         delete(h2.elementEdgesHandle)
%         delete(h2.controlPolygonHandle)
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test C1 NURBS curve
% close all
% clear all
% 
% syms theta x
% % theta = 30*vpa('pi')/180;
% a = cos(theta/2);
% PI = vpa('pi');
% phi = (PI-theta)/2;
% w = cos(theta/2);
% third = int(x^2,0,1);
% fift = int(x^4,0,1);
% sixt = int(x^5,0,1);
% coeffs = [cos(phi),sin(phi),1;
%           0,1/w,-w;
%           -cos(phi),sin(phi),1].';
% knots = {[0,0,0,1,1,1]};
% nurbs = createNURBSobject(coeffs,knots);
% % nurbs = insertKnotsInNURBS(nurbs,{[third,2*third]});
% % nurbs = insertKnotsInNURBS(nurbs,{[fift,2*fift,3*fift,4*fift]});
% nurbs = insertKnotsInNURBS(nurbs,{[sixt,2*sixt,3*sixt,4*sixt,5*sixt]});
% 
% % theta_num = 30*pi/180;
% theta_num = 1*pi/180;
% nurbs{1}.coeffs = subs(nurbs{1}.coeffs,'theta',theta_num);
% nurbs{1}.coeffs = subs(nurbs{1}.coeffs,'pi',pi);
% nurbs{1}.knots{1} = subs(nurbs{1}.knots{1},'theta',theta_num);
% nurbs{1}.knots{1} = subs(nurbs{1}.knots{1},'pi',pi);
% nurbs{1}.coeffs = double(nurbs{1}.coeffs);
% nurbs{1}.knots{1} = double(nurbs{1}.knots{1});
% xi = linspace(0,1,1000).';
% plotNURBS(nurbs,'plotControlPolygon',true,'resolution',numel(xi))
% axis equal
% 
% C = evaluateNURBS(nurbs{1},xi);
% NURBS_error = abs(norm2(C)-1);
% figure
% semilogy(xi,NURBS_error)
% return
% w = [1.5,0.5,0.5,1/3];
% if 1
%     options = optimset('Display','iter','TolX',10*eps,'TolFun',10*eps,'MaxFunEvals',1e5,'MaxIter',1e5);
%     w = fminsearchbnd(@(w) f_obj(w,xi),w,[0,0,0,0.1],[10,10,10,0.4],options);
% end
% [coeffs,knots] = getNURBScoeffAndKnots(w);
% nurbs = createNURBSobject(coeffs,knots);
% plotNURBS(nurbs,'plotControlPolygon',true,'resolution',numel(xi))
% axis equal
% 
% C = evaluateNURBS(nurbs{1},xi);
% NURBS_error = abs(norm2(C)-1);
% figure
% semilogy(xi,NURBS_error)
% 
% function [coeffs,knots] = getNURBScoeffAndKnots(w)
% 
% x = w(1);
% y = -(x^2+1)/(x^2-1);
% coeffs = [0,1,1;
%           -x,1,w(2);
%           0,y,w(3);
%           x,1,w(2);
%           0,1,1].';
% knots = {[0,0,0,w(4),1-w(4),1,1,1]};
% 
% end
% 
% function NURBS_error = f_obj(w,xi)
% [coeffs,knots] = getNURBScoeffAndKnots(w);
% nurbs = createNURBSobject(coeffs,knots);
% C = evaluateNURBS(nurbs{1},xi);
% NURBS_error = norm(abs(norm2(C)-1));
% end

%% NURBS parametrization of circle
% close all
% knots = {[0,0,0,1/2,1,1,1]};
% coeffs = [1,0,1;
%           1,1,0.5;
%           -1,1,0.5;
%           -1,0,1].';
% nurbs = createNURBSobject(coeffs,knots);
% nurbs = elevateNURBSdegree(nurbs,1);
% nurbs = insertKnotsInNURBS(nurbs,6);
% figure(1)
% plotNURBS(nurbs,'resolution',1000,'plotControlPolygon',true);
% axis equal
% nurbs{1}.knots{1}
% 
% rotationMatrix(2*atan(2/10),'Zaxis')*[-2/10,1,0].'
% rotationMatrix(2*atan(1/4),'Zaxis')*[-13/15,6/10,0].'


% 
% 
% % nurbs{1}.coeffs(:,6) = [-13/15,6/10,0.5];
% % nurbs{1}.coeffs(:,7) = [-79/75,-4/100,0.5];
% figure(1)
% plotNURBS(nurbs,'resolution',1000,'plotControlPolygon',true);
% 
% figure(2)
% C = abs(norm2(evaluateNURBS(nurbs{1},xi))-1);
% semilogy(xi,C);
% nurbs = getArcData('theta',120*pi/180,'Xi',[0,0,0,1,1,1]);
% plotNURBS(nurbs,'resolution',10000,'plotControlPolygon',true);
% axis equal
% nurbs{1}.coeffs(3,2) = -nurbs{1}.coeffs(3,2);
% plotNURBS(nurbs,'resolution',10000,'plotControlPolygon',true);
% figure(2)
% xi = linspace(0,1,1000).';
% C = abs(norm2(evaluateNURBS(nurbs{1},xi))-1);
% semilogy(xi,C);
% hold on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test meanRatioJacobian
% close all
% theta = pi/4;
% % theta = pi/20;
% nurbs = getRectangleData('L',[1,4]);
% nurbs = insertKnotsInNURBS(nurbs,[0,3]);
% % nurbs = getHalfDiskData('theta3',theta,'parm',3);
% % nurbs = ensure3DNURBS(nurbs);
% % nurbs = makeUniformNURBSDegree(nurbs);
% figure
% [~,maxC_patch,minC_patch] = plotNURBS(nurbs,'plotControlPolygon',0,'plotParmDir',0,'colorFun',@(xi,nurbs,b,c) meanRatioJacobian(nurbs,xi),'resolution',[100,100],'coarseLinearSampling',false);
% colorbar
% axis equal
% 
% clim([0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate and save BeTSSi BCA meshes
% % close all
% % nurbs = getBCAData(2,false);
% % plotNURBS(nurbs)
% % axis equal
% % view(getView())
% % camlight
% % return
% for p = 2:15
%     nurbs = getBCAData(p,false);
%     save(['NURBSgeometries\BCAdata\BeTSSi_BCA_p' num2str(p) '_unRefined.mat'],'nurbs')
%     fprintf('Completed p = %d\n',p)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test normalBasedSurface2volume
% close all
% C = 1;
% C = [1,2,3];
% nurbs = getEllipsoidData('C',C,'parm',2);
% M = 1;
% nurbs = insertKnotsInNURBS(nurbs,(2^(M-1)-1)*[1,1]);
% t_PML = 2;
% 
% I = 1:6;
% % I = 1;
% prePlot.colorFun = @(v) log10(abs(norm2(v./C)-1));
% prePlot.resolution = [50,50,0];
% prePlot.plotControlPolygon = true;
% % plotNURBS(nurbs(I), prePlot);
% C = C+t_PML;
% [nurbsVol,nurbs2] = normalBasedSurface2volume(nurbs,t_PML);
% prePlot.colorFun = @(v) log10(abs(norm2(v./C)-1));
% plotNURBS(nurbs2(I), prePlot);
% camlight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test autorefine
% % close all
% R = 5;
% % nurbs = getDiskData('R',R,'parm',1,'t',0.3*R);
% % nurbs = getDiskData('R',R,'parm',2,'t',0.3*R);
% % nurbs = getDiskData('R',R,'parm',2);
% % nurbs = getEllipsoidData('C',R,'parm',1);
% % nurbs = getEllipsoidData('C',R,'parm',2);
% % nurbs = getEllipsoidData('C',R,'parm',2,'t',0.3*R);
% % nurbs = explodeNURBS(getEllipsoidData('C',R,'parm',1,'t',0.3*R));
% % nurbs = getBeTSSiM1Data('parm',1);
% % nurbs = getBeTSSiM1Data('parm',2);
% % nurbs = getBeTSSiM2Data('parm',1,'t',0.5);
% % nurbs = getBeTSSiM3Data('parm',2,'t',2);
% % nurbs = getBeTSSiM4Data('parm',1,'t',0.5);
% nurbs = getBeTSSiM4Data('parm',2,'t',0.5);
% % nurbs = getBeTSSiM5Data('parm',1);
% % nurbs = getBeTSSiM5Data('parm',2);
% % nurbs = getBCAData(2);
% % nurbs = getBC_modData2();
% % load('BeTSSi_BCA_p2.mat')
% % nurbs = getLshapeData('d',3);
% if 0 % Alter L-shape to test refinement algorithm
%     nurbs{2}.coeffs(3,2,:,2) = 1.5;
%     nurbs{3}.coeffs(1,1,2,:) = -1;
%     nurbs{4}.coeffs(2,:,2,2) = 2.5;
% end
% if 1 % Test Model4 (with parm==2) with knot insertions
%     xi = [0.1211111142,0.621111142];
%     newKnots = cell(1,16);
%     newKnots{1} = {[], xi, []};
%     newKnots{3} = {[], xi, []};
%     newKnots{5} = {xi, [], []};
%     newKnots{7} = {[], 1-xi, []};
%     newKnots{8} = {xi, [], []};
%     nurbs = insertKnotsInNURBS(nurbs,newKnots);
% end
% if 1
%     geometry = getTopology(nurbs);
%     refLength = R*pi/2;
%     % refLength = R*pi/2/50;
%     M = 7;
%     nurbs = autoRefineNURBS(nurbs,geometry.topology.connection,refLength/2^(M-1));
% end
% close all
% plotNURBS(nurbs,'plotParmDir',0,'plotControlPolygon',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test surfaceToVolume
% 
% nurbs_surf = getLshapeData('d',3);
% geometry = getTopology(nurbs_surf);
% t = 1;
% nurbs_vol = surfaceToVolume(nurbs_surf,geometry.topology.connection,t);
% close all
% plotNURBS(nurbs_surf,'plotParmDir',0,'plotControlPolygon',0);
% % plotNURBS(nurbs_vol,'plotParmDir',0,'plotControlPolygon',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Gordon Hall

% nurbs = getPrismData();
% % nurbs = getRectangleData();
% nurbs = elevateNURBSdegree(nurbs,1);
% nurbs{1}.coeffs = nurbs{1}.coeffs - 0.3*rand([size(nurbs{1}.coeffs,1),nurbs{1}.number]);
% nurbs = insertKnotsInNURBS(nurbs,[1,2,3]);
% nurbs = repeatNURBSknots(nurbs);
% 
% subnurbs = subNURBS(nurbs);
% nurbsGH = GordonHall(subnurbs);
% close all
% plotNURBS(explodeNURBS(nurbs),'plotParmDir',0,'plotControlPolygon',0,'plotObject',0);
% axis equal
% % plotNURBS(subnurbs,'plotParmDir',1,'plotControlPolygon',0);
% plotNURBS(explodeNURBS(nurbsGH),'plotParmDir',0,'plotControlPolygon',0,'plotObject',0);
% % plotNURBS(nurbs_vol,'plotParmDir',0,'plotControlPolygon',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find NURBS representation of part of sphere
% p = 4;
% nurbs = getOctSphereData('parm',2);
% nurbs = nurbs(1);
% plotNURBSoptions = struct('plotControlPolygon', 1, 'resolution',[100,100],'setViewAngle',[116,23]);
% nurbs = elevateNURBSdegree(nurbs,(p-4)*[1,1]);
% nurbs{1}.coeffs(4,:,:) = nurbs{1}.coeffs(4,:,:)/nurbs{1}.coeffs(4,end,end);
% plotNURBSoptions.colorFun = @(v) log10(abs(norm(v)-1));
% % plotNURBS(nurbs,plotNURBSoptions)
% % return
% nurbsArc = rotateNURBS(getArcData('Xi',[0,0,0,1,1,1],'theta',pi/4),'rotAxis',[0,0,1],'theta',pi/4);
% nurbsArc = elevateNURBSdegree(nurbsArc,p-2);
% nurbsArc(2) = rotateNURBS(nurbsArc,'rotAxis',[0,1,0],'theta',-pi/2);
% % plotNURBS(nurbsArc,plotNURBSoptions)
% 
% nurbsArc2 = rotateNURBS(getArcData('Xi',[0,0,0,1,1,1],'theta',pi/2-acos(1/sqrt(3))),'rotAxis',[0,0,1],'theta',pi/4);
% nurbsArc2 = flipNURBSparametrization(rotateNURBS(nurbsArc2,'rotAxis',[1,1,0],'theta',pi/2));
% nurbsArc2(2) = rotateNURBS(nurbsArc2,'rotAxis',[1,1,1],'theta',2*pi/3);
% nurbsArc2 = elevateNURBSdegree(nurbsArc2,p-2);
% % plotNURBS(nurbsArc2,plotNURBSoptions)
% nurbs{1}.coeffs(:,1,:) = nurbsArc2{2}.coeffs;
% nurbs{1}.coeffs(:,:,1) = nurbsArc2{1}.coeffs;
% nurbs{1}.coeffs(:,end,:) = nurbsArc{1}.coeffs;
% nurbs{1}.coeffs(:,:,end) = nurbsArc{2}.coeffs;
% 
% % plotNURBS(nurbs,plotNURBSoptions)
% % return
% xi = linspace(0,1,20);
% eta = xi;
% [XI,ETA] = meshgrid(xi,eta);
% XI = XI(:);
% ETA = ETA(:);
% allFree = 1;
% if allFree
%     x = nurbs{1}.coeffs(:);
% else
%     x = nurbs{1}.coeffs(:,2:end-1,2:end-1);
% end
% x = x(:);
% LB = zeros(size(x));
% UB = ones(size(x));
% UB(4:4:end) = 10;  
% % MaxIter = 1e4;
% MaxIter = 1e10;
% 
% options = optimset('Display','iter','TolX',eps,'TolFun',eps,'MaxFunEvals',MaxIter,'MaxIter',MaxIter);
% % x = [   0.511783921293191
% %    0.732174580986216
% %    0.507851118495701
% %    0.999953475604074
% %    0.576604895873008
% %    0.796718569515584
% %    0.332998545529378
% %    0.995278789296761
% %    0.566086367169920
% %    0.847830751471326
% %    0.145952998414980
% %    0.999768117724820
% %    0.325689942387339
% %    0.812808478090256
% %    0.570583027036215
% %    0.992393941771703
% %    0.396578985426332
% %    0.932723266024324
% %    0.372291974798927
% %    0.998459231423948
% %    0.398441669524650
% %    0.940455188274207
% %    0.134915783665834
% %    0.997781780019961
% %    0.126693188694058
% %    0.838329185251049
% %    0.576127300822254
% %    0.988156803122064
% %    0.170228254227707
% %    0.958374754025724
% %    0.389085706764339
% %    0.998085288701667
% %    0.212073720558690
% %    0.999947746649475
% %    0.173933595031326
% %    0.990111221560513];
% x = fminsearchbnd(@(x) objFun(nurbs,x,XI,ETA,allFree),x,LB,UB,options);
% 
% nurbs{1}.coeffs(:,2:end-1,2:end-1) = reshape(x,size(nurbs{1}.coeffs(:,2:end-1,2:end-1)));
% plotNURBS(nurbs,plotNURBSoptions)
% axis on
% xlabel x
% 
% function f = objFun(nurbs,x,XI,ETA,allFree)
% if allFree
%     nurbs{1}.coeffs(:,:,:) = reshape(x,size(nurbs{1}.coeffs));
% else
%     nurbs{1}.coeffs(:,2:end-1,2:end-1) = reshape(x,size(nurbs{1}.coeffs(:,2:end-1,2:end-1)));
% end
% X = evaluateNURBS(nurbs{1},[XI,ETA]);
% f = sum(abs(norm2(X)-1));
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot mesh from g2-files
% close all
% % plotType = '3D';
% plotType = 'surf';
% % plotType = 'algorithm';
% % plotType = 'algorithm2';
% % plotType = 'crossSection';
% % plotType = 'none';
% 
% % model = 'HuntHill/HuntHill';
% model = 'Vigra/Vigra50m';
% model = 'Vigra/VigraFree10m';
% model = 'Vigra/test';
% model = 'Sula/Sula50m';
% % model = 'test';
% % nurbs = read_g2('/home/zetison/OneDrive/work/openfoam/cylinder/cylinder.g2');
% % nurbs = read_g2('/home/zetison/temp/Re100/1/chorin/cyl2D.g2');
% % nurbs = read_g2('/home/zetison/OpenFOAM/OpenFOAM-v1912/run/cylinder/cylinder.g2');
% % nurbs = read_g2('/home/zetison/results/IFEM/cylinder/MESH7/1/chorin/cyl2D.g2');
% % resolution = 400*[1,1,1];
% % close all
% % nurbs1 = read_g2('/home/zetison/results/simra/out.g2');
% % nurbs2 = read_g2('/home/zetison/results/simra/out2.g2');
% % nurbs3 = read_g2('/home/zetison/results/simra/out3.g2');
% % nurbs4 = read_g2('/home/zetison/results/simra/out4.g2');
% % nurbsVol = read_g2('/home/zetison/results/simra/out3D.g2');
% % [h,maxC,minC] = plotNURBS(nurbsVol,'plotControlPolygon',false, 'plotObject', true, 'plotElementEdges', false, 'resolution',resolution, ...
% %     'colorFun', @(v) v(3), 'alphaValue',0.8);
% % 
% % plotNURBS(nurbsVol,'plotControlPolygon',false, 'plotObject', false, 'plotElementEdges', true, 'resolution',resolution);
% % 
% % [h,maxC,minC] = plotNURBS(nurbs,'plotControlPolygon',true, 'plotObject', false, 'resolution',resolution,'colorFun', @(v) v(3));
% % plotNURBS(nurbs,'plotControlPolygon',true, 'plotObject', false, 'plotElementEdges', true, 'resolution',resolution,'colorFun', @(v) v(3));
% % camlight
% % view(getView(1))
% % export_fig('/home/zetison/OneDrive/work/graphics/cylinder/M7', '-png', '-transparent', '-r200')
% % export_fig('/home/zetison/OneDrive/work/graphics/cylinder/M7_zoomed', '-png', '-transparent', '-r200')
% % nurbs = read_g2('/home/zetison/Downloads/TeslaValve/multipatch.g2');
% % plotNURBS(nurbs)
% % view(getView(1))
% % figure
% % nurbs = read_g2('/home/zetison/Downloads/TeslaValve/reparameterized.g2');
% % plotNURBS(nurbs)
% % view(getView(1))
% switch plotType
%     case 'crossSection'
%         resolution = 400*[1,1,1];
%         nurbs = read_g2(['/home/zetison/results/simra/' model '_3D.g2']);
%         nurbs = insertKnotsInNURBS(nurbs,{[], 0.5*ones(1,2), []});
%         nurbs = explodeNURBS(nurbs);
%         nurbs = subNURBS(nurbs(1),'at',[0,0;0,1;0,0]);
%         plotNURBS(nurbs,'plotControlPolygon',1, 'plotObject', 1, 'plotElementEdges', 1, 'resolution',resolution, 'color', getColor(1))
% %         export_fig('/home/zetison/results/simra/HuntHill/HuntHill', '-png', '-transparent', '-r200')
%         view(0,0)
%     case '3D'
%         resolution = 10*[1,1,1];
%         nurbs = read_g2(['/home/zetison/results/simra/' model '_3D.g2']);
%         plotNURBS(nurbs,'plotControlPolygon',0, 'plotObject', 1, 'plotElementEdges', 1, 'resolution',resolution, 'color', getColor(1))
%     case 'surf'
%         resolution = 8*[1,1,1];
%         nurbs = read_g2(['/home/zetison/results/simra/' model '.g2']);
%         [h,maxC,minC] = plotNURBS(nurbs,'plotControlPolygon',1, 'plotObject', 0, 'plotElementEdges', 0, 'resolution',resolution, ...
%             'colorFun', @(v) v(:,3), 'alphaValue',1);
% %         nurbs = read_g2(['/home/zetison/results/simra/' model '_3D.g2']);
% %         [h,maxC,minC] = plotNURBS(nurbs,'plotControlPolygon',1, 'plotObject', 0, 'plotElementEdges', 0, 'resolution',resolution,'alphaValue',1,'markerEdgeColor',[0,0,1]);
%         camlight
%         camproj('perspective')
%         colormap(getColorMap('MPL_terrain'))
%     case 'algorithm'
%         resolution = [3992,2726,1];
% %         resolution = [399,272,1];
% %         resolution = [40,27,1];
%         nurbs = read_g2(['/home/zetison/results/simra/' model '_3D.g2']);
%         P_sw = nurbs{1}.coeffs(1:2,1,1);
%         P_se = nurbs{1}.coeffs(1:2,end,1);
%         P_nw = nurbs{1}.coeffs(1:2,1,end);
%         P_ne = nurbs{1}.coeffs(1:2,end,end);
%         P_sw = P_sw + (P_sw-P_ne)/10;
%         P_se = P_se + (P_se-P_nw)/10;
%         P_nw = P_nw + (P_nw-P_se)/10;
%         P_ne = P_ne + (P_ne-P_sw)/10;
%         nurbs = subNURBS(nurbs(1),'at',[0,0;0,0;1,0]);
%         nurbs = explodeNURBS(nurbs);
%         nurbs1 = nurbs(18);
%         nurbs2 = glueNURBS(nurbs([14,19]),1);
%         nurbs3 = glueNURBS(nurbs([8,9]),2);
%         nurbs4 = glueNURBS(nurbs([7,12,17]),1);
%         nurbs5 = glueNURBS(nurbs([22,23,24]),2);
%         nurbs6 = glueNURBS(nurbs([10,15,20,25]),1);
%         nurbs7 = glueNURBS(nurbs([2,3,4,5]),2);
%         nurbs8 = glueNURBS(nurbs([1,6,11,16,21]),1);
%         viewAngle = [-20,35];
%         figure('Color','white')
%         [h,maxC,minC] = plotNURBS(nurbs(13),'plotControlPolygon',0, 'plotObject', 1, 'plotElementEdges', 1, 'resolution',resolution, ...
%             'colorFun', @(v) v(:,3), 'alphaValue',1,'view',viewAngle);
%         caxis([0,900])
%         plot3(P_sw(1)*[1,1+eps],P_sw(2)*[1,1+eps],-300*[1,1+eps])
%         plot3(P_sw(1)*[1,1+eps],P_sw(2)*[1,1+eps],1000*[1,1+eps])
%         plot3(P_se(1)*[1,1+eps],P_se(2)*[1,1+eps],-300*[1,1+eps])
%         plot3(P_se(1)*[1,1+eps],P_se(2)*[1,1+eps],1000*[1,1+eps])
%         plot3(P_nw(1)*[1,1+eps],P_nw(2)*[1,1+eps],-300*[1,1+eps])
%         plot3(P_nw(1)*[1,1+eps],P_nw(2)*[1,1+eps],1000*[1,1+eps])
%         plot3(P_ne(1)*[1,1+eps],P_ne(2)*[1,1+eps],-300*[1,1+eps])
%         plot3(P_ne(1)*[1,1+eps],P_ne(2)*[1,1+eps],1000*[1,1+eps])
%         cmap = jet(8);
%         cmap([8,4]) = cmap([4,8]);
%         camproj('perspective')
%         colormap(getColorMap('MPL_terrain'))
%         colorbar off
%         
%         frame_h = get(handle(gcf),'JavaFrame');
%         set(frame_h,'Maximized',1);
%         camlight
%         cropOption = '-c300,500,200,700';
%         export_fig('/home/zetison/OneDrive/work/graphics/simramesh/nurbs0', '-png', '-r200', cropOption)
%         endVec = [100,resolution(2);
%                   resolution(1),100];
%                   
%         nurbsCol = {nurbs1,nurbs2,nurbs3,nurbs4,nurbs5,nurbs6,nurbs7,nurbs8};
%         for i = 1:8
%             [h,maxC,minC] = plotNURBS(nurbsCol{i},'plotControlPolygon',0, 'plotObject', 1, 'plotElementEdges', 1, 'resolution',endVec(mod(i-1,2)+1,:),'color',cmap(i,:),'view',viewAngle);
%             caxis([0,900])
%             frame_h = get(handle(gcf),'JavaFrame');
%             set(frame_h,'Maximized',1);
%             export_fig(['/home/zetison/OneDrive/work/graphics/simramesh/nurbs' num2str(i)], '-png', '-r200', cropOption)
%         end
%     case 'algorithm2'
%         resolution = [3992,2726,1];
% %         resolution = [399,272,1];
%         resolution = [40,27,1];
%         nurbs = read_g2(['/home/zetison/results/simra/' model '_3D.g2']);
%         nurbs = subNURBS(nurbs(1),'at',[0,0;0,0;1,0]);
%         nurbs = explodeNURBS(nurbs);
%         nurbs1 = nurbs(18);
%         nurbs2 = glueNURBS(nurbs([14,19]),1);
%         nurbs3 = glueNURBS(nurbs([8,9]),2);
%         nurbs4 = glueNURBS(nurbs([7,12,17]),1);
%         nurbs5 = glueNURBS(nurbs([22,23,24]),2);
%         nurbs6 = glueNURBS(nurbs([10,15,20,25]),1);
%         nurbs7 = glueNURBS(nurbs([2,3,4,5]),2);
%         nurbs8 = glueNURBS(nurbs([1,6,11,16,21]),1);
%         viewAngle = [-20,35];
%         figure('Color','white')
%         [h,maxC,minC] = plotNURBS(nurbs(13),'plotControlPolygon',1, 'plotObject', 1, 'plotElementEdges', 1, 'resolution',resolution, ...
%             'colorFun', @(v) v(:,3), 'alphaValue',1,'view',viewAngle);
%         caxis([0,900])
%         camproj('perspective')
%         colormap(getColorMap('MPL_terrain'))
%         colorbar off
%         P1 = nurbs{13}.coeffs(:,end-1,12);
%         P2 = nurbs{13}.coeffs(:,end,12);
%         
%         
%         pts2 = nurbs5{1}.coeffs(:,:,12);
%         nurbs1D1 = createNURBSobject(pts,[0,0,0]);
%         frame_h = get(handle(gcf),'JavaFrame');
%         set(frame_h,'Maximized',1);
%         camlight
% %         export_fig('/home/zetison/OneDrive/work/graphics/simramesh/controlPts', '-png', '-r200')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Normal mesh
% nurbs = getLineData('x',1);
% p = 3;
% nurbs = elevateNURBSdegree(nurbs,p-1);
% % nurbs = insertKnotsInNURBS(nurbs,{linspace2(0,1,20)});
% nurbs = insertKnotsInNURBS(nurbs,{0.25});
% aveknt(nurbs{1}.knots{1},p+1).'
% nurbs{1}.coeffs(1,:)'
% a = unique(nurbs{1}.knots{1}).';
% [nurbs{1}.coeffs(1,[1,(p+3)/2:end-(p+1)/2,end]).'- a]
% close all
% resolution = [20,20];
% X = zeros(4,3);
% X([2,4],1) = 3;
% X([3,4],2) = 1;
% nurbs = getQuadrilateralData('X',X);
% nurbs = elevateNURBSdegree(nurbs,2);
% nurbs = insertKnotsInNURBS(nurbs,[5,1]);
% nurbs{1}.coeffs(2,5,[1,2]) = nurbs{1}.coeffs(2,5,[1,2]) + 0.3;
% nurbs{1}.coeffs(1,5,2) = nurbs{1}.coeffs(1,5,2) - 3.3;
% plotNURBS(nurbs,'plotControlPolygon',true, 'plotObject', 1, 'plotElementEdges', true, 'resolution',resolution);
% view(0,90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot bjelkeprofil
% close all
% nurbs = read_g2('/home/zetison/OneDrive/SINTEF/RaPiD/bjelkeprofil.g2');
% plotNURBS(nurbs, 'resolution',[100,100],'plotControlPolygon',1, 'color', getColor(1))
% view(0,90)
% export_fig('/home/zetison/OneDrive/SINTEF/RaPiD/bjelkeprofil', '-r200')
% 
% figure
% nurbs = read_g2('/home/zetison/OneDrive/SINTEF/RaPiD/volPatch.g2');
% nurbs = rotateNURBS(nurbs,'rotAxis','Xaxis');
% plotNURBS(nurbs, 'resolution',[100,100,100],'plotControlPolygon',0, 'color', getColor(1))
% camproj('perspective')
% camlight
% export_fig('/home/zetison/OneDrive/SINTEF/RaPiD/volPatch', '-r200')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of Lagrange-Hermite polynomials
% useHP = 1;
% x = [0;0.5;1];
% if useHP
%     x = mp(x);
% end
% X = linspace(x(1),x(end),1000);
% for n_D = 1:5
%     figure(n_D)
%     Y = getInterpolatingHermite(x,X,n_D,true);
%     for i = 1:size(Y,1)
%         printResultsToFile2(['../plotData/LagrangeHermite/n_D' num2str(n_D) 'i' num2str(i)], X.', Y(i,:).');
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of legendre polynomials
% N = 3;
% x_arr = linspace(-1,1,500)';
% P_t = zeros(length(x_arr),N);
% for n = 0:N
% %     fid = fopen(['../plotData/fundamentalFunctions/legendre_n' num2str(n) '.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
%     P_t(:,n+1) = legendre(n,x_arr);
% end
% plot(x_arr, P_t)

%% Derivatives
% hold on
% N = 4;
% x_arr = linspace(-1,1,1000)';
% % P_t = zeros(N,length(x_arr));
% col = hsv(N+1);
% for n = 0:N
% %     fid = fopen(['../plotData/fundamentalFunctions/legendre_n' num2str(n) '.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
%     for m = 1:length(x_arr)
%         t = x_arr(m);
% %         P_t(n+1,m) = legendreDeriv(n,t);
% %         fprintf(fid,'%1.15f\t%1.15f\n',t,P_t(n+1,m));
%     end
% %     fclose(fid);
% %     plot(x_arr, legendreDeriv(n,x_arr),'color',col(n+1,:))
%     plot(x_arr, legendreDeriv2(n,x_arr),'--', 'color',col(n+1,:))
% end

% figure(2)
% N = 5;
% theta = linspace(0,pi,500);
% P_t = zeros(N,length(theta));
% for n = 0:N
% %     fid = fopen(['../plotData/fundamentalFunctions/legendre_n' num2str(n) '_cosTheta.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
%     for m = 1:length(theta)
%         t = cos(theta(m));
%         P_t(n+1,m) = legendre(n,t);
% %         fprintf(fid,'%1.15f\t%1.15f\n',theta(m),legendre3(n,t));
%     end
% %     fclose(fid);
% end
% plot(theta, P_t)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of Spherical Bessel functions
% close all
% figure
% no2Dpoints = 900;
% % x = linspace(0,100,no2Dpoints);
% x = linspace2(0,1/2,no2Dpoints);
% x = x(end);
% % x = 8.04/9;
% % x = 8.04;
% no2Dpoints = numel(x);
% noBesselFunctions = 100000;
% J = zeros(noBesselFunctions,no2Dpoints);
% j = zeros(noBesselFunctions,no2Dpoints);
% dj = zeros(noBesselFunctions,no2Dpoints);
% Y = zeros(noBesselFunctions,no2Dpoints);
% H = zeros(noBesselFunctions,no2Dpoints);
% Hr = zeros(noBesselFunctions,no2Dpoints);
% dHr = zeros(noBesselFunctions,no2Dpoints);
% % DJ = zeros(noBesselFunctions,no2Dpoints);
% N_arr = 0:noBesselFunctions-1;
% % N_arr = 1:noBesselFunctions;
% for n = N_arr
%     z = x*n;
% %     z = x;
%     J(n+1,:) = besselj(n,z);
%     Y(n+1,:) = bessely(n,z);
%     j1 = sqrt(pi/2)./sqrt(z).*besselj(n+0.5,z);
%     j2 = sqrt(pi/2)./sqrt(z).*besselj(n+0.5+1,z);
%     y1 = sqrt(pi/2)./sqrt(z).*bessely(n+0.5,z);
%     y2 = sqrt(pi/2)./sqrt(z).*bessely(n+0.5+1,z);
%     h1 = j1 + 1i*y1;
%     h2 = j2 + 1i*y2;
%     j(n+1,:) = abs(j1);
%     dj_temp = n./z.*j1 - j2;
%     dy_temp = n./z.*y1 - y2;
%     dj(n+1,:) = abs(dj_temp);
%     H(n+1,:) = abs(h1);
%     Hr(n+1,:) = abs(h2./h1).*z/n;
% %     dHr(n+1,:) = abs(h1)./abs(n./z.*h1 - h2);
%     indices = n > z;
%     dHr(n+1,indices) = abs(h1(indices))./abs(n./z(indices).*h1(indices) - h2(indices));
% %     dHr(n+1,:) = abs(h1)./abs(dj_temp + 1i*dy_temp);
% end
% % plot(x,J.')
% % plot(N_arr,Y.')
% % semilogy(N_arr,-Y.','color','red')
% % hold on
% % semilogy(N_arr,J.','color','blue')
% % semilogy(N_arr,Hr.','color','blue')
% semilogy(N_arr,dHr.')
% % semilogy(ones(2,1).*ceil(x(1)+eps),ylim.'.*ones(1,1),'--','color','black')
% % figure
% % plot(x,dHr(end,:))
% % hold on
% % plot(x,-log(1-x))
% % % plot(x,-0.911665108328680*log(0.995188245896694 -1.005057731715838*x))
% % plot(x,1./(1-x).^0.5)
% % if true
% % %     semilogy(N_arr,j.','color','red')
% % %     hold on
% % %     semilogy(N_arr,dj.','color','blue')
% % %     semilogy(x,j,'color','red')
% % %     hold on
% % %     semilogy(x,dj,'color','blue')
% %     semilogy(x,H)
% % else
% %     semilogy(N_arr(1:end-1),j(1:end-1,:).','color','red')
% %     hold on
% %     semilogy(N_arr(2:end),dj(2:end,:).','color','blue')
% % end
% % ylim([-1.5 1.1])
% % legend('j_0','j_1','Dj_0','Dj_1','Location','Best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of Bessel functions of first kind
% x = linspace(0,20,no2Dpoints);
% noBesselFunctions = 5;
% J = zeros(noBesselFunctions,no2Dpoints);
% % DJ = zeros(noBesselFunctions,no2Dpoints);
% for n = 0:noBesselFunctions-1
%     J(n+1,:) = besselj(n,x);
%     DJ(n+1,:) = besseljDeriv(n,x);
% %     fid = fopen(['../plotData/fundamentalFunctions/bessel1_n' num2str(n) '.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
% %     for m = 1:length(x)
% %         fprintf(fid,'%1.15f\t%1.15f\n',x(m),J(n+1,m));
% %     end
% %     fclose(fid);
% end
% plot(x,J)
% ylim([-1.5 1.1])
% legend('j_0','j_1','Dj_0','Dj_1','Location','Best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of Bessel functions of second kind
% x = linspace(0.001,20,no2Dpoints);
% noBesselFunctions = 5;
% Y = zeros(noBesselFunctions,no2Dpoints);
% DY = zeros(noBesselFunctions,no2Dpoints);
% for n = 0:noBesselFunctions-1
%     Y(n+1,:) = bessely(n,x);
%     DY(n+1,:) = besselyDeriv(n,x);
%     fid = fopen(['../plotData/fundamentalFunctions/bessel2_n' num2str(n) '.txt'],'wt+','b');
%     fprintf(fid,'x\t\t\ty\n');
%     for m = 1:length(x)
%         if -2 < Y(n+1,m) && Y(n+1,m) < 1
%             fprintf(fid,'%1.15f\t%1.15f\n',x(m),Y(n+1,m));
%         end
%     end
%     fclose(fid);
% end
% plot(x,Y,x,DY,'--')
% ylim([-1 0.6])
% % ylim([-5 1])
% % xlim([0 5])
% 
% legend('y_0','y_1','Dy_0','Dy_1','Location','Best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of scaled Bessel functions
% x = linspace(0.01,200,no2Dpoints);
% noBesselFunctions = 5;
% sJ = zeros(noBesselFunctions,no2Dpoints);
% sY = zeros(noBesselFunctions,no2Dpoints);
% for n = 0:noBesselFunctions-1
%     sJ(n+1,:) = besselj(n,x)./x.^n;
% %     sY(n+1,:) = bessely(n,x).*x.^(n-1);
% end
% % plot(x,sJ,x,sY)
% plot(x,sJ)
% % ylim([-1 0.6])
% % ylim([-5 1])
% % xlim([0 5])
% 
% % legend('y_0','y_1','Dy_0','Dy_1','Location','Best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of spherical Bessel functions of first kind
% x = linspace(0,20,no2Dpoints);
% noBesselFunctions = 5;
% j = zeros(noBesselFunctions,no2Dpoints);
% Dj = zeros(noBesselFunctions,no2Dpoints);
% for n = 0:noBesselFunctions-1
%     j(n+1,:) = bessel1(n,x);
%     Dj(n+1,:) = bessel1Deriv(n,x);
%     fid = fopen(['../plotData/fundamentalFunctions/sphericalBessel1_n' num2str(n) '.txt'],'wt+','b');
%     fprintf(fid,'x\t\t\ty\n');
%     for m = 1:length(x)
%         fprintf(fid,'%1.15f\t%1.15f\n',x(m),j(n+1,m));
%     end
%     fclose(fid);
% end
% plot(x,j,x,Dj,'--')
% ylim([-1.5 1.1])
% legend('j_0','j_1','Dj_0','Dj_1','Location','Best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of spherical Bessel functions of second kind
% x = linspace(0.1,20,no2Dpoints);
% noBesselFunctions = 5;
% y = zeros(noBesselFunctions,no2Dpoints);
% Dy = zeros(noBesselFunctions,no2Dpoints);
% for n = 0:noBesselFunctions-1
%     y(n+1,:) = bessel2(n,x);
%     Dy(n+1,:) = bessel2Deriv(n,x);
%     fid = fopen(['../plotData/fundamentalFunctions/sphericalBessel2_n' num2str(n) '.txt'],'wt+','b');
%     fprintf(fid,'x\t\t\ty\n');
%     for m = 1:length(x)
%         if -2 < y(n+1,m) && y(n+1,m) < 1
%             fprintf(fid,'%1.15f\t%1.15f\n',x(m),y(n+1,m));
%         end
%     end
%     fclose(fid);
% end
% plot(x,y,x,Dy,'--')
% ylim([-1.5 0.4])
% 
% legend('y_0','y_1','Dy_0','Dy_1','Location','Best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of sperical Hankel functions
% N = 5;
% M = 40;
% cc = [1 0 0;
%       1 1 154/255;
%       0 1 0;
%       0 0 1;
%       1 0 1;
%       0 0 0];
% k_arr = linspace(0,5,100);
% H_arr = zeros(size(k_arr));
% dH_arr = zeros(size(k_arr));
% for n = 0:N
%     for j = 1:length(k_arr)
%         k = k_arr(j);
%         H_arr(j) = sphericalHankel1(n,k,M);
%         dH_arr(j) = sphericalHankel1Deriv(n,k,M);
%     end
%     if true % Plot the real part of Hankel function (Bessel functions of first kind)
%         plot(k_arr, real(H_arr),'color',cc(n+1,:))
%         xlim([0 5])
%         hold on
%         plot(k_arr, real(dH_arr),'color',cc(n+1,:),'Linestyle','--')
%     else % Plot the real part of Hankel function (Bessel functions of second kind)
%         plot(k_arr, imag(H_arr),'color',cc(n+1,:))
%         ylim([-5 1])
%         xlim([1 5])
%         hold on
%         plot(k_arr, imag(dH_arr),'color',cc(n+1,:),'Linestyle','--')
%     end
%     %legend('Hankel function', 'Derivative of Hankel function')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot B-spline basis functions
% Xi = [0 0 0 1 2 2 2];
% Xi = Xi/Xi(end);
% p = 2;
% n = length(Xi)-(p+1);
% figure(1)
% set(0, 'defaultTextInterpreter', 'latex'); 
% xlabel('$\xi$')
% % Plot basis functions of second order
% for i = 1:n
% %     fid = fopen(['../plotData/Bsplines/order2number' num2str(i) '.txt'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% %     fprintf(fid,'xi\t\t\tBspline\n');
%     h = plotBspline(i,p,n,Xi,10000);
%     hold on
%     xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i))
%     plot([xi xi],[0,1],'--','color','black')
% %     xi = 0.32
% %     2/(Xi(i+2)-Xi(i))*(Xi(i+2)-xi)/(Xi(i+2)-Xi(i+1))-2/(Xi(i+3)-Xi(i+1))*(xi-Xi(i+1))/(Xi(i+2)-Xi(i+1))
%     
% %     xData = get(h,'XData');
% %     yData = get(h,'YData');
% %     for j = 1:length(xData)
% %         fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% %     end
% %     fclose(fid);
% end
% 
% % Plot basis functions of first order
% p = p-1;
% n = n+1;
% for i = 1:n
%     if Xi(i) == Xi(i+p+1)
%         continue
%     end
%     fid = fopen(['../plotData/Bsplines/order1number' num2str(i) '.txt'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%     fprintf(fid,'xi\t\t\tBspline\n');
%     h = plotBspline(i,p,n,Xi,no2Dpoints);
%     xData = get(h,'XData');
%     yData = get(h,'YData');
%     for j = 1:length(xData)
%         fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
%     end
%     fclose(fid);
% end
% 
% % Plot basis functions of zeroth order
% p = p-1;
% n = n+1;
% for i = 1:n
%     if Xi(i) == Xi(i+p+1)
%         continue
%     end
%     fid = fopen(['../plotData/Bsplines/order0number' num2str(i) '.txt'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%     fprintf(fid,'xi\t\t\tBspline\n');
%     h = plotBspline(i,p,n,Xi,no2Dpoints);
%     xData = get(h,'XData');
%     yData = get(h,'YData');
%     for j = 1:length(xData)
%         fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
%     end
%     fclose(fid);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot B-spline basis functions and derivatives
% p = 3;
% % Xi = [zeros(1,p+1) 0.1 0.25*ones(1,p)  0.5*ones(1,p) 0.75*ones(1,p) ones(1,p+1)]; % Compute a uniform knot vector
% Xi = [zeros(1,p+1) 0.5 ones(1,p+1)]; % Compute a uniform knot vector
% n = length(Xi)-(p+1);
% figure(1)
% set(0, 'defaultTextInterpreter', 'latex'); 
% xlabel('$\xi$')
% % Plot basis functions of second order
% for i = 2
% %     fid = fopen(['../plotData/Bsplines/order2number' num2str(i) '.txt'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% %     fprintf(fid,'xi\t\t\tBspline\n');
%     h = plotBspline(i,p,n,Xi,1000,1);
%     hold on
% %     xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i))
% %     plot([xi xi],[0,1])
% %     xi = 0.32
% %     2/(Xi(i+2)-Xi(i))*(Xi(i+2)-xi)/(Xi(i+2)-Xi(i+1))-2/(Xi(i+3)-Xi(i+1))*(xi-Xi(i+1))/(Xi(i+2)-Xi(i+1))
%     
% %     xData = get(h,'XData');
% %     yData = get(h,'YData');
% %     for j = 1:length(xData)
% %         fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% %     end
% %     fclose(fid);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot B-spline basis functions showing k_refinement
% noPts = 200;
% for runC0 = 0
%     for p = 2:5
%         close all
%         % Xi = [zeros(1,p+1) 0.1 0.25*ones(1,p)  0.5*ones(1,p) 0.75*ones(1,p) ones(1,p+1)]; % Compute a uniform knot vector
%         if ~runC0
%             Xi = [zeros(1,p+1) 1 2 3  4*ones(1,p+1)]/4; % Compute a uniform knot vector
%         else
%             Xi = [zeros(1,p+1) ones(1,p) 2*ones(1,p) 3*ones(1,p) 4*ones(1,p+1)]/4; % Compute a uniform knot vector
%         end
%         n = length(Xi)-(p+1);
%         figure(1)
%         set(0, 'defaultTextInterpreter', 'latex'); half_shel
%         xlabel('$\xi$')
%         % Plot basis functions of second order
%         for i = 1:n
% %             if i == p+1 && p == 5 && runC0
% %                 keyboard
% %             end
%             if ~runC0
%                 fid = fopen(['../plotData/Bsplines2/Cpm1p' num2str(p) 'i' num2str(i) '.dat'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%             else
%                 fid = fopen(['../plotData/Bsplines2/C0p' num2str(p) 'i' num2str(i) '.dat'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%             end
%             fprintf(fid,'x\t\t\ty\n');
%             h = plotBspline(i,p,n,Xi,noPts);
%             hold on
%         %     xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i))
%         %     plot([xi xi],[0,1])
%         %     xi = 0.32
%         %     2/(Xi(i+2)-Xi(i))*(Xi(i+2)-xi)/(Xi(i+2)-Xi(i+1))-2/(Xi(i+3)-Xi(i+1))*(xi-Xi(i+1))/(Xi(i+2)-Xi(i+1))
% 
%             xData = get(h,'XData');
%             yData = get(h,'YData');
%             for j = 1:length(xData)
%                 fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
%             end
%             fclose(fid);
%         end
%     end
% end
% p = 3;
% for m = 2:p+1
%     close all
%     Xi = [zeros(1,p+1) 1 2 3*ones(1,m) 4*ones(1,p+1)]/4; % Compute a uniform knot vector
%     n = length(Xi)-(p+1);
%     figure(1)
%     set(0, 'defaultTextInterpreter', 'latex'); 
%     xlabel('$\xi$')
%     for i = 1:n
%         fid = fopen(['../plotData/Bsplines2/Cpm' num2str(m) 'p' num2str(p) 'i' num2str(i) '.dat'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%         fprintf(fid,'x\t\t\ty\n');
%         h = plotBspline(i,p,n,Xi,noPts);
%         hold on
%         xData = get(h,'XData');
%         yData = get(h,'YData');
%         for j = 1:length(xData)
%             if j~=1 && xData(j) == xData(j-1)
%                 fprintf(fid,'\n');
%             end   
%             fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));             
%         end
%         fclose(fid);
%     end
% end
% 
% figure(42)
% for p = 1:5
%     counter = 1;
%     close all
%     for j = 1:4
%         for i = 1:p+1
%             xi = linspace((j-1)/4,j/4,noPts).';
%             Xi = linspace(xi(1),xi(end),p+1);
%             y = lagrangePolynomials(xi,i,p+1,Xi);
%             if i == p+1 && j ~= 4
%                 y = [lagrangePolynomials(xi,p+1,p+1,Xi); lagrangePolynomials(xi(2:end),1,p+1,Xi)];
%                 xiTemp = linspace(j/4,(j+1)/4,noPts).';
%                 xi = [xi; xiTemp(2:end)];
%             end
%             if j > 1
%                 xi = [0; xi];
%                 y = [0; y];
%             end
%             if j < 4
%                 xi = [xi; 1];
%                 y = [y; 0];
%             end
%             if j > 1 && i == 1
%                 continue
%             end
%             printResultsToFile2(['plotData/Bsplines2/FEMp' num2str(p) 'i' num2str(counter)], xi, y);
%             hold on
%             plot(xi,y)
%             counter = counter + 1;
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot B-spline basis functions showing k_refinement
% noPts = 200;
% for runC0 = [0,1]
%     for p = 1:5
%         close all
%         % Xi = [zeros(1,p+1) 0.1 0.25*ones(1,p)  0.5*ones(1,p) 0.75*ones(1,p) ones(1,p+1)]; % Compute a uniform knot vector
%         if ~runC0
%             Xi = [zeros(1,p+1) 1 2 3  4*ones(1,p+1)]/4; % Compute a uniform knot vector
%         else
%             Xi = [zeros(1,p+1) ones(1,p) 2*ones(1,p) 3*ones(1,p) 4*ones(1,p+1)]/4; % Compute a uniform knot vector
%         end
%         n = length(Xi)-(p+1);
%         figure(1)
%         set(0, 'defaultTextInterpreter', 'latex');
%         xlabel('$\xi$')
%         % Plot basis functions of second order
%         for i = 1:n
%             if ~runC0
%                 fid = fopen(['../plotData/Bsplines3/Cpm1p' num2str(p) 'i' num2str(i) '.dat'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%             else
%                 fid = fopen(['../plotData/Bsplines3/C0p' num2str(p) 'i' num2str(i) '.dat'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%             end
%             fprintf(fid,'x\t\t\ty\n');
%             h = plotBspline(i,p,n,Xi,noPts);
%             hold on
%         %     xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i))
%         %     plot([xi xi],[0,1])
%         %     xi = 0.32
%         %     2/(Xi(i+2)-Xi(i))*(Xi(i+2)-xi)/(Xi(i+2)-Xi(i+1))-2/(Xi(i+3)-Xi(i+1))*(xi-Xi(i+1))/(Xi(i+2)-Xi(i+1))
% 
%             xData = get(h,'XData');
%             yData = get(h,'YData');
%             for j = 1:length(xData)
%                 fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
%             end
%             fclose(fid);
%         end
%     end
% end
% p = 1;
% close all
% Xi = [0 0 1 2 3 4 5 6 7 8 8]/8; % Compute a uniform knot vector
% n = length(Xi)-(p+1);
% figure(1)
% set(0, 'defaultTextInterpreter', 'latex');
% xlabel('$\xi$')
% % Plot basis functions of second order
% for i = 1:n
%     fid = fopen(['../plotData/Bsplines3/C0_h_p' num2str(p) 'i' num2str(i) '.dat'],'wt+','b');
%     fprintf(fid,'x\t\t\ty\n');
%     h = plotBspline(i,p,n,Xi,noPts);
%     hold on
% %     xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i))
% %     plot([xi xi],[0,1])
% %     xi = 0.32
% %     2/(Xi(i+2)-Xi(i))*(Xi(i+2)-xi)/(Xi(i+2)-Xi(i+1))-2/(Xi(i+3)-Xi(i+1))*(xi-Xi(i+1))/(Xi(i+2)-Xi(i+1))
% 
%     xData = get(h,'XData');
%     yData = get(h,'YData');
%     for j = 1:length(xData)
%         fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
%     end
%     fclose(fid);
% end
% p = 5;
% for k = p-1:-1:-1
%     close all
%     Xi = [zeros(1,p+1) 1 2 3*ones(1,p-k)  4*ones(1,p+1)]/4; % Compute a uniform knot vector
%     n = length(Xi)-(p+1);
%     figure(1)
%     set(0, 'defaultTextInterpreter', 'latex');
%     xlabel('$\xi$')
%     % Plot basis functions of second order
%     for i = 1:n
%         fid = fopen(['../plotData/Bsplines3/Ck' num2str(k) 'p' num2str(p) 'i' num2str(i) '.dat'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
%         fprintf(fid,'x\t\t\ty\n');
%         h = plotBspline(i,p,n,Xi,noPts);
%         hold on
%     %     xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i))
%     %     plot([xi xi],[0,1])
%     %     xi = 0.32
%     %     2/(Xi(i+2)-Xi(i))*(Xi(i+2)-xi)/(Xi(i+2)-Xi(i+1))-2/(Xi(i+3)-Xi(i+1))*(xi-Xi(i+1))/(Xi(i+2)-Xi(i+1))
% 
%         xData = get(h,'XData');
%         yData = get(h,'YData');
%         for j = 1:length(xData)
%             if isnan(yData(j))
%                 fprintf(fid,'\n',xData(j),yData(j));
%             else
%                 fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
%             end
%         end
%         fclose(fid);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot a spline curve
% P = [5 3 1 4 7 7;
%      2 2 6 3 6 1];
% Xi = [0 0 0 1 2 2 3 3 3];
% Xi = Xi/Xi(end);
% 
% % Plot spline curve
% splineObj = createSplineObject(P,Xi);
% 
% fid = fopen('plotData/splineCurve/degreeElevationSplineCurve.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% h = plotSpline(splineObj, no2Dpoints);
% xData = get(h,'XData');
% yData = get(h,'YData');
% for j = 1:length(xData)
%     fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% end
% fclose(fid);
% 
% % Plot control polygon
% fid = fopen('plotData/splineCurve/controlPolygon.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(P)
%     fprintf(fid,'%1.15f\t%1.15f\n',P(1,j),P(2,j));
% end
% fclose(fid);
% 
% % Plot knot locations
% elements = unique(Xi);
% fid = fopen('plotData/splineCurve/knotLocations.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(elements)
%     xi = elements(j);
%     E = evaluateSpline(splineObj, xi);
%     fprintf(fid,'%1.15f\t%1.15f\n',E(1),E(2));
% end
% fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot a spline volume for thesis (control points)
% t = 0.5;
% P = [5 2; 3 2; 1 6; 4 3; 7 6; 7 1];
% Q = [5 2-t; 3-t 2-t; 1-2*t+0.5 6+3*t-0.8; 4 3+1.2*t; 7+t-0.2 7+t+1.2*t-1-0.7; 7+t 1];
% controlPts = zeros(3,6,2,2);
% Xi = [0 0 0 1 2 2 3 3 3];
% Eta = [0 0 1 1];
% Zeta = [0 0 1 1];
% y = 0.5;
% z = 5;
% 
% controlPts(:,1,1,1) = [  P(1,1)   P(1,2)   0];
% controlPts(:,2,1,1) = [  P(2,1)   P(2,2)   0];
% controlPts(:,3,1,1) = [  P(3,1)   P(3,2)   0];
% controlPts(:,4,1,1) = [  P(4,1)   P(4,2)   0];
% controlPts(:,5,1,1) = [  P(5,1)   P(5,2)   0];
% controlPts(:,6,1,1) = [  P(6,1)   P(6,2)   0];
% 
% controlPts(:,1,2,1) = [  Q(1,1)   Q(1,2)   0];
% controlPts(:,2,2,1) = [  Q(2,1)   Q(2,2)   0];
% controlPts(:,3,2,1) = [  Q(3,1)   Q(3,2)   0];
% controlPts(:,4,2,1) = [  Q(4,1)   Q(4,2)   0];
% controlPts(:,5,2,1) = [  Q(5,1)   Q(5,2)   0];
% controlPts(:,6,2,1) = [  Q(6,1)   Q(6,2)   0];
% 
% controlPts(:,1,1,2) = [  P(1,1)   P(1,2)   z];
% controlPts(:,2,1,2) = [  P(2,1)   P(2,2)   z];
% controlPts(:,3,1,2) = [  P(3,1)   P(3,2)   z];
% controlPts(:,4,1,2) = [  P(4,1)   P(4,2)   z];
% controlPts(:,5,1,2) = [  P(5,1)   P(5,2)   z];
% controlPts(:,6,1,2) = [  P(6,1)   P(6,2)   z];
% 
% controlPts(:,1,2,2) = [  Q(1,1)   Q(1,2)   z];
% controlPts(:,2,2,2) = [  Q(2,1)   Q(2,2)   z];
% controlPts(:,3,2,2) = [  Q(3,1)   Q(3,2)   z];
% controlPts(:,4,2,2) = [  Q(4,1)   Q(4,2)   z];
% controlPts(:,5,2,2) = [  Q(5,1)   Q(5,2)   z];
% controlPts(:,6,2,2) = [  Q(6,1)   Q(6,2)   z];
% 
% 
% 
% splineObj = createSplineObject(controlPts,{Xi, Eta, Zeta});
% plotSpline(splineObj, [1000 2 2]);
% grid off
% axis off
% view(-7,80)
% 
% % export_fig('../graphics/splineVolume', '-png', '-transparent', '-r200')
% 
% hold off
% splineObj = insertKnotsInSplineObjects(splineObj,{[0.5 1.5]/3 [0.5] [0.25 0.5 0.75]});
% 
% plotSpline(splineObj, [1000 2 2]);
% grid off
% axis off
% view(-7,80)
% 
% % export_fig('../graphics/splineVolumeKnotIns4', '-png', '-transparent', '-r200')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Knot insertion
% P = [5 3 1 4 7 7;
%      2 2 6 3 6 1];
% Xi = [0 0 0 1 2 2 3 3 3];
% 
% newKnots = [0.5 1.5]/Xi(end);
% Xi = Xi/Xi(end);
% p = 2;
% n = length(Xi) - (p+1);
% [P, Xi] = insertKnotsInBsplines(n, p, Xi, newKnots, P);
% 
% splineObj = createSplineObject(P,Xi);

% % Plot spline curve
% fid = fopen('plotData/splineCurve/knotInsertionSplineCurve.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% h = plotSpline(splineObj, no2Dpoints);
% xData = get(h,'XData');
% yData = get(h,'YData');
% for j = 1:length(xData)
%     fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% end
% fclose(fid);
% 
% % Plot control polygon
% fid = fopen('plotData/splineCurve/knotInsertionControlPolygon.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(P)
%     fprintf(fid,'%1.15f\t%1.15f\n',P(1,j),P(2,j));
% end
% fclose(fid);
% 
% % Plot knot locations
% elements = unique(Xi);
% fid = fopen('plotData/splineCurve/knotInsertionKnotLocations.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(elements)
%     xi = elements(j);
%     E = evaluateSpline(splineObj, xi);
%     fprintf(fid,'%1.15f\t%1.15f\n',E(1),E(2));
% end
% fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Degree elevation
% m = 1;
% P = [5 3 1 4 7 7;
%      2 2 6 3 6 1];
% Xi = [0 0 0 1 2 2 3 3 3];
% 
% newKnots = [0.5 1.5]/Xi(end);
% Xi = Xi/Xi(end);
% p = 2;
% n = length(Xi) - (p+1);
% [P, Xi] = elevateBsplinesDegree(n, p, Xi, P, m);
% 
% splineObj = createSplineObject(P,Xi);
% 
% % Plot spline curve
% fid = fopen('plotData/splineCurve/degreeElevationSplineCurve.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% h = plotSpline(splineObj, no2Dpoints);
% xData = get(h,'XData');
% yData = get(h,'YData');
% for j = 1:length(xData)
%     fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% end
% fclose(fid);
% 
% % Plot control polygon
% fid = fopen('plotData/splineCurve/degreeElevationControlPolygon.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(P)
%     fprintf(fid,'%1.15f\t%1.15f\n',P(1,j),P(2,j));
% end
% fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Circular arc
% nurbs = getArcData;
% 
% % Plot NURBS curve
% % fid = fopen('../plotData/NURBScircle/NURBScircle.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% % fprintf(fid,'x\t\t\ty\n');
% nurbs = elevateDegreeInPatches(nurbs,4);
% h = plotNURBS(nurbs{1}, {'resolution',no2Dpoints});
% hold on
% plotControlPts(nurbs)
% axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of NURBS curve (circle)
% 
% nurbs = getCircleData;
% if true
%     nurbs.knots{1} = [-1/4*ones(1,3), nurbs.knots{1}(2:end-1), 5/4*ones(1,3)];
%     nurbs.coeffs = [repmat(nurbs.coeffs(:,1),1,2), nurbs.coeffs, repmat(nurbs.coeffs(:,1),1,2)];
%     nurbs.number = nurbs.number+4;
% %     coef
% end
% % Plot NURBS curve
% % fid = fopen('../plotData/NURBScircle/NURBScircle.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% % fprintf(fid,'x\t\t\ty\n');
% nurbs = elevateDegreeInPatches(nurbs,1);
% h = plotNURBS(nurbs{1}, {'resolution',no2Dpoints});
% hold on
% plotControlPts(nurbs)
% axis equal

% xData = get(h,'XData');
% yData = get(h,'YData');
% for j = 1:length(xData)
%     fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% end
% fclose(fid);

% % Plot control polygon
% P = nurbs.coeffs(1:2,:);
% fid = fopen('../plotData/NURBScircle/controlPolygon.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(P)
%     fprintf(fid,'%1.15f\t%1.15f\n',P(1,j),P(2,j));
% end
% fclose(fid);

% Plot NURBS curve after knot insertion
% nurbs = insertKnotsInNURBS(nurbs,[0.5 0.5 1.5 2.25 2.5 2.75]/4);
% 
% fid = fopen('../plotData/NURBScircle/knotInsertionNURBScircle.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% h = plotNURBS(nurbs, no2Dpoints);
% xData = get(h,'XData');
% yData = get(h,'YData');
% for j = 1:length(xData)
%     fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% end
% fclose(fid);


% % Plot control polygon after knot insertion
% P = nurbs.coeffs(1:2,:);
% fid = fopen('../plotData/NURBScircle/knotInsertionControlPolygon.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(P)
%     fprintf(fid,'%1.15f\t%1.15f\n',P(1,j),P(2,j));
% end
% fclose(fid);
% 
% % Plot knot locations after knot insertion
% Xi = nurbs.knots;
% elements = unique(Xi);
% 
% fid = fopen('../plotData/NURBScircle/knotInsertionKnotLocations.txt','wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% fprintf(fid,'x\t\t\ty\n');
% for j = 1:length(elements)
%     xi = elements(j);
%     E = evaluateNURBS(nurbs, xi);
%     fprintf(fid,'%1.15f\t%1.15f\n',E(1),E(2));
% end
% fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Spline basis functions in complex plane
% close all
% Xi = [0 0 0 1 2 2 2];
% Xi = Xi/Xi(end);
% p = 2;
% n = length(Xi)-(p+1);
% figure(1)
% set(0, 'defaultTextInterpreter', 'latex'); 
% xlabel('$\xi$')
% % Plot basis functions of second order
% % xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i));
% for i = 2
% %     fid = fopen(['../plotData/Bsplines/order2number' num2str(i) '.txt'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
% %     fprintf(fid,'xi\t\t\tBspline\n');
%     figure(i)
%     hold on
%     h = plotBspline(i,p,n,Xi,100,false,1);
% %     xi = linspace(0,1,1000);
% %     xi = xi + 1i*xi.';
% %     i1 = findKnotSpan(n, p, xi, Xi);
% %     N = BsplineBasis(i1, xi, p, Xi, 0);
% %     h = plotBspline(p,n,Xi,100);
% %     hold on
% %     xi = (Xi(i+2)*(Xi(i+3)-Xi(i+1))+Xi(i+1)*(Xi(i+2)-Xi(i)))/(Xi(i+3)-Xi(i+1)+Xi(i+2)-Xi(i))
% %     surf()
% %     xi = 0.32
% %     2/(Xi(i+2)-Xi(i))*(Xi(i+2)-xi)/(Xi(i+2)-Xi(i+1))-2/(Xi(i+3)-Xi(i+1))*(xi-Xi(i+1))/(Xi(i+2)-Xi(i+1))
%     
% %     xData = get(h,'XData');
% %     yData = get(h,'YData');
% %     for j = 1:length(xData)
% %         fprintf(fid,'%1.15f\t%1.15f\n',xData(j),yData(j));
% %     end
% %     fclose(fid);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solid cylinder
% L = 5;  % Height
% R = 1;    % Inner surface radius
% R_o = 2;    % Outer surface radius
% 
% solid = getSolidCylinderData(R, R_o, L);
% 
% plotNURBS(solid,[100 0 0])
% grid off
% axis off
% view(50,30)
% drawnow
% export_fig('../graphics/solidCylinder', '-png', '-transparent', '-r200')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pinched hemisphere
% for M = 6 %[1 2 3 6]
%     t = 0.04;  % Thickness
%     R = 10; % Mid surface radius
%     R = R-t/2;
%     R_o = R+t/2;
%     solid = getPinchedHemisphereData(R, R_o);
%     solid = insertKnotsInNURBS(solid,{linspace2(0,1,2^(M-1)-1) linspace2(0,1,2^(M-1)-1) 0.5});
%     plotNURBS(solid,[100 100 0])
%     
%     axis equal
%     view(40-240,30)
%     camlight;
%     set(gca, 'Color', 'none');
%     grid off
%     axis off
% 
%     export_fig(['../graphics/pinchedHemisphereMesh' num2str(M)], '-png', '-transparent', '-r300')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pinched cylinder

% for M = 7 %[1 2 3 7]
%     t = 3;  % Thickness
%     R = 300; % Mid surface radius
%     R = R-t/2;
%     R_o = R+t/2;
%     L = 600;    % Length
%     solid = getPinchedCylinderData(R, R_o, L/2);
% 
%     solid = insertKnotsInNURBS(solid,{linspace2(0,1,2^(M-1)-1) linspace2(0,1,2^(M-1)-1) 0.5});
%     
%     plotNURBS(solid,[100 0 0])
%     view(40,30)
%     axis equal
%     camlight;
%     set(gca, 'Color', 'none');
%     grid off
%     axis off
% 
%     export_fig(['../graphics/pinchedCylinderMesh' num2str(M)], '-png', '-transparent', '-r300')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scordelis-Lo Roof
% for M = 6 %[1 2 3 6]
%     R = 25; % Mid surface radius
%     L = 50;
%     phi = 40*pi/180;
%     t = 0.25;
%     
%     R = R-t/2;
%     R_o = R+t/2;
% 
%     solid = getScordeliLoRoofData(R, R_o, L/2, phi);
%     solid = insertKnotsInNURBS(solid,{linspace2(0,1,2^(M-1)-1) linspace2(0,1,2^(M-1)-1) 0.5});
%     
%     plotNURBS(solid,[100 0 0])
%     view(40,30)
%     axis equal
%     camlight;
%     set(gca, 'Color', 'none');
%     grid off
%     axis off
% %     xlabel x
% %     ylabel y
% %     zlabel z
% 
% %     export_fig(['../graphics/scordelisLoRoofMesh' num2str(M)], '-png', '-transparent', '-r300')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rectangular Prism
% for M = 5 %[1 2 3 5]
%     wx = 4;
%     wy = 2;
%     wz = 1;
%     solid = getRectangularPrismData(wx, wy, wz, [wx/2, wy/2, wz/2]);
% 
%     solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 3)...
%                                       insertUniform2(solid.knots{2}, 1) ...
%                                       insertUniform2(solid.knots{3}, 0)});
%                                   
%     solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 2^(M-1)-1)...
%                                       insertUniform2(solid.knots{2}, 2^(M-1)-1)...
%                                       insertUniform2(solid.knots{3}, 2^(M-1)-1)});
%                                   
%     plotNURBS(solid,[0 0 0])
%     axis equal
%     camlight;
%     set(gca, 'Color', 'none');
%     grid off
%     axis off
% 
%     export_fig(['../graphics/rectangularPrismMesh' num2str(M)], '-png', '-transparent', '-r300')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solid circular cylinder
% for M = 1
%     L = 5;  % Height
%     R = 1;    % Inner surface radius
%     R_o = 2;    % Outer surface radius
% 
%     solid = getSolidCylinderData(R, R_o, L);
% 
%     solid = insertKnotsInNURBS(solid,{[] [] 0.5});
%     solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 2^(M-1)-1)...
%                                       insertUniform2(solid.knots{2}, 2^(M-1)-1)...
%                                       insertUniform2(solid.knots{3}, 2^(M-1)-1)});
%     plotNURBS(solid,[100 0 0])
%     
%     axis equal
%     camlight;
%     set(gca, 'Color', 'none');
%     grid off
%     axis off
%     
%     
%     export_fig(['../graphics/solidCircularCylinderMesh' num2str(M)], '-png', '-transparent', '-r300')
% %     export_fig('../graphics/solidCircularCylinder', '-png', '-transparent', '-r300')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kneaded cylinder
% for M = 4
%     L = 5;  % Height
%     R = 1;    % Inner surface radius
%     R_o = 2;    % Outer surface radius
% 
%     solid = getSolidCylinderData(R, R_o, L);
% 
%     solid = insertKnotsInNURBS(solid,{[] [1 2 3 4]/5 []});
%     solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 2^(M-1)-1)...
%                                       insertUniform2(solid.knots{2}, 2^(M-1)-1)...
%                                       insertUniform2(solid.knots{3}, 2^(M-1)-1)});
%     plotNURBS(solid,[200 0 0])
% 
%     axis equal
%     camlight;
%     set(gca, 'Color', 'none');
%     grid off
%     axis off
% 	
%     
%     export_fig(['../graphics/kneadedCylinderMesh' num2str(M)], '-png', '-transparent', '-r300')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS circular plate
% t = 0.02;  % Thickness
% R = 2;    % Radius
% solid = getCircularPlateData(R, t);
% solid = elevateNURBSdegree(solid,[0 1 1]);
% solid = insertKnotsInNURBS(solid,{[] 0.5 []});
% plotNURBS(solid,[100 0 0])
% grid off
% axis off
% daspect([1 1 1])
% 
% lighting phong;
% camlight('right')
% set(gca, 'Color', 'none');
% view(20, 25)
% 
% export_fig('../graphics/circularPlate_2', '-png', '-transparent', '-r200')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS circular plate 2
% t = 0.02;  % Thickness
% R = 2;    % Radius
% solid = getCircularPlateData2(R, t);
% solid = elevateNURBSdegree(solid,[0 0 1]);
% solid = insertKnotsInNURBS(solid,{[0.25 0.5 0.75] [0.25 0.5 0.75] []});
% 
% plotNURBS(solid,[100 100 0])
% grid off
% axis off
% daspect([1 1 1])
% 
% lighting phong;
% camlight('right')
% set(gca, 'Color', 'none');
% view(20, 20)
% export_fig('../graphics/circularPlate2_2', '-png', '-transparent', '-r200')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS horse shoe
% R = 1; 
% L = 4; % Must be chosen such that R < L
% H = 7+sqrt(2);
% D = R;
% solid = getHorseShoeData(R, L, H, D);
% M = 4;
% %     solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 2^(M-1)-1)...
% %                                       insertUniform2(solid.knots{2}, 2^(M-1)-1)...
% %                                       insertUniform2(solid.knots{3}, 2^(M-1)-1)});
% plotNURBS(solid,[100 2 100])
% % controlPts = solid.coeffs;
% % scatter3(reshape(controlPts(1,:,:,:), 1, 50), ...
% %          reshape(controlPts(2,:,:,:), 1, 50), ...
% %          reshape(controlPts(3,:,:,:), 1, 50),'MarkerEdgeColor', 'red')
% % set(gca,'zdir','reverse')
% % grid off
% % axis off
% % camlight(50,30)
% view(135,0)
% % % export_fig('../graphics/horseShoe', '-png', '-transparent', '-r200')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spherical shell
% 
% R = 5;  % Mid surface radius
% t = 0.15;   % Thickness
% R = R-t/2; % Inner radius
% R_o = R+t/2; % Outer radius
% degenerateNURBS = 1;
% M = 3;
% degElevArr = 0;
% coreMethod = 'linear_FEM';
% for M = M
%     for degElev = degElevArr
%         p = degElev + 2;
%         close all
%         figure(1)
% %         set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
%         solid = getSphericalShellData(R, R_o);
%         fluid = extractOuterSurface(solid);
%         fluid = elevateNURBSdegree(fluid,degElev*[1 1]);
% 
%         newKnots = 2^(M-1)-1;
%         
%         fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, newKnots)...
%                                           insertUniform2(fluid.knots{2}, newKnots)});
%         fluid = repeatKnots(fluid,coreMethod);
%         fluid = degenerateIGAtoFEM(fluid,coreMethod);
%         if strcmp(coreMethod,'linear_FEM')
%             plotNURBS(fluid,[2 2])
%         else
%             plotNURBS(fluid,[200 200])
%         end
%         camlight
%         grid off
%         axis off
%         axis equal
%         drawnow
% 
% %         export_fig(['../graphics/sphericalShell/sphericalShellMesh' num2str(M) '_' num2str(p) coreMethod], '-png', '-transparent', '-r300')
% 
% %         
% %         figure(2)
% %         plotNURBS(nurbs,[100 100 0])
%         % view(50,30)
% %         camlight
% %         grid off
% %         axis off
% %         axis equal
% %         drawnow
% %         return
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Half spherical shell
% R = 5;  % Mid surface radius
% t = 0.15;   % Thickness
% R = R-t/2; % Inner radius
% R_o = R+t/2; % Outer radius
% 
% solid = getHalfSphericalShellData(R, R_o,'Xaxis');
% 
% plotNURBS(solid,[200 200 0], 1, getColor(1), 1)
% view(-40+180,30)    
% camlight
% grid off
% axis off
% axis equal
% export_fig('../graphics/sphericalShell/halfSphericalShellMesh', '-png', '-transparent', '-r300')
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS mock shell
% t = 0.0406;   % Thickness
% R_o = 4.686; % Outer radius
% R = R_o-t; % Inner radius
% L = 68;
% 
% solid = getMockShellData(R, R_o, L,'Xaxis');
% plotNURBS(solid,[100 100 0])
% view(40,30)    
% camlight
% grid off
% axis off
% axis equal
% % export_fig('../graphics/mockShell', '-png', '-transparent', '-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS model 1
% 
% t = 0.02;   % Thickness
% R_o = 3;
% R = R_o - t;
% L = 43-R_o;
% 
% % t = 0.5;   % Thickness
% % R_o = 10;
% % R = R_o - t;
% % L = 10;
% solid = getModel1Data(R, R_o, L);
% fluid = extractOuterSurface(solid);
% plotNURBS(fluid,[200 200], 1, getColor(1), 1);
% % plotNURBS(solid,[100 100 0])
% 
% view(-40+90+180, 30)
% % view(-40+90, 30)
% % view(-90, 0)
% axis equal
% grid off
% axis off
% camlight
% % n = nurbs.number(1);
% % m = nurbs.number(2);
% % l = nurbs.number(3);
% % controlPts = nurbs.coeffs;
% % scatter3(reshape(controlPts(1,:,:,1), 1, n*m), ...
% %          reshape(controlPts(2,:,:,1), 1, n*m), ...
% %          reshape(controlPts(3,:,:,1), 1, n*m),'MarkerEdgeColor', 'red')
% % xlabel x
% % ylabel y
% % zlabel z
% export_fig('../graphics/BeTSSi_M1_2', '-png', '-transparent', '-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS model 3
% setBeTSSi_M3Parameters
% if false
%     solid = getModel3Data(R_o1, R_o2, t, L);
%     fluid = extractOuterSurface(solid);
% else
%     fluid = getModel3Data2(R_o1, t, L);
% end
% plotNURBS(fluid{1},{'resolution',[20 20]});
% 
% view(-40+90, 30)
% axis equal
% grid off
% axis off
% camlight
% % 
% % export_fig('../graphics/BeTSSi_M3', '-png', '-transparent', '-r300')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS model 4
% setBeTSSi_M4Parameters
% fluid = getModel4Data(R, t);
% plotNURBS(fluid,{'resolution',[200 2]});
% 
% view(52, 32)
% axis equal
% grid off
% axis off
% camlight
% 
% export_fig('../../graphics/BeTSSi_M4_T4', '-png', '-transparent', '-r300')
% close all
% x = [0; 0; -t];
% A = [1 0 0;
%      0 1 0;
%      0 0 1];
% y = [2*R; 2*R; 0];
% y = [0; 0; 0];
% wedge1 = getWedgeData(R, t, pi/4, x+y, A);
% x = [0; -t; 0];
% A = [1 0 0;
%      0 0 1;
%      0 1 0];
% wedge2 = getWedgeData(R, t, pi/4, x+y, A);
% x = [-t; 0; 0];
% A = [0 0 1;
%      0 1 0;
%      1 0 0];
% wedge3 = getWedgeData(R, t, pi/4, x+y, A);
% A = [1 0 0;
%      0 1 0;
%      0 0 1];
% x = [(R-t)/2; -t/2; -t/2];
% prism1 = getRectangularPrismData(R+t, t, t, x+y);
% x = [-t/2; R/2; -t/2];
% prism2 = getRectangularPrismData(t, R, t, x+y);
% x = [-t/2; -t/2; R/2];
% prism3 = getRectangularPrismData(t, t, R, x+y);
% 
% hold on
% plotNURBS(wedge1,{'resolution',[200 2 2]});
% plotNURBS(wedge2,{'resolution',[200 2 2]});
% plotNURBS(wedge3,{'resolution',[200 2 2]});
% plotNURBS(prism1,{'resolution',[2 2 2]});
% plotNURBS(prism2,{'resolution',[2 2 2]});
% plotNURBS(prism3,{'resolution',[2 2 2]});
% 
% view(52, 32)
% axis equal
% grid off
% axis off
% camlight
% 
% export_fig('../../graphics/BeTSSi_M4_T3', '-png', '-transparent', '-r300')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot disk parametrizations and SINTEF/NTNU logos
% options = optimset('Display','iter','TolFun',eps,'TolX',eps);
% s_trans = fminsearchbnd(@(s)objFun(s),0.886,0,1,options);
% function f = objFun(s_trans)
% [~,Ae] = measureNURBS(getDiskData('R',1,'parm',2,'s_trans',s_trans));
% f = abs(sum(Ae(1:4))-Ae(end));
% end
% R = 1;
% for parm = 2
%     close all
%     for s_trans = 0.825354941454944 %[0.5,0.6,0.7,0.8,0.9]
%         nurbs = getDiskData('R',R,'parm',parm,'s_trans',s_trans);
%         switch parm
%             case 1
%                 colors = [31,119,180]/255;
%                 nurbs = rotateNURBS(nurbs,'theta',-pi/4,'rotAxis','Zaxis');
%                 noNewKnots = 2;
%             case 2
%                 if 0
%                     noNewKnots = 12;
%                     for i = 1:5
%                         if i == 5
%                             nurbs(i) = insertKnotsInNURBS(nurbs(i),{insertUniform2(nurbs{i}.knots{1}, noNewKnots) ...
%                                                                    insertUniform2(nurbs{i}.knots{2}, noNewKnots)}); 
%                         else
%     %                         nurbs(i) = insertKnotsInNURBS(nurbs(i),{insertUniform2(nurbs{i}.knots{1}, noNewKnots) []}); 
%                             nurbs(i) = insertKnotsInNURBS(nurbs(i),{insertUniform2(nurbs{i}.knots{1}, noNewKnots) insertUniform2(nurbs{i}.knots{2}, noNewKnots)}); 
%                         end
%                     end
%                     M = 1;
%                     plotElementEdges = true;
%                 else
% %                     noNewKnots = 1;
% %                     for i = 1:5
% %                         if i == 5
% %                             nurbs(i) = insertKnotsInNURBS(nurbs(i),{insertUniform2(nurbs{i}.knots{1}, noNewKnots) ...
% %                                                                    insertUniform2(nurbs{i}.knots{2}, noNewKnots)}); 
% %                         else
% %     %                         nurbs(i) = insertKnotsInNURBS(nurbs(i),{insertUniform2(nurbs{i}.knots{1}, noNewKnots) []}); 
% %                             nurbs(i) = insertKnotsInNURBS(nurbs(i),{[] insertUniform2(nurbs{i}.knots{2}, noNewKnots)}); 
% %                         end
% %                     end
%                     plotElementEdges = false;
% %                     plotElementEdges = true;
%                     M = 1;
%                 end
%                 nurbs = rotateNURBS(nurbs,'theta',-pi/4,'rotAxis','Zaxis');
%     %             colors = [214,39,40;
%     %                       31,119,180;
%     %                       255,127,14;
%     %                       44,160,44;
%     %                       148,103,189]/255;
%                 colors = [getColor(7); getColor(7); getColor(7); getColor(7); getColor(6)];
%                 noNewKnots = 2^(M-1)-1;
%             case 3
%                 colors = [31,119,180]/255;
%                 nurbs = rotateNURBS(nurbs,'theta',-pi/4,'rotAxis','Zaxis');
%                 noNewKnots = 4;
%         end
%         nurbs = explodedNURBS(nurbs,'L',0.02);
% %         resolution = 10*[1,1];
%         resolution = 100*[1,1];
%         % noNewKnots = 12;
%         write_g2(nurbs,'NURBSgeometries/g2files/disk')
%         nurbs = insertKnotsInNURBS(nurbs,[noNewKnots,noNewKnots,noNewKnots]);
% 
%         plotNURBS(nurbs,{'resolution',resolution,'plotControlPolygon',false,'color',colors,'setViewAngle',getView(1),'lineColor',[1,1,1],'plotElementEdges',plotElementEdges,'LineWidth',4});
%         axis equal
%         grid off
%         axis off
%         % camlight
%         figureFullScreen(gcf)
%         export_fig(['../../graphics/logo/IFEMlogo'], '-png', '-transparent', '-r600')
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS model 5
% R = 0.6/2; % Outer radius
% L = 4.8;
% l = 1;
% model = 'M5B';
% switch model
%     case 'M5A'
%         alignWithAxis = 'Xaxis';
%     case 'M5B'
%         alignWithAxis = 'Zaxis';
% end
% 
% fluid1 = getModel5Data_1(R, 0.2, 0.8, L, l, alignWithAxis);
% fluid2 = getModel5Data_2(R, 0.2, 0.8, L, l, alignWithAxis);
% plotNURBS(fluid1,[200 0], 1, getColor(1), 1);
% hold on
% plotNURBS(fluid2,[200 0], 1, getColor(1), 1);
% 
% switch model
%     case 'M5A'
%         view(60, 25)
%     case 'M5B' 
%         view(45, 32)
% end
% axis equal
% grid off
% axis off
% camlight
% 
% % export_fig(['../graphics/BeTSSi_' model], '-png', '-transparent', '-r300')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot barrel
% E = 210e9;
% rho = 7850;
% nu = 0.3;
% lossFact = 0.001;
% 
% t = 0.008;   % Thickness
% R_o = 2.5; % Outer radius
% R = R_o-t; % Outer radius
% L = 20;
% nurbs = getBarrelData(R_o,R,L,'Xaxis', [0 0 0]);
% 
% plotNURBS(nurbs,[100 100 0])
% 
% view(-40+90, 30)
% axis equal
% grid off
% axis off
% camlight
% 
% export_fig('../graphics/model3', '-png', '-transparent', '-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot NURBS spherical shell in water
% % R = 4;
% % R_o = 5;
% R = 0.8;
% R_o = 0.9;
% M = 1;
% 
% nurbs = getSphericalShellData(R, R_o,'Xaxis');
% nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, 2^(M-1)-1) ...
%                                   insertUniform2(nurbs.knots{2}, 2^(M-1)-1) ...
%                                   insertUniform2(nurbs.knots{3}, 0)});
% 
% plotNURBS(nurbs,[ceil(100/M) ceil(100/M) 0], 1, getColor(1), 1);
% camlight
% 
% hold on 
% % % 
% % R = 5;
% % R_o = 6;
% R = 0.9;
% R_o = 1;
% 
% nurbs = getSphericalShellData(R, R_o,'Xaxis');
% nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, 2^(M-1)-1) ...
%                                   insertUniform2(nurbs.knots{2}, 2^(M-1)-1) ...
%                                   insertUniform2(nurbs.knots{3}, 0)});
%                               
% h2 = plotNURBS(nurbs, [ceil(100/M) ceil(100/M) 0], 1, [187 217 238]/255, 0.4);
% 
% grid off
% axis off
% axis equal
% 
% % set(h2,'FaceAlpha',0.5)
% % ,'LineStyle','none','EdgeLighting','flat','FaceAlpha',0.5
% % testt = {'EdgeColor','none'};
% % plot([3 2],[3 54],testt)
% % % alpha(0.5)
% % view(50,30)
% % export_fig('../graphics/sphere_in_water', '-png', '-transparent', '-r400')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wine glass
% solid = getWineGlassData();
% n = solid.number(1);
% m = solid.number(2);
% l = solid.number(3);
% plotAt = [0 0;
%                  1 1;
%                  1 1];
% % plotNURBS(solid,[7 4 4], 0, [200,200,200]/255, 0.2,NaN,varCol)
% grid off
% axis off
% axis equal
% view(25,35)
% camlight
% lighting phong
% 
% drawnow
% solid.coeffs(1:3,:,:,:) = 100*solid.coeffs(1:3,:,:,:);
% nurbs = extractSurface(solid,  'zeta', 'outer');
% printNURBSToFile(nurbs{1}, ['../rhinoceros/wineGlassdata/nurbs' num2str(1) '.txt'])
% plotNURBS(nurbs{1},[10 10], 1)
% nurbs = extractSurface(solid,  'zeta', 'inner');
% printNURBSToFile(nurbs{1}, ['../rhinoceros/wineGlassdata/nurbs' num2str(2) '.txt'])
% plotNURBS(nurbs{1},[10 10], 1)
% nurbs = extractSurface(solid,  'eta', 'outer');
% printNURBSToFile(nurbs{1}, ['../rhinoceros/wineGlassdata/nurbs' num2str(3) '.txt'])
% plotNURBS(nurbs{1},[10 10], 1)
% nurbs = extractSurface(solid,  'eta', 'inner');
% printNURBSToFile(nurbs{1}, ['../rhinoceros/wineGlassdata/nurbs' num2str(4) '.txt'])
% plotNURBS(nurbs{1},[10 10], 1)
% export_fig('../graphics/wineGlass', '-png', '-transparent', '-r200')
% % 
% Xi = [0 0 0 0.5 1 1 1.1 1.9 2 2 2.2 2.5 3.2 4 4 4]/4;
% 
% 
% t = 0.002;
% 
% k = 0.001;
% 
% h1 = 0.005;
% h2 = 0.115;
% h6 = 0.125;
% h7 = 0.15; 
% h8 = 0.22;
% 
% x_1 = 1.8*t;
% y_1 = 1.4*t+h1;
% 
% x_2 = x_1;
% y_2 = h2-2*t;
% 
% 
% R_0i = 0.075/2;
% R_1i = 0.02;
% R_2i = 0.002;
% R_5i = 0.005;   % bottom of cup
% R_6i = 0.035;
% R_7i = 0.05;
% R_8i = 0.0375;
% 
% R_0o = R_0i;
% R_1o = 0.6*R_0o;
% R_4o = 0.4*x_1; %bending radius of shaft
% R_6o = R_6i+t;
% R_7o = R_7i+t;
% R_8o = R_8i+t;

% 
% p = 2;
% n = length(Xi)-p-1;
% controlPts = zeros(3,13);
% controlPts(:,1) = [ R_0i   0    1           ];
% controlPts(:,2) = [ R_1i   0   1   ];
% controlPts(:,3) = [ R_2i   h1   1   ];
% controlPts(:,4) = [ 0   h1   1           ]; %interpolation
% controlPts(:,5) = [ 0   (h1+h2)/4   1   ];
% controlPts(:,6) = [ 0   2*(h1+h2)/4   1   ];
% controlPts(:,7) = [ 0   3*(h1+h2)/4   1   ];
% controlPts(:,8) = [ 0   h2   1           ]; %interpolation
% controlPts(:,9) = [ R_5i   h2   1   ];
% controlPts(:,10) = [ R_6i   h6+t/2   1           ];
% controlPts(:,11) = [ R_7i   h7   1   ];
% controlPts(:,12) = [ R_8i   0.9*h8   1           ];
% controlPts(:,13) = [ R_8i   h8   1           ];
% 
% nurbs = createNURBSobject(controlPts,Xi);
% 
% % Plot NURBS curve
% plotNURBS(nurbs, 10000);
% hold on
% P = nurbs.coeffs(1:2,:);
% % plot(P(1,:),P(2,:), '*-','color','red')
% 
% controlPts(:,1) = [ R_0o   t    1           ];
% controlPts(:,2) = [ R_1o   t   1   ];
% controlPts(:,3) = [ x_1+k   y_1-k   1   ];
% controlPts(:,4) = [ x_1     y_1   1           ]; %interpolation
% controlPts(:,5) = [ x_1-k   y_1+k   1   ];
% controlPts(:,6) = [ R_4o   (h1+h2)/2   1           ];
% controlPts(:,7) = [ x_2-k   y_2-k   1   ];
% controlPts(:,8) = [ x_2     y_2   1           ];    %interpolation
% controlPts(:,9) = [ x_2+k   y_2+k   1   ];
% controlPts(:,10) = [ R_6o   h6   1           ];
% controlPts(:,11) = [ R_7o   h7   1           ];
% controlPts(:,12) = [ R_8o   0.9*h8   1           ];
% controlPts(:,13) = [ R_8o   h8   1           ];
% nurbs = createNURBSobject(controlPts,Xi);
% 
% 
% plotNURBS(nurbs, 10000);
% 
% % Plot control polygon
% P = nurbs.coeffs(1:2,:);
% % plot(P(1,:),P(2,:), '*-')
% set(gca,'xdir','reverse')
% 
% % xlim([-0.01 0.11])
% % ylim([0 0.12])
% axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot illustration figure 2D
% 
% Xi = [0 0 0 1 2 3 4 5 6 6 6];
% Eta = [0 0 0 1 1 1];
% 
% t = 0.5;
% R = 1;
% R_o = R+t/2;
% R = R-t/2;
% controlPts = zeros(3,8,3);
% controlPts(:,1,1) = [R 0 1];
% controlPts(:,2,1) = [R R 1];
% controlPts(:,3,1) = [0 R 1];
% controlPts(:,4,1) = [R/2 R 1];
% controlPts(:,5,1) = [-R R 1];
% controlPts(:,6,1) = [-R/2 0 1];
% controlPts(:,7,1) = [R -R 1];
% controlPts(:,8,1) = [R 0 1];
% 
% controlPts(:,1,2) = [R 0 1];
% controlPts(:,2,2) = [R R 1];
% controlPts(:,3,2) = [0 R 1];
% controlPts(:,4,2) = [R/2 R 1];
% controlPts(:,5,2) = [-R R 1];
% controlPts(:,6,2) = [-R/2 0 1];
% controlPts(:,7,2) = [R -R 1];
% controlPts(:,8,2) = [R 0 1];
% 
% controlPts(:,1,3) = [R_o 0 1];
% controlPts(:,2,3) = [R_o R_o 1];
% controlPts(:,3,3) = [0 R_o 1];
% controlPts(:,4,3) = [R_o/2 R_o 1];
% controlPts(:,5,3) = [-R_o R_o 1];
% controlPts(:,6,3) = [-R_o/2 0 1];
% controlPts(:,7,3) = [R_o -R_o 1];
% controlPts(:,8,3) = [R_o 0 1];
% 
% nurbs = createNURBSobject(controlPts,{Xi, Eta});
% 
% M = 1;
% 
% nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, 2^(M-1)-1)...
%                                   insertUniform2(nurbs.knots{2}, 2^(M-1)-1)});
% 
% 
% plotNURBS(nurbs, [100 100]);
% 
% 
% n = nurbs.number(1);
% m = nurbs.number(2);
% 
% controlPts = nurbs.coeffs;
% scatter(reshape(controlPts(1,:,:), 1, n*m), ...
%         reshape(controlPts(2,:,:), 1, n*m),'MarkerEdgeColor', 'red')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Infinite Plate With Circular Hole
% R = 1;
% L = 4;
% nurbs = getInfinitePlateWithCircularHoleData(R, L, 1);
% M = 1;
% nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, 2^(M-1)-1)...
%                                   insertUniform2(nurbs.knots{2}, 2^(M-1)-1)});
% plotNURBS(nurbs, [100 100]);
% 
% 
% n = nurbs.number(1);
% m = nurbs.number(2);
% 
% controlPts = nurbs.coeffs;
% scatter(reshape(controlPts(1,:,:), 1, n*m), ...
%         reshape(controlPts(2,:,:), 1, n*m),'MarkerEdgeColor', 'red')
% 
% 
% axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D example geometry for thesis
% nurbs = getExampleGeometry();
% M = 2;
% nurbs = insertKnotsInNURBS(nurbs,{linspace2(0,1,15)...
%                                   [0.75 0.90]});
% nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, 2^(M-1)-1)...
%                                   []});
% plotExample_InfiniteElements(nurbs, [40 40], 1);
% 
% 
% n = nurbs.number(1);
% m = nurbs.number(2);
% % 
% % controlPts = nurbs.coeffs;
% % scatter(reshape(controlPts(1,:,:), 1, n*m), ...
% %         reshape(controlPts(2,:,:), 1, n*m),'MarkerEdgeColor', 'red')
% 
% axis equal
% axis off
% grid off
% y1 = ylim;
% height = y1(2) - y1(1);
% x1 = xlim;
% width = x1(2) - x1(1);

% matlab2tikz('plotData/infiniteElements/example.tikz', 'height', '\figureheight', 'width', '\figurewidth', ...
%     'showInfo',false, 'standalone',false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of NURBS ellipse (quarter (circle)
% no2Dpoints  = 1000;
% Xi = [0 0 0 1 1 1];
% p = 2;
% n = length(Xi)-p-1;
% controlPts = zeros(3,3);
% 
% a = 1;
% b = 1;
% 
% controlPts(:,1) = [ a   0   1           ];
% controlPts(:,2) = [ a   b   0.5   ];
% controlPts(:,3) = [ 0   b   1           ];
% 
% nurbs = createNURBSobject(controlPts,Xi);
% 
% h = plotNURBS(nurbs, no2Dpoints);
% xData = get(h,'XData');
% yData = get(h,'YData');
% 
% axis equal
% 
% hold on
% theta = linspace(0,pi/2,no2Dpoints);
% 
% plot(a*cos(theta),b*sin(theta),'color','red')
% ticklabelformat(gca,'x','%1.10f')
% ticklabelformat(gca,'y','%1.10f')


% error_curves = max(abs(xData-

%% Plot NURBS model 1
% E = 210e9;
% rho = 7850;
% nu = 0.3;
% lossFact = 0.001;
% 
% t = 0.02;   % Thickness
% R_o = 3;
% R = R_o - t;
% L = 43-R_o;
% L1 = 1;
% % t = 2;   % Thickness
% % R_o = 10;
% % R = R_o - t;
% % L = 10-R_o;
% solid = getModel1Data(R, R_o, L, [0 0 L/2], 'Zaxis');
% 
% intermediateLayer = getModel1Data2(R_o+L1,L, [0 0 0],  'Xaxis', L1);
% 
% c_x = 25;
% R = 6;
% s = 100;
% % 
% % solidAndEllipsoid = embedSolidInProlateSpheroid(solid,c_x,R,s,'Zaxis');
% % 
% % solidAndEllipsoid = insertKnotsInNURBS(solidAndEllipsoid,{insertUniform2(solidAndEllipsoid.knots{1}, 1) ...
% %                                             insertUniform2(solidAndEllipsoid.knots{2}, 1) ...
% %                                             insertUniform2(solidAndEllipsoid.knots{2}, 1)});
% % solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 1) ...
% %                                   insertUniform2(solid.knots{2}, 1) ...
% %                                   insertUniform2(solid.knots{2}, 1)});
% % 
% % 
% % plotNURBS(solid,[20 20 0], 1, getColor(1), 1);
% % hold on
% % plotNURBS(solidAndEllipsoid,[20 20 20], 1, [187 217 238]/255, 0.4)
% plotNURBS(intermediateLayer,[20 20], 1, getColor(1), 1);
% 
% n = intermediateLayer.number(1);
% m = intermediateLayer.number(2);
% controlPts = zeros(n*m,3);
% count = 0;
% for j=1:m
%     controlPts(n*count+1:n*(count+1),:) = intermediateLayer.coeffs(1:3,:,j)';
%     count = count+1;
% end
% scatter3(controlPts(:,1), controlPts(:,2), controlPts(:,3),'fill', 'MarkerEdgeColor', 'red')
% camlight
% grid off
% axis off
% axis equal
% % export_fig('../graphics/model1_in_water', '-png', '-transparent', '-r400')
% 


%% Plot NURBS Model 3 in water
% 
% f_arr = 1e3; %[1 3 10 30]*1e3;
% M = 1;
% 
% eta2 = 0.8;
% eta1 = 0.265;
% 
% omega_arr = 2*pi*f_arr;
% c_f = 1500;
% t = 0.008;   % Thickness
% R_o1 = 5; % Outer radius
% R_o2 = 3; % Outer radius
% R1 = R_o1-t; % Inner radius
% R2 = R_o2-t; % Inner radius
% L = 49-R_o1-R_o2;
% 
% 
% x_0 = [-L/2-(R_o1-R_o2)/2, 0, 0]; % The origin of the model
% alignWithAxis = 'Xaxis';
% 
% s1 = 9;
% s2 = 2;
% c_z = 30; % 30
% c_xy = c_z/2.8; % 12
% L1 = c_z/19;
% nurbs = getModel3Data(R_o1, R_o2, R1, R2, L, eta1, eta2);
% nurbs2 = getModel3Data(L1/2+R_o1, L1/2+R_o2, R1, R2, L, eta1, eta2);
% nurbs3 = getModel3Data(L1+R_o1, L1+R_o2, R1, R2, L, eta1, eta2);
% nurbs4 = getModel3Data(s2*L1+R_o1, s2*L1+R_o2, R1, R2, L, eta1, eta2);
% 
% f = sqrt(c_z^2-c_xy^2);
% [R_a, ~, ~] = evaluateProlateCoords(0,c_xy,0,f);
% 
% 
% water = embedSolidInProlateSpheroid3(nurbs,nurbs3,c_z,c_xy,alignWithAxis,x_0);
% 
% %         noNewEtaKnots = 10; % 4
% %             M = 8; % 13 -> 4,918,188 elements with full refinement
% 
% x_vec1 = evaluateNURBS(water, [0, eta1, 1]);
% [~, theta_eta1, ~] = evaluateProlateCoords(x_vec1(2),x_vec1(3),x_vec1(1),f);
% 
% totArcLength1 = findArcLength(R_a,f,theta_eta1,pi);
% 
% x_vec3 = evaluateNURBS(water, [0, eta2, 1]);
% [~, theta_eta3, ~] = evaluateProlateCoords(x_vec3(2),x_vec3(3),x_vec3(1),f);
% 
% x_vec2 = evaluateNURBS(water, [0, 0.5, 1]);
% [~, theta_eta2, ~] = evaluateProlateCoords(x_vec2(2),x_vec2(3),x_vec2(1),f);
% 
% totArcLength2 = findArcLength(R_a,f,0,theta_eta3);
% 
% totArcLengthIntermediate1 = findArcLength(R_a,f,theta_eta2,theta_eta1);
% totArcLengthIntermediate2 = findArcLength(R_a,f,theta_eta3,theta_eta2);
% %         4*ceil(0.5*R_o1*pi/2*M)
% noNewXiKnots = ceil(0.55*c_xy*pi/2*M);  
% noNewZetaKnots = ceil(0.6364*(c_z-L/2-R_o2+(R_o2-R_o1)/2)*M);
% mesh1 = ceil(0.55*totArcLength1*M);
% mesh2_1 = ceil(0.55*totArcLengthIntermediate1*M);
% mesh2_2 = ceil(0.55*totArcLengthIntermediate2*M);
% mesh3 = ceil(0.55*totArcLength2*M);
% 
% findEtaKnots2
% 
% water = insertKnotsInNURBS(water,{insertUniform2(water.knots{1}, noNewXiKnots) ...
%                                   [newEtaValues1 ...
%                                   insertUniform2([eta1 0.5], mesh2_1)' insertUniform2([0.5 eta2], mesh2_2)' newEtaValues2]' ...
%                                   insertUniform2(water.knots{3}, noNewZetaKnots).^0.80});                
% plotNURBS(nurbs,[100 100 0], 1, getColor(1), 1);
% hold on
% plotNURBS(water,[20 20 20], 1, [187 217 238]/255, 0.4);
% view(30,25)
% camlight
% grid off
% axis off
% axis equal
% % export_fig(['../graphics/model3_3Dmesh_' num2str(M)], '-png', '-transparent', '-r400')
% 
% % close all
% % plotCuttingPlaneXi0_xy(water, [10 10]);
% % xlabel x
% % % zlabel z
% % axis equal
% % axis off
% % set(gca, 'Color', 'none');


%% Degenerate IGA to FEM
% 
% for M = 2 %:4
%     R = 5.075;
%     R_o = 6;
%     alignWithAxis = 'Zaxis';
% 
%     nurbs = getSphericalShellData(R, R_o, alignWithAxis);
%     elevateFEM = true;
%     if elevateFEM
%         nurbs = elevateNURBSdegree(nurbs,[1 1 2]);
%     end
%     nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, (2^(M-1)-1)) ...
%                                       insertUniform2(nurbs.knots{2}, (2^(M-1)-1)) ...
%                                       insertUniform2(nurbs.knots{3}, 0)});
% %     plotNURBS(nurbs,[ceil(2000/(2*(2^(M-1)-1))) ceil(2000/(2*(2^(M-1)-1))) 0], 1, getColor(1), 1);
% %     camlight
% %     grid off
% %     axis off
% %     axis equal
% %     export_fig(['../graphics/sphericalShellMesh' num2str(M) 'IGAequaivalent'], '-png', '-transparent', '-r300')
% %     
%     nurbs = repeatKnots(nurbs);
%     nurbs = degenerateIGAtoFEM4(nurbs);
%     close all
% 
%     plotNURBS(nurbs,[ceil(10/M) ceil(10/M) 0], 1, getColor(1), 1);
%     camlight
%     grid off
%     axis off
%     axis equal
% 
% 
%     export_fig(['../graphics/sphericalShellMesh' num2str(M) '_0_3_FEM'], '-png', '-transparent', '-r300')
% end



%% Plot ellipsoid
% 
% nurbs = getEllipsoidalData(1.2,2,4,'Xaxis');
% 
% plotNURBS(nurbs,[ceil(100) ceil(100) 0], 1, getColor(1), 1);
% camlight
% grid off
% % axis off
% axis equal
% axis on
% xlabel x
% ylabel y
% zlabel z
% % export_fig('../graphics/ellipsoid', '-png', '-transparent', '-r300')



%% Plot Torus

% nurbs = getTorusData(5,1);
% 
% plotNURBS(nurbs,[ceil(100) ceil(100) 0], 1, getColor(1), 1);
% camlight
% grid off
% axis off
% axis equal
% % export_fig('../graphics/torus', '-png', '-transparent', '-r300')



%% Plot Horseshoe

% nurbs = getHorseShoeData;
% 
% plotNURBS(nurbs,[ceil(100) 0 ceil(100)], 1, getColor(1), 1);
% camlight
% grid off
% axis off
% axis equal
% export_fig('../graphics/horseShoe', '-png', '-transparent', '-r300')

%% Plot Wineglass
% 
% for M = 6
%     close all
%     nurbs = getWineGlassData();
%     nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, 2^(M-1)-1) ...
%                                       insertUniform2(nurbs.knots{2}, 2^(M-1)-1) ...
%                                       insertUniform2(nurbs.knots{3}, 0)});
% 
%     plotNURBS(nurbs,[100 150 5], 1, getColor(1), 1);
%     camlight
%     grid off
%     axis off
%     axis equal
%     drawnow
%     export_fig(['../graphics/wineGlass_' num2str(M)], '-png', '-transparent', '-r300')
% end

%% Huang2004fde
% Xi = [0.,0.,0.,0.,0.5,0.5,1.,1.,1.,1.,];
% controlpoints = [ 260,   100, 1;
%                   100,   260, 1;
%                   260,   420, 1;
%                   420,   420, 1;
%                   580,   260, 1;
%                   420,   100, 1];
% nurbs = createNURBSobject(controlpoints.',{Xi});
% 
% plotNURBS(nurbs);
% nurbs = elevateNURBSdegree(nurbs,1);
% plotNURBS(nurbs,{'colorControlPolygon','blue'});


%% Plot Zeta knots distribution
% paramList = [0.06	0.8	0.8]; 
% PGC
% paramLists = [0.088374187687555	0.844620195968228	0.799818126847321;
%               0.068873130451583	0.869262418421568	0.809435594886491;
%               0.091622548456167	0.849129962426039	0.819526173746566;
%               0.087663923129889	0.865130371656120	0.756651847113639;
%               0.047708572324096	0.731219160949416	0.555701215384014;
%               0.030875475001067	0.766941100998776	0.625552673843656;
%               0.046840315955893	0.784405669677714	0.651131593732605;
%               0.059152175662640	0.828533898177390	0.726769976563761;
%               0.036938189461371	0.731726765315975	0.518762399991244;
%               0.022650941310735	0.702281351902585	0.431841107126532];
% paramLists = [0.0258	0.904	0.905];
% % PGU
% paramLists = [0.073772875512035	0.902047412685032	0.772462540308722;
%               0.141052138361215	0.922930845354778	0.823643667193636;
%               0.071311503068557	0.852046889428049	0.779427183767731;
%               0.035801972256923	0.911553923067910	0.857923839816364;
%               0.032729166685189	0.934741996968893	0.901149820471724;
%               0.025832537156827	0.914032663864238	0.865319269672358;
%               0.017042956112358	0.682225274920377	0.493208358502151;
%               0.013135353010521	0.643733308897193	0.407931712430845;
%               0.011955434109506	0.679689221116941	0.487539416520375;
%               0.009435630480941	0.658633686954660	0.416774132047842];
% legendArr = cell(0,1);
% col=hsv(size(paramLists,1));
% counter = 1;
% for i = 1:size(paramLists,1)
%     paramList = paramLists(i,:);
%     % paramList = [0.063501163854092	0.970932081261572	0.935889417188447]; % 200
%     % paramList = [0.022306515951600	0.781151111765398	0.691283388693613]; % 300
%     % paramList = [0.025347132391263	0.807558917125858	0.698344143696279]; % 400
%     % paramList = [0.022596067540252	0.768481490810300	0.5895451517437387]; % 500
%     % paramList = [0.021100977305892	0.764924116976478	0.560018251482435]; % 600
%     % paramList = [0.034106190380163	0.829548827790361	0.722001554899507]; % 700
%     % paramList = [0.019584555124128	0.759257551400163	0.519110346776298]; % 800
%     % paramList = [1/200, 2/7, 0.7];
%     c_1 = paramList(1);
%     c_2 = paramList(2);
%     zeta_u4 = paramList(3);
%     freq = 100*i;
% 
%     coeffs = [0, c_1, c_2, 1];
%     Zeta_u = [0 0 0 zeta_u4 1 1 1];
% 
%     zeta_vec = linspace2(0,1,1000);
%     D = ZetaKnotsDistribution(coeffs, Zeta_u, zeta_vec);
% 
%     plot(zeta_vec, D,'color',col(counter,:))
%     legendArr{counter} = ['Freq ' num2str(freq)];           
%     counter = counter + 1;
%     hold on
%     fid = fopen(['../plotData/radialKnotDistributions/f_' num2str(freq) '_PGU.txt'],'wt+','b');
%     fprintf(fid,'x\t\t\ty\n');
%     for m = 1:length(zeta_vec)
%         fprintf(fid,'%1.15f\t%1.15f\n',zeta_vec(m),D(m));
%     end
%     fclose(fid);
% 
% end
% % legend(legendArr)
% n = 4;
% controlPolygon = zeros(2,n);
% p = 2;
% for j = 1:n
%     controlPolygon(1,j) = sum(Zeta_u(j+1:j+p))/p;
% end
% controlPolygon(2,:) = coeffs;
% 
% plot(controlPolygon(1,:), controlPolygon(2,:),'*-','color','red')
% 
% freq = 300;
% % plot([0 0.5 1], [0 0.1875 1])
% % fid = fopen(['../plotData/radialKnotDistributions/f_' num2str(freq) '.txt'],'wt+','b');
% % fprintf(fid,'x\t\t\ty\n');
% % for m = 1:length(zeta_vec)
% %     fprintf(fid,'%1.15f\t%1.15f\n',zeta_vec(m),D(m));
% % end
% % fclose(fid);


%% Plot exponential integral function
% npts = 100;
% p = 3;
% k = 40;
% c_z = 20;
% c_xy = 4; % 2.5, 3.75; 
% Upsilon = sqrt(c_z^2-c_xy^2);
% chimin = c_z;
% chimax = 1.1*c_z;
% f = @(chi) radialIntegral2(4, chi, k, Upsilon, 'PGU', 1);
% 
% n_arr = ceil(2*10.^(linspace(1,3,10)));
% p_arr = 7;
% if max(p_arr) >= min(n_arr)
%     error('stop')
% end
% errArr = zeros(length(n_arr),length(p_arr));
% for j = 1:length(p_arr)
%     p = p_arr(j);
%     for i = 1:length(n_arr)
%         n = n_arr(i);
%         chiarr = linspace(chimin,chimax,n);
% 
%         splineStrct = splineInterpolation(f,chimin,chimax,n,p);
%         fapprox = @(chi)evaluateNURBS(splineStrct,(chi-chimin)/(chimax-chimin));
% 
% 
%         chi_t = (chiarr(1:end-1) + chiarr(2:end))/2;
% 
%         errArr(i,j) = sqrt(sum(abs(arrayfun(@(chi)real(f(chi)),chi_t)-arrayfun(@(chi)real(fapprox(chi)),chi_t)).^2))...
%                         /sqrt(sum(abs(arrayfun(@(chi)real(f(chi)),chi_t)).^2));
%     end
% end
% loglog(n_arr,errArr,'*-');
% 
% figure(2)
% chiarrplot = linspace(chimin,chimax,npts);
% plot(chiarrplot,arrayfun(@(chi)real(f(chi)),chiarrplot))
% hold on
% plot(chiarrplot,arrayfun(@(chi)imag(f(chi)),chiarrplot))
% plot(chiarrplot,arrayfun(@(chi)real(fapprox(chi)),chiarrplot))
% plot(chiarrplot,arrayfun(@(chi)imag(fapprox(chi)),chiarrplot))
% 
% legend({'Exact real', 'Exact imag', 'Approx real', 'Approx imag'})

%% Chebyshev approximation
% 
% f = @(x)exp(x);
% 
% B = adaptiveChebychevInterp(f,1e-14);
% 
% x = linspace(-1,1,1000);
% fappr = 0;
% for j = 1:length(B)
%     if j == 1 || j == length(B)
%         fappr =fappr + 0.5*B(j)*chebyshevT(j-1,x);
%     else
%         fappr = fappr + B(j)*chebyshevT(j-1,x);
%     end
% end
% plot(x,f(x),x,fappr)
% 

%% Test LR-Bsplines
% close all
% p_xi = 2;
% p_eta = 2;
% 
% Xi = [0, 0, 0, 1, 2, 4, 5, 6, 6, 6]/6;
% Eta = Xi;
% 
% n_xi = length(Xi)-p_xi-1;
% n_eta = length(Eta)-p_eta-1;
% controlPts = zeros(3,n_xi,n_eta);
% for j = 1:n_eta
%     for i = 1:n_xi
%         controlPts(:,i,j) = [Xi(i+p_xi); Eta(j+p_eta); 1];
%     end
% end
% 
% nurbs = createNURBSobject(controlPts,{Xi, Eta});
% 
% % plotNURBS(nurbs,[40 40]);
% % nurbs = createLRobject(nurbs);
% 
% 
% LRobj = createLRobject(nurbs);
% 
% plotLRinParameterSpace(LRobj, [100 100], 5)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot half of S1, S3, S5
% setS3Parameters
% R = R_o - t;
% % resolution = [20 20,0];
% resolution = [200 200,0];
% solidColor = [173, 178, 189]/255;
% R_a = 1.2*R_o;
% alignWithAxis = 'Yaxis';
% fluid = getHalfSphericalShellData(R_o, R_a,alignWithAxis);
% solid = getHalfSphericalShellData(R, R_o,alignWithAxis);
% fluid_i = getHalfSolidSphereData(R,alignWithAxis);
% 
% plotAt = [0, 0;
%                  0, 1;
%                  1, 1];
% plotNURBS(solid,resolution, 0, solidColor, 1, NaN, varCol)  
% view(-40+180,30)    
% camlight
% hold on
% plotAt = [0, 0;
%                  0, 1;
%                  0, 1];
% plotNURBS(fluid,resolution, 0, [33 76 161]/255, 0.8, NaN, varCol)
% plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
% plotNURBS(fluid_i,resolution, 0, [187 217 238]/255, 0.3, NaN, varCol)
% % plotNURBS(fluid_i,resolution, 0, [33 76 161]/255, 0.8, NaN, varCol)
% % view(140,30)  
% grid off
% axis off
% axis equal
% % export_fig('../graphics/sphericalShell/half_S3', '-png', '-transparent', '-r600')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot half of S13, S35, S15
% setS13Parameters
% R = R_o-t;
% % resolution = [20 20,0];
% resolution = [200 200,0];
% solidColor = [173, 178, 189]/255;
% R_a = 1.2*R_o(1);
% alignWithAxis = 'Yaxis';
% fluid = getHalfSphericalShellData(R_o(1), R_a,alignWithAxis);
% solid = getHalfSphericalShellData(R(1), R_o(1),alignWithAxis);
% 
% fluid_2 = getHalfSphericalShellData(R_o(2), R(1),alignWithAxis);
% solid_2 = getHalfSphericalShellData(R(2), R_o(2),alignWithAxis);
% fluid_i = getHalfSolidSphereData(R(2),alignWithAxis);
% 
% plotAt = [0, 0;
%                  0, 1;
%                  1, 1];
% plotNURBS(solid,resolution, 0, solidColor, 1, NaN, varCol)  
% hold on
% plotNURBS(solid_2,resolution, 0, solidColor, 1, NaN, varCol)  
% view(-40+180,30)    
% camlight
% plotAt = [0, 0;
%                  0, 1;
%                  0, 1];
% plotNURBS(fluid,resolution, 0, [33 76 161]/255, 0.8, NaN, varCol)
% plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
% % plotNURBS(fluid_2,resolution, 0,     [33 76 161]/255, 0.8, NaN, varCol)
% plotNURBS(fluid_2,resolution, 0, [187 217 238]/255, 0.3, NaN, varCol)
% 
% plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
% plotNURBS(fluid_i,resolution, 0, [187 217 238]/255, 0.3, NaN, varCol)
% % view(140,30)  
% grid off
% axis off
% axis equal
% % export_fig('../graphics/sphericalShell/half_S13', '-png', '-transparent', '-r600')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot half of S135
% 
% % setS3Parameters
% % resolution = [20 20,0];
% resolution = [200 200,0];
% % solidColor = getColor(1);
% solidColor = [173, 178, 189]/255;
% setS135Parameters
% R = R_o-t;
% R_a = 1.2*R_o(1);
% alignWithAxis = 'Yaxis';
% fluid = getHalfSphericalShellData(R_o(1), R_a,alignWithAxis);
% solid = getHalfSphericalShellData(R(1), R_o(1),alignWithAxis);
% % 
% fluid_2 = getHalfSphericalShellData(R_o(2), R(1),alignWithAxis);
% solid_2 = getHalfSphericalShellData(R(2), R_o(2),alignWithAxis);
% fluid_3 = getHalfSphericalShellData(R_o(3), R(2),alignWithAxis);
% solid_3 = getHalfSphericalShellData(R(3), R_o(3),alignWithAxis);
% fluid_i = getHalfSolidSphereData(R(3),alignWithAxis);
% 
% plotAt = [0, 0;
%                  0, 1;
%                  1, 1];
% plotNURBS(solid,resolution, 0, solidColor, 1, NaN, varCol)  
% hold on
% plotNURBS(solid_2,resolution, 0, solidColor, 1, NaN, varCol)  
% plotNURBS(solid_3,resolution, 0, solidColor, 1, NaN, varCol)  
% view(-40+180,30)    
% camlight
% plotAt = [0, 0;
%                  0, 1;
%                  0, 1];
% plotNURBS(fluid,resolution, 0, [33 76 161]/255, 0.8, NaN, varCol)
% plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
% plotNURBS(fluid_2,resolution, 0, [33 76 161]/255, 0.8, NaN, varCol)
% plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
% plotNURBS(fluid_3,resolution, 0, [187 217 238]/255, 0.3, NaN, varCol)
% 
% plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
% plotNURBS(fluid_i,resolution, 0, [187 217 238]/255, 0.3, NaN, varCol)
% % view(140,30)  
% grid off
% axis off
% axis equal
% % export_fig('../graphics/sphericalShell/half_S135', '-png', '-transparent', '-r600')
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot half of S13_ESBC and S15_ESBC

% setS15Parameters
% R = R_o-t;
% % resolution = [20 20,0];
% resolution = [200 200,0];
% solidColor = [173, 178, 189]/255;
% R_a = 1.2*R_o(1);
% alignWithAxis = 'Yaxis';
% fluid = getHalfSphericalShellData(R_o(1), R_a,alignWithAxis);
% solid = getHalfSphericalShellData(R(1), R_o(1),alignWithAxis);
% 
% fluid_2 = getHalfSphericalShellData(R_o(2), R(1),alignWithAxis);
% solid_2 = getHalfSolidSphereData(R_o(2),alignWithAxis);
% 
% plotAt = [0, 0;
%                  0, 1;
%                  1, 1];
% plotNURBS(solid,resolution, 0, solidColor, 1, NaN, varCol)  
% hold on
% plotAt = [0, 0;
%                  1, 1;
%                  1, 1];
% plotNURBS(solid_2,resolution, 0, solidColor, 1, NaN, varCol)  
% view(-40+180,30)    
% camlight
% grid off
% axis off
% axis equal
% plotAt = [0, 0;
%                  0, 1;
%                  0, 1];
% plotNURBS(fluid,resolution, 0, [33 76 161]/255, 0.8, NaN, varCol)
% plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
% plotNURBS(fluid_2,resolution, 0,     [33 76 161]/255, 0.8, NaN, varCol)
% % plotNURBS(fluid_2,resolution, 0, [187 217 238]/255, 0.3, NaN, varCol)
% % view(140,30)  
% grid off
% axis off
% axis equal
% % export_fig('../graphics/sphericalShell/S15_ESBC', '-png', '-transparent', '-r600')
% %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot half of many spherical shells
% k=NaN;
% task.R_far = NaN;
% % setS3Parameters
% % resolution = [20 20,0];
% resolution = [200 200,0];
% % solidColor = getColor(1);
% solidColor = [173, 178, 189]/255;
% % R = [5.8, 3.2, 1.9];
% % R_o = [6, 4, 2];
% setTripleShellParameters
% R = R_o-t;
% R_a = 7;
% for withSolidSphere = 1
%     close all
%     alignWithAxis = 'Yaxis';
%     fluid = getHalfSphericalShellData(R_o(1), R_a,alignWithAxis);
%     solid = getHalfSphericalShellData(R(1), R_o(1),alignWithAxis);
%     % setS2Parameters
%     fluid_2 = getHalfSphericalShellData(R_o(2), R(1),alignWithAxis);
%     solid_2 = getHalfSphericalShellData(R(2), R_o(2),alignWithAxis);
%     fluid_3 = getHalfSphericalShellData(R_o(3), R(2),alignWithAxis);
%     if withSolidSphere
%         solid_3 = getHalfSolidSphereData(R_o(3), alignWithAxis);
%     else
%         solid_3 = getHalfSphericalShellData(R(3), R_o(3), alignWithAxis);
%     end
%     fluid_i = getHalfSolidSphereData(R(3),alignWithAxis);
% 
%     plotAt = [0, 0;
%              0, 1;
%              1, 1];
%     plotNURBS(solid,{'resolution',resolution,'plotElementEdges',false, 'alphaValue',1,'color',solidColor,'plotAt',plotAt});
%     hold on
%     plotNURBS(solid_2,{'resolution',resolution,'plotElementEdges',false, 'alphaValue',1,'color',solidColor,'plotAt',plotAt});
%     plotNURBS(solid_3,{'resolution',resolution,'plotElementEdges',false, 'alphaValue',1,'color',solidColor,'plotAt',plotAt});
%     view(-40+180,30) 
%     camlight 
%     if withSolidSphere
%         view(40+180,30)    
%     end
%     plotAt = [0, 0;
%              0, 1;
%              0, 1];
%     plotNURBS(fluid,{'resolution',resolution,'plotElementEdges',false, 'alphaValue',0.8,'color',[33 76 161]/255,'plotAt',plotAt});
%     plotAt = [0, 0;
%              0, 1;
%              0, 0];
%     plotNURBS(fluid_2,{'resolution',resolution,'plotElementEdges',false, 'alphaValue',0.6,'color',[187 217 238]/255,'plotAt',plotAt});
%     plotAt = [0, 0;
%              0, 1;
%              0, 0];
%     plotNURBS(fluid_3,{'resolution',resolution,'plotElementEdges',false, 'alphaValue',0.5,'color',[33 76 161]/255,'plotAt',plotAt});
% 
%     if ~withSolidSphere
%         plotAt = [0, 0;
%                  0, 1;
%                  0, 0];
%         plotNURBS(fluid_i,{'resolution',resolution,'plotElementEdges',false, 'alphaValue',0.3,'color',[187 217 238]/255,'plotAt',plotAt});
%     end
%     % view(140,30)  
%     grid off
%     axis off
%     axis equal
%     figureFullScreen(gcf)
%     if false
%         if withSolidSphere
%             export_fig('../../graphics/sphericalShell/half_shells_2', '-png', '-transparent', '-r600')
%         else
%             export_fig('../../graphics/sphericalShell/half_shells', '-png', '-transparent', '-r600')
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Freud polynomials
% close all
% N = 27;
% useSym = true;
% r = 1;
% if r == 0
%     xmax = 3.8736*N+11;
% elseif r == 1
%     xmax = 0.335*N+6.7299;
% else
%     xmax = 3.5;
% end
% x = linspace(0,xmax,1000);
% legendArr = cell(0,1);
% if useSym
%     x = vpa(x);
%     r = vpa(r);
% end
% % ff = @(n) 0.007166666666667*n + 0.8255;
% if r < 4
%     regF = @(n) 0.9347+(-0.1461+0.008928*n)*(exp(-0.01558*(n-4)^2)*heaviside(n-4)+1-heaviside(n-4));
% else
%     regF = @(n) 0.948+(-0.1461+0.008928*n)*(exp(-0.01558*(n-4)^2)*heaviside(n-4)+1-heaviside(n-4));
% end
% if 1
%     for n = 2
%         [P, S] = FreudP(n,r,x);
%         plot(x,P)
%         hold on
%         legendArr{1} = ['n = ' num2str(n)];
%         n = n+1;
%         [P, S] = FreudP(n,r,x);
%         plot(x,P)
%         hold on
%         legendArr{2} = ['n = ' num2str(n)];
%         ylim([-4, 5])
%         xlim(double([min(x),max(x)]))
%         legend(legendArr)
%         plot(x,zeros(size(x)),'black')
%         drawnow
%         pause(5)
%         hold off
%     end
% else
%     for n = N:N
%         [P, S] = FreudP(n,r,x);
%         plot(x,P)
%         hold on
%         legendArr{end+1} = ['n = ' num2str(n)];
%     end
%     ylim([-4, 5])
%     xlim(double([min(x),max(x)]))
%     legend(legendArr)
%     plot(x,zeros(size(x)),'black')
%     hold off
% end
% % ftType = fittype('b*exp(c*x)','coeff',{'b','c'});
% % f = fit(x',Q,ftType);
% return
% figure(2)
% Q = [5,  0.835;
%      9   0.89;
%      15, 0.933;
%      20, 0.935];
% % f = fit(Q(:,1),Q(:,2),'linear');
% % 
% % plot(f,Q(:,1),Q(:,2))
% P = polyfit(Q(:,1),Q(:,2),1);
% x = linspace(Q(1,1),Q(end,1),1000);
% plot(x,P(1)*x + P(2),Q(:,1),Q(:,2),'*')
% 
% ftType = fittype('c+(a+b*x)*exp(-d*(x-4)^2)','coeff',{'a','b','c','d'});
% f = fit(Q(:,1),Q(:,2),ftType);
% % f = fit(x',Q,'exp1');
% % a = f.a;
% % b = f.b;
% plot(f,Q(:,1),Q(:,2))
% % figure(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot max Laguerre quadrature point
% N = 63;
% useSym = true;
% r = 0;
% legendArr = cell(0,1);
% N_arr = 1:N;
% maxQ = zeros(1,N);
% for n = N_arr
%     [Q, W] = gaussLaguerreQuad(n, r);
%     maxQ(n) = max(Q);
% end
% a = maxQ(end)-maxQ(end-1);
% plot(N_arr, maxQ, N_arr, a*(N_arr-N_arr(end))+maxQ(end)+11)
% figure(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot max Freud quadrature point
% N = 7;
% useSym = true;
% r = 4;
% legendArr = cell(0,1);
% N_arr = 1:N;
% maxQ = zeros(1,N);
% for n = N_arr
%     [Q, W] = gaussFreudQuad(n, r);
%     maxQ(n) = max(Q);
% end
% if true
%     a = maxQ(end)-maxQ(end-1);
%     b = -a*N_arr(end)+maxQ(end);
% else
%     a = 0.190548016121440;
%     b = 2.659598064635782;
% end
% plot(N_arr, maxQ, N_arr, a*N_arr + b)
% % figure(4)
% a
% b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Laguerre quadrature points
% close all
% N = 63;
% useSym = true;
% r = 0;
% legendArr = cell(0,1);
% N_arr = N:N;
% maxQ = zeros(1,N);
% for n = N_arr
%     [Q, W] = gaussLaguerreQuad(n, r);
%     plot(Q,'-*')
%     hold on
% end
% % a = maxQ(end)-maxQ(end-1);
% % plot(N_arr, maxQ, N_arr, a*(N_arr-N_arr(end))+maxQ(end)+11)
% figure(1)
% Q = gaussLaguerreQuad(N, r);
% x = 1:length(Q);
% n = 2;
% % p = polyfit(x,Q', n);
% % fitOptions = fitoptions('General model');
% % fitOptions.StartPoint = [1,1,1,1,1,1];
% 
% % ftType = fittype('b*exp(c*x)','coeff',{'b','c'});
% % f = fit(x',Q,ftType);
% % f = fit(x',Q,'exp1');
% % a = f.a;
% % b = f.b;
% c = 0;
% maxQ = 3.873575844419008*N+11;
% plot(x',Q,x',maxQ*exp(0.04*(x'-N)))
% % plot(f,x',Q,x')
% % 
% % P = 0;
% % for m = 0:n
% %     P = P + p(m+1)*x.^(n-m);
% % end
% % hold on
% % plot(x,P)
% %     
% % 
% % 
% % 
% % plot(Q,zeros(size(Q)),'*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot p_inc in time domain
% close all
% T = 0.03;
% % N = 200;
% c = 1500;
% npts = 100;
% P_inc = 1;
% 
% p_inc = @(z,t,N) P_inc*(1-exp(2*pi*1i*N*(t+z/c)/T))./(1-exp(2*pi*1i*(t+z/c)/T))/N;
% dt = T/npts;
% t_arr = linspace(0,T-dt,npts);
% z = linspace(0,5,10*npts);
% for i = 1:npts
%     t = t_arr(i);
%     plot(z,real(p_inc(z,t,200)),z,real(p_inc(z,t,1000)))
%     ylim([-1,1])
%     pause(0.01)
% end
% 
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of infinite element polynomials (chebychev and multipole)
% 
% E = 210e9;
% rho = 7850;
% nu = 0.3;
% lossFact = 0.001;
% 
% t = 0.008;   % Thickness
% R_o1 = 5; % Outer radius
% R_o2 = 3; % Outer radius
% R1 = R_o1-t; % Inner radius
% R2 = R_o2-t; % Inner radius
% L = 49-R_o1-R_o2;
% 
% 
% r_a = 26; % 30
% c_xy = r_a/2.9; % 2.5, 3.75
% L1 = 10*(r_a/30)^7/19;
% 
% N = 4;
% r = linspace(r_a,20*r_a,500)';
% varCol.IEbasis = 'Shifted Chebyshev';
% varCol.formulation = 'PGU';
% varCol.N = N;
% varCol = generateCoeffMatrix(varCol);
% D = varCol.D;
% 
% figure(1)
% x = linspace(0,1,500).';
% P = zeros(length(x),N);
% legendArr = cell(0,1);
% for n = 1:N
% %     fid = fopen(['../plotData/fundamentalFunctions/legendre_n' num2str(n) '.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
%     for nt = 1:N
%         P(:,n) = P(:,n) + D(n,nt)*x.^(nt-1);
%     end
%     if n > 1
%         P(:,n) = P(:,n) + 1;
%     end
%     legendArr{end+1} = sprintf('n = %d', n-1);
% end
% plot(x, P)
% legend(legendArr)
% hold on 
% plot(xlim,[0,0],'color','black')

% 
% figure(2)
% x = r_a./r;
% P = zeros(length(r),N);
% legendArr = cell(0,1);
% for n = 1:N
% %     fid = fopen(['../plotData/fundamentalFunctions/legendre_n' num2str(n) '.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
%     for nt = 1:N
%         P(:,n) = P(:,n) + D(n,nt)*x.^nt;
%     end
%     legendArr{end+1} = sprintf('n = %d', n);
% end
% plot(r, P)
% legend(legendArr)
% hold on 
% plot(xlim,[0,0],'color','black')
% 
% 
% 
% figure(3)
% varCol.IEbasis = 'Multipole';
% varCol = generateCoeffMatrix(varCol);
% D = varCol.D;
% x = r_a./r;
% P = zeros(length(r),N);
% legendArr = cell(0,1);
% for n = 1:N
% %     fid = fopen(['../plotData/fundamentalFunctions/legendre_n' num2str(n) '.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
%     for nt = 1:N
%         P(:,n) = P(:,n) + D(n,nt)*x.^nt;
%     end
%     legendArr{end+1} = sprintf('n = %d', n);
% end
% plot(r, P)
% legend(legendArr)
% hold on 
% plot(xlim,[0,0],'color','black')
% 
% 
% figure(4)
% varCol.IEbasis = 'Bernstein';
% varCol = generateCoeffMatrix(varCol);
% D = varCol.D;
% x = r_a./r;
% P = zeros(length(r),N);
% legendArr = cell(0,1);
% for n = 1:N
% %     fid = fopen(['../plotData/fundamentalFunctions/legendre_n' num2str(n) '.txt'],'wt+','b');
% %     fprintf(fid,'x\t\t\ty\n');
%     for nt = 1:N
%         P(:,n) = P(:,n) + D(n,nt)*x.^nt;
%     end
%     legendArr{end+1} = sprintf('n = %d', n);
% end
% plot(r, P)
% legend(legendArr)
% hold on 
% plot(xlim,[0,0],'color','black')





%% Plot infinite element figure from prolate spheroid (BeTSSi submarine)
% setBCParameters
% L_gamma = L+a+g2+g3;
% 
% 
% c_z = 35; % 26
% c_xy = c_z/2.9*5/6; % c_z/2.9
% 
% eta2 = 0.8*(c_z/30)^(1-1.0);
% eta1 = 0.265*(c_z/30)^(1-1.6);
% 
% x_0 = [-L_gamma/2+a; 0; 0]; % The origin of the model
% M = 4;
% open('../results/geometries/BeTSSi.fig')
% hold on
% 
% solid = getEllipsoidalData(c_z,c_xy,c_xy,'Xaxis',x_0);
% % solid = getEllipsoidalData(c_xy,c_xy,c_z,'Zaxis',x_0);
% 
% 
% solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 2^(M-1)-1) ...
%                                   insertUniform2(solid.knots{2}, 2*2^(M-1)-1)});
% 
% % solid = repeatKnots(solid);
% 
% % h1 = plotNURBSinfElement(solid,[10 10], 1, [0,80,158]/256, 0.8, c_xy,c_z);
% h1 = plotNURBSinfElement(solid,[20 20], 1, [187 217 238]/255, 0.5, c_xy,c_z,x_0);
% hold off
% 
% view(18,10)
% camlight
% grid off
% axis off
% axis equal
% % export_fig('../../graphics/BeTSSi_in_waterInf', '-png', '-transparent', '-r600')
% % axis on
% 
% %% Plot error on sphere discretization
% % R = 1;
% % alignWithAxis = 'Xaxis';
% % fluid = getEllipsoidalData(R,R,R,alignWithAxis);
% % 
% % plotNURBS(fluid,[100 100], 1, getColor(1), 1, NaN);
% % 
% % % camlight
% % grid off
% % axis off
% % axis equal
% % colorbar
% % 
% 


%% Plot infinite element figure from prolate spheroid (BeTSSi model 3)
% setBeTSSi_M3Parameters
% % L_gamma = L+a+g2+g3;
% 
% 
% R_max = max(R_o1,R_o2);
% s = 0.25 + 0.05*(L-R_max)/R_max;
% c_z = (L+R_o1+R_o2)/2 + s*R_max;
% c_xy = R_max + s*R_max;
% % c_z = 35; % 26
% % c_xy = c_z/2.9*5/6; % c_z/2.9
% 
% eta2 = 0.8*(c_z/30)^(1-1.0);
% eta1 = 0.265*(c_z/30)^(1-1.6);
% 
% x_0 = [-L/2; 0; 0]; % The origin of the model
% nurbs = getModel3Data(R_o1, R_o2, t, L);
% fluid = extractOuterSurface(nurbs);
% plotNURBS(fluid{1},{'resolution',[10 10], 'elementBasedSamples',true,'samplingDistance',0.1});
% hold on
% 
% solid = getEllipsoidalData(c_z,c_xy,c_xy,'Xaxis',x_0);
% % solid = getEllipsoidalData(c_xy,c_xy,c_z,'Zaxis',x_0);
% 
% 
% M = 4;
% solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 2^(M-1)-1) ...
%                                   insertUniform2(solid.knots{2}, 2*2^(M-1)-1)});
% 
% % solid = repeatKnots(solid);
% 
% % h1 = plotNURBSinfElement(solid,[10 10], 1, [0,80,158]/256, 0.8, c_xy,c_z);
% h1 = plotNURBSinfElement(solid,[20 20], 1, [187 217 238]/255, 0.5, c_xy,c_z,x_0);
% hold off
% 
% view(18,10)
% camlight
% grid off
% axis off
% axis equal
% figureFullScreen(gcf)
% export_fig('../../graphics/model3_in_waterInf', '-png', '-transparent', '-r300')
% % axis on
% 

%% Plot infinite element figure from prolate spheroid (Non Sparable Geometries)

% setBeTSSi_M3Parameters
% 
% c_z = (L+R_o1+R_o2)/2;
% %     c_xy = (R_o1+R_o2)/2; % 2.5, 3.75; 
% c_xy = 3*sqrt(665)*(1/20); % 2.5, 3.75; 
% % c_z = 24.5;
% % c_xy = 4; % 2.5, 3.75
% 
% totLength = R_o1*pi/2 + sqrt(L^2+(R_o1-R_o2)^2)+ R_o2*pi/2;
% eta1 = R_o1*pi/2/totLength;
% eta2 = (R_o1*pi/2 + sqrt(L^2+(R_o1-R_o2)^2))/totLength;
% 
% x_0 = -[L/2+(R_o1-R_o2)/2; 0; 0]; % The origin of the model
% M = 3;
% nn = 2^(M-1)-1;
% solid = getModel3Data(R_o1,R_o2,t,L);
% % solid = translateNURBS(solid,x_0);
% solid = extractOuterSurface(solid);
% solid = insertKnotsInNURBS(solid{1},{[] linspace2(eta1, eta2, 5) []});
% solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
%                                   insertUniform2(solid.knots{2}, nn) []});
% 
% plotNURBS(solid,{'resolution',[100 100], 'elementBasedSamples',true,'samplingDistance',0.1,'alphaValue',1});
% hold on
% % ellipsoid = getEllipsoidalData(c_z,c_xy,c_xy,'Xaxis', x_0);
% % plotNURBS(ellipsoid,{'resolution',[20 40],'alphaValue',0.6,'color','blue'});
% 
% hold on
% 
% solid = repeatKnots(solid,'C0_IGA');
% endpoint = evaluateNURBS(solid{1},[0,1]);
% 
% h1 = plotNURBSinfElement2(solid{1},[40 40], 1, getColor(1), 1, c_xy,c_z,x_0);
% hold off
% 
% view(18,10)
% camlight
% grid off
% axis off
% axis equal
% figureFullScreen(gcf)
% export_fig('../../graphics/model3_in_waterInf2', '-png', '-transparent', '-r300')
% 
%% Plot error on sphere/ellipsoid discretization
% R = [1,2,3];
% % R = ones(1,3);
% x_0 = rand(1,3);
% % x_0 = zeros(1,3);
% alpha = 0.432;
% % alpha = 0;
% fluid = getEllipsoidalData({'R',R, ...
%                              'alignWithAxis', 'Yaxis', ...
%                              'x_0',x_0,...
%                              'alpha', alpha,...
%                              'parm', 2, ...
%                              't', 0.1, ...
%                              'iXi', [1,1,2,2]/3,...
%                              'iEta', [1,1]/2});
% for i = 1:numel(fluid)
%     plotNURBS(fluid{i},{'resolution',[20 20 20], 'plotControlPolygon', false, 'colorFun', @(x) log10(abs(norm((x-x_0)./R)-1))});
% end
% 
% % camlight
% grid off
% axis off
% axis equal
% view(18,30)
% colorbar
% % addpath(genpath('../export_fig/matlab2tikz'))
% % addpath(genpath('../export_fig/matlab2tikz/src/private'))
% % 
% % extraAxisOptions = {...
% %     'axis equal image', ...
% %     'xtick={-3.141592654,-1.570796327,0,1.570796327,3.141592654}', ...
% %     'xticklabels={$-\PI$,$-\PI/2$,0,$\PI/2$,$\PI$}', ...
% %     'ytick={-1.570796326794897,0,1.570796326794897}', ...
% %     'yticklabels={$-\PI/2$,0,$\PI/2$}', ...
% %     'xlabel=$\xi$', ...
% %     'ylabel=$\eta$'};
% % matlab2tikz('test.tex', 'height', '4in', ...
% %         'extraAxisOptions', extraAxisOptions, 'relativeDataPath', 'contents/kirchhoff/') % , 'imagesAsPng', false
% % matlab2tikz('../test.tex') % , 'imagesAsPng', false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getSphericalShellDataPatchedQuarter Errors
% R_o = sym('1');
% t = sym('0.01');
% R_o = 1;
% t = 0.01;
% 
% nurbs = getSphericalShellDataPatchedQuarter(R_o);
% nurbs = elevateDegreeInPatches(nurbs,[1 2]);
% npts = 29;
% xi = linspace(0,R_o,npts);
% eta = linspace(0,R_o,npts);
% [XI,ETA] = meshgrid(xi,eta);
% XI = reshape(XI,npts^2,1);
% ETA = reshape(ETA,npts^2,1);
% X = evaluateNURBS_2ndDeriv(nurbs{1}, [XI,ETA]);
% C = abs(norm2(X)-1);
% % return
% X = double(X);
% C = double(C);
% surf(reshape(X(:,1),npts,npts),reshape(X(:,2),npts,npts),reshape(X(:,3),npts,npts),reshape(C,npts,npts))
% colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solid sphere new parametrization
% R = 1;
% npts = 200;
% colorFun = @(x) log10(abs(norm(x)-R));
% nurbs = getSphericalShellDataPatched(R,0.1);
% nurbs = insertKnotsInPatches(nurbs,1,2,3);
% nurbs = elevateDegreeInPatches(nurbs,[3 1 2]);
% fluid = extractOuterSurface(nurbs);
% for i = 1:6
%     plotNURBS(fluid{i},{'resolution', [npts npts],'colorFun',colorFun});
% end
% % plotNURBS(fluid{1},{'resolution', [npts npts]});
% plotControlPts(fluid)
% view(100,20)
% camlight
% grid off
% axis off
% axis equal
% hold on
% % export_fig('../graphics/MA8502/Sphere2controlPolygon', '-png', '-transparent', '-r300')
% return
% solid = getSolidSphereData2(R, 'Zaxis');
% plotNURBS(solid,{'resolution',[npts npts npts], 'colorFun', @(x) log10(abs(norm(x)-1))});
% view(100,20)
% camlight
% grid off
% axis off
% axis equal
% colorbar
% % export_fig('../graphics/sphericalShell/Sphere2error', '-png', '-transparent', '-r600')
% % 
% if false
%     figure(2)
%     plotNURBS(solid,{'resolution', [npts npts npts]});
%     view(100,20)
%     camlight
%     grid off
%     axis off
%     axis equal
%     % export_fig('../graphics/sphericalShell/Sphere2', '-png', '-transparent', '-r600')
% 
%     figure(3)
%     plotNURBS(solid,{'resolution', [npts npts npts]});
%     view(100,20)
%     camlight
%     grid off
%     axis off
%     axis equal
%     hold on
%     plotControlPts(nurbs)
%     % export_fig('../graphics/sphericalShell/Sphere2controlPolygon', '-png', '-transparent', '-r600')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check ellipsoidal parametrization
% R = 1;
% npts = 100;
% x_0 = [1,1,exp(1)]*pi;
% a = 5.2;
% b = 1.6;
% c = 1.2;
% 
% varCol.colorFun = @(v) abs((v(1)-x_0(1))^2/a^2 + (v(2)-x_0(2))^2/b^2 + (v(3)-x_0(3))^2/c^2 - 1);
% 
% alignWithAxis = 'Xaxis';
% % nurbs = getEllipsoidalData3(a,c,alignWithAxis,x_0);
% % nurbs = getEllipsoidalData(a,b,c,alignWithAxis,x_0, 0,[1,1,2,2,3,3]/4, [1,1]/2);
% % nurbs = getEllipsoidalData(a,b,c,alignWithAxis,x_0, 0,[1,1,2,2,2.5,3,3]/4, [0.3,0.3,0.5,0.7,0.7]);
% nurbs = getEllipsoidalData(a,b,c,alignWithAxis,x_0, 0, [1,1,2,2]/3, [1,1]/2);
% % nurbs = insertKnotsInNURBS(nurbs,{[] linspace2(0, 0.5, 3)});
% plotNURBS(nurbs,[npts npts], 1, getColor(1), 1, NaN, varCol);
% view(100,20)
% camlight
% grid off
% % axis off
% xlabel x
% ylabel y
% zlabel z
% axis equal
% colorbar
% export_fig('../graphics/sphericalShell/Ellipsoiderror', '-png', '-transparent', '-r600')

% figure(2)
% plotNURBS(nurbs,[npts npts npts]);
% view(100,20)
% camlight
% grid off
% axis off
% axis equal
% % export_fig('../graphics/sphericalShell/Ellipsoid', '-png', '-transparent', '-r600')
% 
% figure(3)
% plotNURBS(nurbs,[npts npts npts]);
% view(100,20)
% camlight
% grid off
% axis off
% axis equal
% hold on
% plotControlPts(nurbs)
% % export_fig('../graphics/sphericalShell/Ellipsoid', '-png', '-transparent', '-r600')
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot BeTSSi from igs file
% addpath '../../../Documents/MATLAB/Add-Ons/Toolboxes/IGES Toolbox/code/'
% close all
% [ParameterData,EntityType,numEntityType,unknownEntityType,numunknownEntityType]=iges2matlab('../rhinoceros/BeTSSi.igs');
% fignr = 1;
% srf = true;
% subd = 100;
% holdoff_flag = false;
% fine_flag = 1; % 3 being the finest (0 is default)
% plotCrvPnts = false;
% srfClr = getColor(1);
% % plotIGES(ParameterData,srf,fignr,subd,holdoff_flag,fine_flag,plotCrvPnts,srfClr);
% 
% ParameterData2 = {};
% counter = 0;
% for i = 1:length(ParameterData)
% %     ParameterData{i}.name
%     if strcmp(ParameterData{i}.name, 'B-NURBS SRF')
%         ParameterData2{counter+1} = ParameterData{i};
%         ParameterData2{counter+1}.nurbs.degree = ParameterData{i}.nurbs.order-1;
%         coeffs = ParameterData{i}.nurbs.coefs;
%         for ii = 1:size(coeffs,2)
%             for jj = 1:size(coeffs,3)
%                 coeffs(1:3,ii,jj) = coeffs(1:3,ii,jj)/coeffs(4,ii,jj);
%             end
%         end
%         ParameterData2{counter+1}.nurbs.coeffs = coeffs;
%         ParameterData2{counter+1}.nurbs.type = '3Dsurface';
%         counter = counter + 1;
%     end
% end
% axis off
% axis equal
% view(63,16)
% camlight
% 
% for i = 1:length(ParameterData)
%     plotNURBS(ParameterData2{i}.nurbs,[100 100], 1, getColor(1), 1);
%     hold on
%     drawnow
% end
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Triangularization of sphere
% addpath otherFunctions/spheretri
% nPoints = 1000; 
% tic;
% [P,tri] = spheretri(nPoints);
% noElems = size(tri,1)
% % return
% toc; 
% trisurf(tri,P(:,1),P(:,2),P(:,3), 'FaceColor', getColor(1))
% axis off
% axis equal
% camlight
% 
% trir = reshape(tri,noElems*3,1);
% indices = and(any(reshape(P(trir,1),noElems,3) > 0,2), any(reshape(P(trir,1),noElems,3) < 0,2));
% indicesP = unique(tri(indices,:));
% indices2 = indicesP(P(indicesP,1) < 0);
% t = atan2(P(indices2,3),P(indices2,2));
% hold on
% % trisurf(tri(indices,:),P(:,1),P(:,2),P(:,3), 'FaceColor', 'red')
% % plot3(P(indices2,1),P(indices2,2),P(indices2,3),'*')
% P2 = P;
% % P2(indices2,:) = [zeros(numel(indices2),1), cos(t), sin(t)];
% figure(42)
% 
% trisurf(tri,P2(:,1),P2(:,2),P2(:,3), 'FaceColor', getColor(1))
% axis off
% axis equal
% camlight
% ax = gca;               % get the current axis
% ax.Clipping = 'off';    % turn clipping off
% figureFullScreen(gcf)
% % export_fig('../../graphics/sphericalShell/triangles', '-png', '-transparent', '-r300')
% return
% % P_inc = 1;
% % varCol.P_inc = P_inc;
% % k = 10.^linspace(-1,4,1000); % 10.^linspace(-1,2,1000)
% % 
% % varCol.k = k;
% % X = [1,0,0];
% % p = kirchApprTri(tri,P,X,varCol);
% % p2 = kirchApprTri(tri,P2,X,varCol);
% % TS = 20*log10(abs(p/P_inc));
% % TS2 = 20*log10(abs(p2/P_inc));
% % figure(2)
% % semilogx(k,TS,k,TS2)
% % p_ref = exactKDT(varCol.k,1,1);
% % figure(3)
% % p3 = readLaTeXFormat('results/_studies/Fillinger/S1_KDT_M6_linear_FEM_error_pVSk.txt');
% % loglog(k,100*abs(p_ref-p)./abs(p_ref),k,100*abs(p_ref-p2)./abs(p_ref),p3(:,1),p3(:,2))
% % legend('unmatched','matched')
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rigid sphere scattering formula limit k -> inf

% k = linspace(0.1,10,1000);
% k = 100000;
% p_0 = zeros(size(k))+1e-300;
% P_0 = 1;
% R_o = 1;
% n = 0;
% temp = inf;
% while max(abs(temp./p_0)) > 1e-12
%     temp = (2*n+1)*(-1)^n./(1+1i*dbessel_s(n,k*R_o,2,NaN)./dbessel_s(n,k*R_o,1,NaN));
%     p_0 = p_0 + temp;
%     n = n + 1;
% end
% p_0 = -P_0*1i^(-1)./k.*p_0;
% abs(p_0)
% return
% plot(k,abs(p_0),'DisplayName','new')
% legend('show')
% 
% 
% P_inc = 1; % Amplitude of incident wave
% rho_f = 1000;
% c_f = 1500;
% 
% omega = c_f*k; % Angular frequency
% 
% options = struct('d_vec', [0,0,1]',... 
%                  'omega', omega, ...
%                  'P_inc', P_inc, ...
%                  'rho_f', rho_f, ...
%                  'calc_farField', 1, ...
%                  'c_f', c_f);
%              
%              
% v = R_o(1)*[0,0,1];
% data = e3Dss(v, options);
% hold on
% plot(k,abs(data.p))
% figure(2)
% semilogy(k,abs(data.p-p_0)./abs(p_0))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rigid sphere scattering Kirchhoff approximation

% R_o = 1;
% P_inc = 1; % Amplitude of incident wave
% rho_f = 1000;
% c_f = 1500;
% if true
%     k = 10.^linspace(-1,3,3000);
% else
%     % compensate for the half size of the sphere)
%     eigenValues = [3.141592653589794 % pi
%                        6.283185307179587 % 2*pi
%                        9.424777960769379 % 3*pi
%                        4.493409457909065
%                        7.725251836937707
%                        5.763459196894550
%                        9.095011330476353
%                        6.987932000500519
%                        8.182561452571243
%                        9.355812111042747]';
% 
%     % k = 10.^linspace(-1,4,1000); % 10.^linspace(-1,2,1000)
%     k = sort([eigenValues, k]);
% end
% omega = c_f*k; % Angular frequency
% options = struct('d_vec', [0,0,-1]',... 
%                  'omega', omega, ...
%                  'P_inc', P_inc, ...
%                  'rho_f', rho_f, ...
%                  'calc_farField', 1, ...
%                  'c_f', c_f);
%              
%              
% v = R_o(1)*[0,0,1];
% data = e3Dss(v, options);
% % % hold on
% % % plot(k,abs(data.p))
% % figure(3)
% % p_0 = exactKDT(k,P_inc,R_o);
% % E = 100*abs(data.p-p_0)./abs(data.p);
% % loglog(k,E)
% % hold on
% % % printResultsToFile2('results/_studies/Fillinger/Kirchhoff', k.', E.', 'k', 'error_p');
% p_0 = P_inc*R_o/2*exp(-2*1i*k*R_o);
% E = 100*abs(data.p-p_0)./abs(data.p);
% loglog(k,E)
% k = k';
% % % printResultsToFile2('results/_studies/Fillinger/asymptotic', k.', E.', 'k', 'error_p');
% % E = abs(data.p-p_0)./abs(p_0);
% % loglog(k,E)
% % printResultsToFile2('results/_studies/Fillinger/asymptotic2', k.', E.', 'k', 'error_p');
% % return
% 
% N_max = 1219;
% tiny = 1e-200;
% Eps = eps;
% nFreqs = numel(k);
% nExtraTerms = 2;
% hasCnvrgd = zeros(nFreqs,nExtraTerms); % matrix of element that "has converged"
% flag = zeros(size(hasCnvrgd,1),1); % Program terminated successfully unless error occurs (for each frequency)
% E = zeros(size(k));
% p_02 = zeros(size(k));
% p2 = zeros(size(k));
% n = 0;
% indices = (1:nFreqs)';
% while n >= 0 && n <= N_max
%             hasCnvrgdTmp = zeros(length(indices),1); % temporary hasCnvrgd matrix
%     k_temp = k(indices);
%     p_0temp = (2*n+1)*(-1)^n*bessel_s(n,2*k_temp,1)*1i^n;
%     ptemp = -(2*n+1)*(-1)^n*2/1i./k_temp./(1+1i*dbessel_s(n,k_temp,2,NaN)./dbessel_s(n,k_temp,1,NaN));
%     E_n = p_0temp - ptemp;
%     p_02(indices) = p_02(indices) + p_0temp;
%     p2(indices) = p2(indices) + ptemp;
%     E(indices) = E(indices) + E_n;
%             hasCnvrgdTmp2 = abs(E_n)./(abs(E(indices))+tiny) < Eps;
%     
%             hasCnvrgdTmp(logical(hasCnvrgdTmp2)) = 1;
%             
%             hasCnvrgd(indices,:) = [hasCnvrgd(indices,2:end), hasCnvrgdTmp];
%             indicesPrev = indices;
%             indices = find(~prod(hasCnvrgd,2));
%             if isempty(indices) % every element has converged
%                 break;
%             end
%             n = n + 1;
% end
% hold on
% loglog(k,abs(E))
%     





