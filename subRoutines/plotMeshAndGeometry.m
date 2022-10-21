function task = plotMeshAndGeometry(task)
prePlot = task.prePlot;
noDomains = numel(task.varCol);
if prePlot.plot3Dgeometry
    task.prePlot.fig3Dplot = figure('Color','white','name',['3D plot of geometry with mesh ' num2str(task.msh.M)]);
    if strcmp(task.misc.method, 'IENSG') && prePlot.plotArtificialBndry
        c_z = task.varCol{1}.c_z;
        c_xy = task.varCol{1}.c_xy;
        alignWithAxis = task.varCol{1}.alignWithAxis;
        x_0 = task.iem.x_0;
        ellipsoid = getEllipsoidData('C',[c_xy,c_xy,c_z]*task.iem.A_2,'alignWithAxis', alignWithAxis, 'x_0', x_0,'Xi',task.msh.Xi);
        if prePlot.alphaValue == 1
            prePlot.alphaValue = 0.8;
        end
        if prePlot.plotFullDomain
            plotNURBSvec(ellipsoid,'resolution',[100 100],'alphaValue',0.6,'color','blue');
        end
        noTopsets = numel(prePlot.plotSubsets);
        for i = 1:noTopsets
            if strcmp(prePlot.plotSubsets{i},'xy')
                theta = linspace(0,2*pi,10000);
                x = c_z*cos(theta);
                y = c_xy*sin(theta);
                plot(x,y,'black','DisplayName','Ellipse')
            end
        end
    end
    if iscell(prePlot.color)
        colorsCell = prePlot.color;
    else
        colorsCell = cell(noDomains,1);
        if isempty(prePlot.color)
            for j = 1:noDomains
                switch task.varCol{j}.media
                    case 'fluid'
                        if task.varCol{j}.boundaryMethod
                            colorsCell{j} = getColor(1);
                        else
                            colorsCell{j} = getColor(10);
                        end
                    case 'solid'
                        colorsCell{j} = getColor(1);
                end
            end
        else
            [colorsCell{:}] = deal(prePlot.color);
        end
    end
            
    for j = 1:noDomains
        nurbs = task.varCol{j}.nurbs;
        prePlot.displayName = ['Domain ' num2str(j)];
        prePlot.color = colorsCell{j};
        if prePlot.plotFullDomain
            plotNURBSvec(nurbs, prePlot);
        end
        if isfield(task.varCol{j},'geometry') && ~isempty(prePlot.plotSubsets)
            topset = task.varCol{j}.geometry.topologysets.set;
            noTopsets = numel(prePlot.plotSubsets);
            colors = jet(noTopsets);
            colors(1:size(prePlot.color,1),:) = prePlot.color;
            for i = 1:noTopsets
                idx = findSet(topset,prePlot.plotSubsets{i});
                if isnan(idx)
                    warning(['Subset ' prePlot.plotSubsets{i} ' does not exist for domain ' num2str(j)])
                    continue
                end
                noPatches = numel(topset{idx}.item);
                nurbs = cell(1,noPatches);
                for ii = 1:numel(topset{idx}.item)
                    at = zeros(2,3,'logical');
                    patch = topset{idx}.item{ii}.Attributes.patch;
                    midx = topset{idx}.item{ii}.Text;
                    at(midx) = true;
                    outwardPointingNormals = strcmp(topset{idx}.Attributes.normal, 'outward');
                    inwardPointingNormals = strcmp(topset{idx}.Attributes.normal, 'inward');
    
                    nurbs(ii) = subNURBS(task.varCol{j}.nurbs(patch),'at',at.','outwardPointingNormals',outwardPointingNormals,'inwardPointingNormals',inwardPointingNormals);
                end
                prePlot.displayName = ['Domain ' num2str(j) ', ' topset{idx}.Attributes.name];
                if strcmp(topset{idx}.Attributes.name,'Gamma')
                    prePlot.color(j,:) = getColor(1);
                else
                    prePlot.color = colors(i,:);
                end
                plotNURBSvec(nurbs, prePlot);
            end
        end
    end
    axis equal
    axis(prePlot.axis)
    if ~strcmp(prePlot.axis,'off')
        xlabel(prePlot.xlabel)
        ylabel(prePlot.ylabel)
        zlabel(prePlot.zlabel)
    end
    set(gca, 'Color', 'none');
    if ~isempty(prePlot.title)
        title(prePlot.title)
    end
    drawnow
    camproj(prePlot.camproj)
    
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    if prePlot.useCamlight
        camlight
        material dull
    end
    figName = prePlot.export_fig_name3D;
    if isempty(figName)
        figName = task.saveName;
    end
    figName = [task.resultsFolder '/' figName];
    if exist('../export_fig', 'dir')
%         ax.SortMethod='ChildOrder';
        export_fig(figName, prePlot.format, '-transparent', prePlot.pngResolution)
    end
    savefig([figName, '.fig'])
end   
if isa(prePlot.addCommands,"function_handle")
    prePlot.addCommands()
end


% plotControlPts(fluid);
% dofsToRemove = task.varCol{1}.dofsToRemove;
% controlPts = task.varCol{1}.controlPts;
% plot3(controlPts(dofsToRemove,1),controlPts(dofsToRemove,2),controlPts(dofsToRemove,3),'o','color','blue','MarkerFaceColor','blue')

