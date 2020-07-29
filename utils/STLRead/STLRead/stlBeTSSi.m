% Copyright 2011 The MathWorks, Inc.
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To extract BeTSSi data run "COMSOL multiphysics with MATLAB" and run the following commands:
% for res = [1,2,4,8,16,32] % [1,2,4,8,16,32,64]
%     filename = ['../rhinoceros/BeTSSi/BeTSSi_mod2_BEM_res' num2str(res)];
%     model = mphload([filename '.mph']);
%     [meshstats,meshdata] = mphmeshstats(model);
%     save([filename '_org'],'meshdata')
% end
%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
% filename = '../../../../../rhinoceros/BeTSSi/BeTSSi_mod_gmsh_res4.stl';
% filename = '../../../../../rhinoceros/BeTSSi/BeTSSi_mod_COMSOL_res1.stl';
% filename = '../../../../../rhinoceros/BeTSSi/BeTSSi_mod_COMSOL_res100.stl';
% filename = '../../../../../rhinoceros/BeTSSi/BeTSSi_mod_10cm.stl';
% filename = '../../../../../rhinoceros/test.stl';
% filename = '../../../../../../FFI/BeTSSiIIb/Model_Files/outer_hull/hull_20cmMesh.stl';
% ../../../rhinoceros/BeTSSi/BeTSSi_mod_COMSOL_res1.stl
% fv = stlread(filename);
alpha_threshold = 125;
plotBetSSi = 0;
resArr = [1,2,4,8,16,32];
% return
for res = resArr
    filename = ['C:/Users/Zetison/Dropbox/work/rhinoceros/BeTSSi/BeTSSi_mod2_BEM_res' num2str(res) '_org'];
    load(filename)
    meshdataOrg = meshdata;
    clear meshdata
    for i = 1:numel(meshdataOrg.types)
        if strcmp(meshdataOrg.types{i},'tri')
            idx = i;
            break
        end
    end
    fv.vertices = meshdataOrg.vertex.';
    fv.faces = uint32(meshdataOrg.elem{idx}.'+1);
%     [noElems, dofs, h_max, alpha_max, alpha_min, aspect, waterTight, area_min, area_max, i, Pi] = computeMeshData(fv);
%     if plotBetSSi
%         plotTriangulartion(fv,alpha_threshold)
%     end
    fv = fixBeTSSiTriangularization(fv,alpha_threshold);
    [noElems, dofs, h_max1, h_max2, alpha_max, alpha_min, aspect, waterTight, area_min, area_max, skewness, i, Pi] = computeMeshData(fv);
    %% Render
    % The model is rendered with a PATCH graphics object. We also add some dynamic
    % lighting, and adjust the material properties to change the specular
    % highlighting.

    counter = 1;
    while alpha_threshold < alpha_max
%         if plotBetSSi
%             plotTriangulartion(fv,alpha_threshold)
%         end
        fv = insertNode(fv,alpha_threshold);
        [noElems, dofs, h_max1, h_max2, alpha_max, alpha_min, aspect, waterTight, area_min, area_max, skewness, i, Pi] = computeMeshData(fv);
        counter = counter + 1;
    end
    if plotBetSSi
        plotTriangulartion(fv,alpha_threshold)
    end
    meshdata.noElems = noElems;
    meshdata.dofs = dofs;
    meshdata.h_max1 = h_max1;
    meshdata.h_max2 = h_max2;
    meshdata.alpha_max = alpha_max;
    meshdata.alpha_min = alpha_min;
    meshdata.aspect = aspect;
    meshdata.waterTight = waterTight;
    meshdata.area_min = area_min;
    meshdata.area_max = area_max;
    meshdata.skewness = skewness;
    meshdata.faces = fv.faces;
    meshdata.vertices = fv.vertices;
    save(filename(1:end-4),'meshdata')
    convertSTLtoBDF(filename(1:end-4), fv)
    TR = triangulation(double(fv.faces),fv.vertices);
    stlwrite(TR,[filename(1:end-4) '_ASCII.stl'],'text')
    stlwrite(TR,[filename(1:end-4) '.stl'],'binary')
    fprintf('Finished res %d\n', res)
end
return
for res = resArr
    filename = ['C:/Users/Zetison/Dropbox/work/rhinoceros/BeTSSi/BeTSSi_mod2_BEM_res' num2str(res)];
    load(filename)
    x = meshdata.aspect;
    [~,edges] = histcounts(log10(x),100);
    histogram(x,10.^edges,'DisplayName',['res = ' num2str(res)])
    hold on
    legend show
    set(gca, 'xscale','log')
    set(gca, 'yscale','log')
%     histogram()
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    savefig('../../../../aspectRates')
    export_fig('../../../../aspectRates', '-pdf', '-transparent')
end
counter = 1;    
for res = resArr
    filename = ['C:/Users/Zetison/Dropbox/work/rhinoceros/BeTSSi/BeTSSi_mod2_BEM_res' num2str(res)];
    load(filename)
    fv.vertices = meshdata.vertices;
    fv.faces = meshdata.faces;
    [noElems, dofs, h_max1, h_max2, alpha_max, alpha_min, aspect, waterTight, area_min, area_max, skewness, i, Pi] = computeMeshData(fv);
    fprintf('\t\t${\\cal M}_{%d}^{\\textsc{COMSOL}}$ & %d & %d & %g & %g & %g & %g & %g & %g\\\\\n', ...
        counter, noElems,dofs,h_max1, h_max2, alpha_min, alpha_max, max(aspect), skewness); 
    counter = counter + 1;
end

function plotTriangulartion(fv,alpha_threshold)
    figure
    P = fv.vertices;
    tri = fv.faces;
    l1 = norm2(P(tri(:,2),:)-P(tri(:,3),:));
    l2 = norm2(P(tri(:,1),:)-P(tri(:,3),:));
    l3 = norm2(P(tri(:,1),:)-P(tri(:,2),:));
    % h_max = max(2*l1.*l2.*l3./sqrt((l1+l2+l3).*(l1+l2-l3).*(l1+l3-l2).*(l2+l3-l1)));
    % l_min = min([l1, l2, l3],[],2);
    % l_max = max([l1, l2, l3],[],2);
    % aspect = max(l_max./l_min);

    alpha1 = 180*acos((l2.^2+l3.^2-l1.^2)./(2*l2.*l3))/pi;
    alpha2 = 180*acos((l1.^2+l3.^2-l2.^2)./(2*l1.*l3))/pi;
    alpha3 = 180*acos((l1.^2+l2.^2-l3.^2)./(2*l1.*l2))/pi;
    alpha = max([alpha1, alpha2, alpha3], [], 2);
    patch(fv,'FaceVertexCData', alpha,'FaceColor','flat');
    colormap('jet')
    %          'EdgeColor',       'none',        ...

    % Add a camera light, and tone down the specular highlighting
%     camlight('headlight');
%     material('dull');
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    axis equal
    axis off
%     view([-135 35]);
    colorbar
    caxis([60,180])
%     caxis([alpha_threshold,180])
%     view([90 0]);
    testView = 4;
    switch testView
        case 1
            campos([-5.7942    4.4273    4.9686])
            camtarget([-5.3799    3.3570    4.6496])
            camup([0, 0, 1])
            camva( 10.9946)
        case 2
            campos([-16.131670968818071   1.447645480546323   4.838310268566696])
            camtarget([-16.099618160375385   1.198874199085474   4.008438867255977])
            camup([0, 0, 1])
            camva( 10.9946)
        case 3
            campos([-15.275521849434936  -3.892238630926054   6.383393316712156])
            camtarget([-15.496186693850142  -1.165099877478600   2.702521925031080])
            camup([0, 0, 1])
            camva( 10.9946)
        case 4
            campos(1.0e+02*[-0.389451895422342  -1.206668187300884   0.575243476587550])
            camtarget([-24.098167184169036  -0.744622631746219   2.602087806370932])
            camup([0, 0, 1])
            camva( 10.9946)
    end
    TR = triangulation(double(fv.faces),fv.vertices);
    hold on  
    P = incenter(TR);
    F = faceNormal(TR);  
    quiver3(P(:,1),P(:,2),P(:,3), ...
         F(:,1),F(:,2),F(:,3),0.5,'color','r');
     
    figureFullScreen(gcf)
    drawnow
        
end