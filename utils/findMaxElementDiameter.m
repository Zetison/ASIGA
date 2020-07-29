function [hmax, hmin, diagsMax] = findMaxElementDiameter(patches)

if ~iscell(patches)
    patches = {patches};
end
H_max = zeros(numel(patches),1);
H_min = zeros(numel(patches),1);
diagsMax = cell(numel(patches),1);
parfor patch = 1:numel(patches)
    h_max = 0;
    h_min = inf;
    nurbs = patches{patch}.nurbs;
    switch nurbs.d_p
        case 3
            Xi = nurbs.knots{1};
            Eta = nurbs.knots{2};
            Zeta = nurbs.knots{3};

            extraXiPts = 0;
            extraEtaPts = 0;
            extraZetaPts = 0;

            noElems = (numel(unique(Xi))-1)*(numel(unique(Eta))-1)*(numel(unique(Zeta))-1);

            [nodes, ~, visElements] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, nurbs);


            for e = 1:noElems
                ind = visElements(e,:);

                diameter = norm(nodes(ind(1),1:3) - nodes(ind(7),1:3));
                if h_max < diameter
                    h_max = diameter;
                end
                if h_min > diameter
                    h_min = diameter;
                end

                diameter = norm(nodes(ind(2),1:3) - nodes(ind(8),1:3));
                if h_max < diameter
                    h_max = diameter;
                end
                if h_min > diameter
                    h_min = diameter;
                end

                diameter = norm(nodes(ind(3),1:3) - nodes(ind(5),1:3));
                if h_max < diameter
                    h_max = diameter;
                end
                if h_min > diameter
                    h_min = diameter;
                end

                diameter = norm(nodes(ind(4),1:3) - nodes(ind(6),1:3));
                if h_max < diameter
                    h_max = diameter;
                end
                if h_min > diameter
                    h_min = diameter;
                end
            end
            diaMax = NaN;
        case 2
            Xi = nurbs.knots{1};
            Eta = nurbs.knots{2};

            extraXiPts = 0;
            extraEtaPts = 0;

            [nodes, ~, visElements] = buildVisualization3dSurfaceMesh(Xi, Eta, extraXiPts, extraEtaPts, nurbs);

            dia = [norm2(nodes(visElements(:,1),:) - nodes(visElements(:,3),:)), ...
                   norm2(nodes(visElements(:,2),:) - nodes(visElements(:,4),:))];
            diaMax = max(dia,[],2);
            diaMin = max(dia,[],2);
            h_max = max(diaMax);
            h_min = min(diaMin);
    end
    H_max(patch) = h_max;
    H_min(patch) = h_min;
    diagsMax{patch} = diaMax;
end
diagsMax = cell2mat(diagsMax);
hmin = min(H_min);
hmax = max(H_max);