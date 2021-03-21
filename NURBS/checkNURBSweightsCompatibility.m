function checkNURBSweightsCompatibility(task)

for i_domain = 1:numel(task.varCol)
    weights = task.varCol{i_domain}.weights;
    controlPts = task.varCol{i_domain}.controlPts;

    Eps = 1e-10;
    equalWeights = true;
    indices = [];
    for i = 1:numel(task.varCol{i_domain}.gluedNodes)
        I = task.varCol{i_domain}.gluedNodes{i}(1);
        for j = 2:numel(task.varCol{i_domain}.gluedNodes{i})
            J = task.varCol{i_domain}.gluedNodes{i}(j);
            if abs(weights(I) - weights(J)) > Eps
                equalWeights = false;
                indices = [indices,I];
                break
            end
        end
    end
    if task.prePlot.plot3Dgeometry
        set(0, 'CurrentFigure', task.prePlot.fig3Dplot)
        plot3(controlPts(indices,1),controlPts(indices,2),controlPts(indices,3),'o','color','black','MarkerFaceColor', 'blue', 'markersize',10, ...
                                    'MarkerEdgeColor', 'black');
        hold on
    end
    if ~equalWeights
        warning('NURBS:weights','Some weights in the geometry are not equal. For geometries containing singularities this might be ok (this warning may then be supressed using the key NURBS:weights).')
        % supress the following warning with warning('off','NURBS:weights')
        % in your getTask_<model> script if the model contains
        % singularities
    end
end