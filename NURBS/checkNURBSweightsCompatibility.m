function equalWeights = checkNURBSweightsCompatibility(nurbs)


varCol.dimension = 1;
varCol.nurbs = nurbs;
varCol = findDofsToRemove(generateIGAmesh(convertNURBS(varCol)));
weights = varCol.weights;
controlPts = varCol.controlPts;

Eps = 1e-10;
equalWeights = true;
indices = [];
for i = 1:numel(varCol.gluedNodes)
    I = varCol.gluedNodes{i}(1);
    for j = 2:numel(varCol.gluedNodes{i})
        J = varCol.gluedNodes{i}(j);
        if abs(weights(I) - weights(J)) > Eps
            equalWeights = false;
            indices = [indices,I];
            break
        end
    end
end
plot3(controlPts(indices,1),controlPts(indices,2),controlPts(indices,3),'o','color','black','MarkerFaceColor', 'blue', 'markersize',10, ...
                            'MarkerEdgeColor', 'black');
hold on
    