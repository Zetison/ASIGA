function plotControlPts2(varCol,dofs,labels)

controlPts = varCol.controlPts;
noTypes = numel(dofs);
cmap = [1,0,0;
        1,1,0;
        0,0,1;
        1,0,1;
        0,0,0.5;
        0,1,1;
        0,1,0];

hold on
for i = 1:noTypes
    indices = dofs{i};
    p(i) = plot3(controlPts(indices,1),controlPts(indices,2),controlPts(indices,3),'o',...
                'color',cmap(i,:),'MarkerFaceColor',cmap(i,:), 'MarkerEdgeColor', cmap(i,:));
end
legend(p,labels)
% legend show
