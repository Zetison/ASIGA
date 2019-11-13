function plotControlPts2(varCol,dofs,labels)

controlPts = varCol.controlPts;
noTypes = numel(dofs);
cmap = [1,0,0;
        1,1,0;
        0,0,1;
        1,0,1;
        0,0,0.8;
        0,1,1;
        0,1,0];

for i = 1:noTypes
    indices = dofs{i};
    p(i) = plot3(controlPts(indices,1),controlPts(indices,2),controlPts(indices,3),'o',...
                'color',cmap(i,:),'MarkerFaceColor',cmap(i,:), 'MarkerEdgeColor', cmap(i,:));
    hold on
end
legend(p,labels)
% legend show
hold off
