function C = findCenterPoints(patches)

if ~iscell(patches)
    patches = {patches};
end
C = cell(numel(patches),1);
for patch = 1:numel(patches)
    elRange = patches{patch}.elRange;
    xi = mean(elRange{1},2);
    eta = mean(elRange{2},2);
    xi1 = copyVector(xi.',numel(eta),1);
    eta1 = copyVector(eta.',numel(xi),2);
    C{patch} = zeros(numel(xi1),3);
    for i = 1:numel(xi1)
        C{patch}(i,:) = evaluateNURBS(patches{patch}.nurbs, [xi1(i),eta1(i)]);
    end
end
C = cell2mat(C);