
counter = 1;
for patch = 1:numel(patches)
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    noElementsXi = length(uniqueXi)-1;
    noElementsEta = length(uniqueEta)-1;
    eMap = zeros(noElementsXi,noElementsEta);
    for i_xi = 1:noElementsXi
        for i_eta = 1:noElementsEta
            eMap(i_xi,i_eta) = counter;
            counter = counter + 1;
        end
    end
    patches{patch}.eMap = eMap;
end

for patch = 1:numel(patches)
    eMap = patches{patch}.eMap;

    eNeighbourPatch = NaN(size(eMap,1),size(eMap,2),4);

    eNeighbourPatch(1:end-1,:,1) = eMap(2:end,:); % east
    eNeighbourPatch(:,1:end-1,2) = eMap(:,2:end); % north
    eNeighbourPatch(2:end,:,3) = eMap(1:end-1,:); % west
    eNeighbourPatch(:,2:end,4) = eMap(:,1:end-1); % south
    
    eNeighbourPatch(end,:,1) = patches{patchTop(1)}.eMap;
    patches{patch}.eNeighbourPatch = eNeighbourPatch;
end

counter = 1;
eNeighbour = NaN(noElems,4);
for patch = 1:numel(patches)
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    noElementsXi = length(uniqueXi)-1;
    noElementsEta = length(uniqueEta)-1;
    eNeighbourPatch = patches{patch}.eNeighbourPatch;
    for i_xi = 1:noElementsXi
        for i_eta = 1:noElementsEta
            eNeighbour(counter,:) = eNeighbourPatch(i_xi,i_eta,:);
            counter = counter + 1;
        end
    end
end