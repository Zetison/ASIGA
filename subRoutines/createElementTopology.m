patchTop = varCol.patchTop;
counter = 1;
for patch = 1:numel(patches)
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    noElementsXi = length(uniqueXi)-1;
    noElementsEta = length(uniqueEta)-1;
    eMap = zeros(noElementsXi,noElementsEta);
    for i_eta = 1:noElementsEta
        for i_xi = 1:noElementsXi
            eMap(i_xi,i_eta) = counter;
            counter = counter + 1;
        end
    end
    patches{patch}.eMap = eMap;
end

for patch = 1:numel(patches)
    eMap = patches{patch}.eMap;

    eNeighbourPatch = NaN(size(eMap,1),size(eMap,2),4);
    eTwistPatch = NaN(size(eMap,1),size(eMap,2),4);

    eNeighbourPatch(1:end-1,:,1) = eMap(2:end,:); % east
    eNeighbourPatch(:,1:end-1,2) = eMap(:,2:end); % north
    eNeighbourPatch(2:end,:,3) = eMap(1:end-1,:); % west
    eNeighbourPatch(:,2:end,4) = eMap(:,1:end-1); % south
    eTwistPatch(1:end-1,:,1) = 0; % east
    eTwistPatch(:,1:end-1,2) = 0; % north
    eTwistPatch(2:end,:,3) = 0; % west
    eTwistPatch(:,2:end,4) = 0; % south
    
    Dir = 1; % east
    twist = patchTop{patch}(Dir,2);
    neighbourPatchIdx = patchTop{patch}(Dir,1);
    eTwistPatch(end,:,Dir) = twist;
    switch twist
        case 0
            eNeighbourPatch(end,:,Dir) = patches{neighbourPatchIdx}.eMap(1,:);
        case 1
            eNeighbourPatch(end,:,Dir) = patches{neighbourPatchIdx}.eMap(:,end);
        case 2
            eNeighbourPatch(end,:,Dir) = patches{neighbourPatchIdx}.eMap(end,end:-1:1);
        case 3
            eNeighbourPatch(end,:,Dir) = patches{neighbourPatchIdx}.eMap(end:-1:1,1);
    end
    Dir = 2; % north
    twist = patchTop{patch}(Dir,2);
    neighbourPatchIdx = patchTop{patch}(Dir,1);
    eTwistPatch(:,end,Dir) = twist;
    switch twist
        case 0
            eNeighbourPatch(:,end,Dir) = patches{neighbourPatchIdx}.eMap(:,1);
        case 1
            eNeighbourPatch(:,end,Dir) = patches{neighbourPatchIdx}.eMap(1,end:-1:1);
        case 2
            eNeighbourPatch(:,end,Dir) = patches{neighbourPatchIdx}.eMap(end:-1:1,end);
        case 3
            eNeighbourPatch(:,end,Dir) = patches{neighbourPatchIdx}.eMap(end,:);
    end
    Dir = 3; % west
    twist = patchTop{patch}(Dir,2);
    neighbourPatchIdx = patchTop{patch}(Dir,1);
    eTwistPatch(1,:,Dir) = twist;
    switch twist
        case 0
            eNeighbourPatch(1,:,Dir) = patches{neighbourPatchIdx}.eMap(end,:);
        case 1
            eNeighbourPatch(1,:,Dir) = patches{neighbourPatchIdx}.eMap(:,1);
        case 2
            eNeighbourPatch(1,:,Dir) = patches{neighbourPatchIdx}.eMap(1,end:-1:1);
        case 3
            eNeighbourPatch(1,:,Dir) = patches{neighbourPatchIdx}.eMap(end:-1:1,end);
    end
    Dir = 4; % south
    twist = patchTop{patch}(Dir,2);
    neighbourPatchIdx = patchTop{patch}(Dir,1);
    eTwistPatch(:,1,Dir) = twist;
    switch twist
        case 0
            eNeighbourPatch(:,1,Dir) = patches{neighbourPatchIdx}.eMap(:,end);
        case 1
            eNeighbourPatch(:,1,Dir) = patches{neighbourPatchIdx}.eMap(end,end:-1:1);
        case 2
            eNeighbourPatch(:,1,Dir) = patches{neighbourPatchIdx}.eMap(end:-1:1,1);
        case 3
            eNeighbourPatch(:,1,Dir) = patches{neighbourPatchIdx}.eMap(1,:);
    end
    patches{patch}.eNeighbourPatch = eNeighbourPatch;
    patches{patch}.eTwistPatch = eTwistPatch;
end

counter = 1;
eNeighbour = NaN(noElems,4,2);
for patch = 1:numel(patches)
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    noElementsXi = length(uniqueXi)-1;
    noElementsEta = length(uniqueEta)-1;
    eNeighbourPatch = patches{patch}.eNeighbourPatch;
    eTwistPatch = patches{patch}.eTwistPatch;
    for i_eta = 1:noElementsEta
        for i_xi = 1:noElementsXi
            eNeighbour(counter,:,1) = eNeighbourPatch(i_xi,i_eta,:);
            eNeighbour(counter,:,2) = eTwistPatch(i_xi,i_eta,:);
            counter = counter + 1;
        end
    end
end