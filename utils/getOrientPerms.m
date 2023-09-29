function [flips,indices,orientMap,noOrientations] = getOrientPerms(d_p)
% Get arrays for all different permutations of a nurbs of dimension d_p
flips = cell(2^d_p,1);
indices = flipud(perms(1:d_p));

if d_p >= 1
    flips{1} = [];
    flips{2} = 1;
end

if d_p >= 2
    flips{3} = 2;
    flips{4} = [1,2];
end

if d_p >= 3
    flips{5} = 3;
    flips{6} = [1,3];
    flips{7} = [2,3]; 
    flips{8} = [1,2,3];
end

if d_p >= 4
    error('not implemented')
end

if nargout > 2
    orientMap = cell(numel(flips)*size(indices,1),3);
    counter = 1;
    for i = 1:size(indices,1)
        for j = 1:numel(flips)
            orientMap{counter,1} = flips{j};
            orientMap{counter,2} = indices(i,:);
            counter = counter + 1;
        end
    end
    noOrientations = size(orientMap,1);
end
