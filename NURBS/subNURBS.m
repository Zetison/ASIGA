function subnurbs = subNURBS(varargin)
options = struct('outwardPointingNormals',false,...
                 'inwardPointingNormals',false,...
                 'useFlipping',true);
nurbs = varargin{1};
if ~iscell(nurbs)
    nurbs = {nurbs};
end
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
noPatches = numel(nurbs);
if isfield(options,'at') && ~iscell(options.at)
    temp = options.at;
    options.at = cell(1,noPatches);
    [options.at(:)] = deal({temp});
end
if options.outwardPointingNormals && options.inwardPointingNormals
    error('Contradicting options')
end
subnurbs = cell(1,6*noPatches);
counter = 1;
for patch = 1:noPatches
    d = nurbs{patch}.d;
    d_p = nurbs{patch}.d_p;
    coeffs = nurbs{patch}.coeffs;
    dimensions = size(coeffs);
    number = nurbs{patch}.number;
    if isfield(options,'at')
        at = options.at{patch};
        if numel(at) == 1
            atv = at;
            at = false(d_p,2).';
            at(atv) = true;
            at = at.';
        end
    else
        at = true(d_p,2);
    end
    for i = 1:d_p
        for j = 1:2
            if at(i,j)
                % Get the index J of the subnurbs of interest
                if j == 1
                    J = 1;
                else
                    J = dimensions(i+1);
                end
                idx = mod(i:i+1+d_p-2,d_p)+1; % cycle index i to dimension d_p
                % permute to obtains normal vectors aligned with the leftover parametric direction
                controlPts = permute(slc(coeffs,J,i+1),[1,idx+1]); % permute the index not part of the subnurbs to the back
                idx = idx(1:end-1); % Skip the index not part of the subnurbs
                if d_p > 1
                    controlPts = reshape(controlPts,[d+1,number(idx)]);
                end
                knots = nurbs{patch}.knots(idx);
                subnurbs(counter) = createNURBSobject(controlPts,knots);
                if isfield(nurbs{patch},'isPML')
                    subnurbs{counter}.isPML = nurbs{patch}.isPML;
                end
                if options.outwardPointingNormals && j == 1
                    if options.useFlipping
                        subnurbs(counter) = flipNURBSparametrization(subnurbs(counter));
                    else
                        subnurbs(counter) = permuteNURBS(subnurbs(counter),[2,1]);
                    end
                end
                if options.inwardPointingNormals && j == 2
                    if options.useFlipping
                        subnurbs(counter) = flipNURBSparametrization(subnurbs(counter));
                    else
                        subnurbs(counter) = permuteNURBS(subnurbs(counter),[2,1]);
                    end
                end
                counter = counter + 1;
            end
        end
    end
end
subnurbs(counter:end) = [];
