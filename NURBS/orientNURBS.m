function nurbs = orientNURBS(nurbs,orient,useInverse)

% Orient the parametrization of nurbs as follows for a d_p=2 patch
% orient=0: no change
% orient=1: flip first parametric direction                                     (this flips the normal vector)
% orient=2: flip second parametric direction                                    (this flips the normal vector)
% orient=3: flip both parametric directions
% orient=4:                                      swaps parametric directions    (this flips the normal vector)
% orient=5: flip first parametric direction  and swaps parametric directions     
% orient=6: flip second parametric direction and swaps parametric directions                                    
% orient=7: flip both parametric directions  and swaps parametric directions    (this flips the normal vector)


if orient == 0
    return
end
if nargin < 3
    useInverse = false;
end
if useInverse
    if orient == 5
        orient = 6;
    elseif orient == 6
        orient = 5;
    end % Otherwise the orient mapping its own inverse nurbs == orientNURBS(orientNURBS(nurbs,orient,~useInverse),orient,useInverse)
end

if isfield(nurbs,'d_p')
    d_p = nurbs.d_p;
    [~,~,orientMap] = getOrientPerms(d_p);
    
    
    flips = orientMap{orient+1,1};
    indices = orientMap{orient+1,2};
    
    for jj = 1:numel(flips)
        flipIdx = flips(jj);
        nurbs.coeffs = flip(nurbs.coeffs,flipIdx+1);
        nurbs.knots{flipIdx} = 1-flip(nurbs.knots{flipIdx});
    end
    if d_p > 0
        nurbs.coeffs = permute(nurbs.coeffs,[1,indices+1]);
    end
    nurbs.knots = nurbs.knots(indices);
    nurbs.degree = nurbs.degree(indices);
    nurbs.number = nurbs.number(indices);
else % "nurbs" is a cell of knots
    d_p = numel(nurbs);
    [~,~,orientMap] = getOrientPerms(d_p);
    
    
    flips = orientMap{orient+1,1};
    indices = orientMap{orient+1,2};
    
    for jj = 1:numel(flips)
        flipIdx = flips(jj);
        nurbs{flipIdx} = 1-flip(nurbs{flipIdx});
    end
    nurbs = nurbs(indices);
end
