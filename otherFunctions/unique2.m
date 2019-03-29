function [a,I] = unique2(a)

m = size(a,1);
        
% groupsSortA indicates the location of non-matching entries.
groupsSortA = a(1:m-1,:) ~= a(2:m,:);
groupsSortA = any(groupsSortA,2);

groupsSortA = [true; groupsSortA];          % First row is always a member of unique list.
            
a = a(groupsSortA,:);                                   % Create unique list by indexing into unsorted a.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indSortC = uint32((1:m)');

[~,ind] = sort(a(:,2));
indSortC = indSortC(ind);
[~,ind] = sort(a(indSortC,1));
indSortC = indSortC(ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthGroupsSortA = diff(find([groupsSortA; true]));    % Determine how many of each of the above indices there are in IC.

diffIndSortC = diff(indSortC);
diffIndSortC = [indSortC(1); diffIndSortC];
clear indSortC

indLengthGroupsSortA = cumsum([1; lengthGroupsSortA]);  % Get the correct amount of each index.
indLengthGroupsSortA(end) = [];

I(indLengthGroupsSortA,1) = diffIndSortC;        % Since indCOrderedBySortA is not already established as a column,

if sum(lengthGroupsSortA) ~= length(I)
    I(sum(lengthGroupsSortA),1) = 0;
end

I = cumsum(I);
