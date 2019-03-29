function [elRange,elConn] = buildConnectivity(p,Xi,noElems)

elRange         = zeros(noElems,2);   
elKnotIndices   = zeros(noElems,2);
elConn          = zeros(noElems,p+1);
ne              = length(unique(Xi))-1;   % number of elements

element         = 1;
previousKnotVal = 0;

for i = 1:length(Xi) 
    currentKnotVal = Xi(i);
    if Xi(i)~=previousKnotVal
        elRange(element,:)=[previousKnotVal currentKnotVal];
        elKnotIndices(element,:)=[i-1 i];
        element=element+1;
    end
    previousKnotVal = currentKnotVal;
end

numRepeatedKnots=0;

for e = 1:ne
    indices=(elKnotIndices(e,1)-p+1):elKnotIndices(e,1);
    previousKnotVals=Xi(indices);
    currentKnotVals=ones(1,p)*Xi(elKnotIndices(e,1));
    if isequal(previousKnotVals,currentKnotVals) && length(nonzeros(previousKnotVals))>1;
        numRepeatedKnots=numRepeatedKnots+1;
    end
    elConn(e,:) = (elKnotIndices(e,1)-p):elKnotIndices(e,1);
end
