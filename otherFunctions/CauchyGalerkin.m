% DATA: Polynomial degree : p [Assumed to be 3 , 5 or 7]
% DATA: Number of control points : ncp [Assumed to be even ]
% RESULTS: Cauchy-Galerkin points : cg (:) [The length of cg is ncp ]
% REQUIREMENTS: Spline Toolbox
function [cg, grev] = CauchyGalerkin(p, n, Xi)
grev = aveknt(Xi, p+1); % Compute Greville points
if ~mod(p,2)
    cg = grev;
    return
end
% cg = NaN;
% return
uniqueXi = unique(Xi);

% Compute CG points from Greville points    
cg = zeros(size(grev));
C0avg = zeros(size(grev));
idx = 1;
for i = 1:n
    if sum(Xi == grev(i)) >= p
        C0avg(idx:i) = (i+idx)/2;
        idx = i;
    end
end
idx = 1;
for i = 1:length(uniqueXi)-1
    ncpInE = 0;
    
    while grev(idx + ncpInE) - uniqueXi(i+1) < 10*eps
        ncpInE = ncpInE+1;
        if idx + ncpInE == n
            ncpInE = ncpInE+1;
            break
        end
    end
    
    if ~(Xi(end) == uniqueXi(i+1)) && abs(grev(idx + ncpInE - 1) - uniqueXi(i+1)) < 10*eps %% Check for C^-1 or if grev is a knot
        ncpInE = ncpInE - 1;
    end
    switch ncpInE
        case 1
            if p == 3
                locOfCG = (1 - 1/sqrt(3))/2; 
            elseif p == 5
                locOfCG = (1 - sqrt(225-30*sqrt(30))/15)/2; 
            elseif mod(p,2) == 1 
                locOfCG = (1 - 0.5049185675126533060793906)/2;
            else
                locOfCG = 0.5;
            end
        case 2
            if p == 3
                locOfCG = [0.147170272421024, 0.701674114600789]; 
            elseif p == 5
                locOfCG = [0.236168492138485, 0.688496784202927]; 
            elseif mod(p,2) == 1 
                locOfCG = [0.255450957692400, 0.652743121420659];
            elseif p == 2
                locOfCG = [0.25, 0.500515036863072];
            elseif p == 4
                locOfCG = [0.129227111965745, 0.5];
            end
        case 3
            if p == 3
                CG1 = 0.147170272421024;
                CG2 = 0.701674114600789;
                locOfCG = [CG1, (CG1+CG2)/2, CG2]; 
            elseif p == 5
                CG1 = 0.237476255750422;
                CG2 = 0.692231057550879;
                locOfCG = [CG1, (CG1+CG2)/2, CG2]; 
            elseif mod(p,2) == 1 
                CG1 = 0.237476255750422;
                CG2 = 0.692231057550879;
                locOfCG = [CG1, (CG1+CG2)/2, CG2]; 
            elseif p == 2
                locOfCG = [0.25, 0.500515036863072, 0.75];
            elseif p == 4
                locOfCG = [0.129397883099115, 0.5, 0.75];
            end
        otherwise
            locOfCG = linspace2(0,1,ncpInE);
    end
    spnLength = uniqueXi(i+1)-uniqueXi(i);
    if idx <= C0avg(idx)
        cg(idx:idx+ncpInE-1) = uniqueXi(i) + spnLength*locOfCG;
    else
        spnLength = uniqueXi(i+1)-uniqueXi(i);
        temp = uniqueXi(i) + spnLength*(1-fliplr(locOfCG));
        cg(idx:idx+ncpInE-1) = temp;
    end
    idx = idx + ncpInE;
end

% j = 1;
% for i = 1:n
%     xi = grev(i);
%     while xi > uniqueXi(j+1)
%         j = j+1;
%     end
%     jj = 0;
%     while grev(i + jj) <= uniqueXi(j+1)
%         jj = jj+1;
%     end
%     if i <= C0avg(i)
%         spnLength = uniqueXi(j+1)-uniqueXi(j);
%         ofst = spnLength*(1-ncgp)/2;
%         cg(i) = grev(i) + ofst;
%         if sum(Xi == grev(i)) >= p
%             cg(i) = cg(i) + spnLength*0.129927739973652;
%         end
%     else
%         ofst = spnLength*(1-ncgp)/2;
%         cg(i) = grev(i) - ofst;
%         if sum(Xi == grev(i)) >= p
%             cg(i) = cg(i) - spnLength*0.129927739973652;
%         end
%     end
% end  
%     else
%         j = 1;
%         for i = 1:n
%             xi = grev(i);
%             while xi > uniqueXi(j+1)
%                 j = j+1;
%             end
%             jj = 0;
%             while grev(i + jj) < uniqueXi(j+1)
%                 jj = jj+1;
%             end
%             
%             if i <= C0avg(i)
%                 spnLength = uniqueXi(j+1)-uniqueXi(j);
%                 ofst = spnLength*(1-ncgp)/2;
%                 cg(i) = uniqueXi(j) + ofst;
%     %             cg(i) = grev(i) + ofst;
%             else
%                 ofst = spnLength*(1-ncgp)/2;
%                 cg(i) = uniqueXi(j) - ofst;
%     %             cg(i) = grev(i) - ofst;
%             end
%         end  
%     end
% end
