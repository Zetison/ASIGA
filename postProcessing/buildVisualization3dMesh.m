% function [node,elementV]=buildVisualization3dMesh(controlPts,weights,...
%                         uKnot,vKnot,wKnot,p,q,r)
%
% build a H8 mesh (eight node brick elements) from image of knots.
% this mesh is used for visualization purpose.
% Standard FE visualization function like plot_field
% can then be reused.
% Vinh Phu Nguyen
% Johns Hopkins University

% get rid of zero measure knot spans

uniqueXiVec = unique(Xi);
uniqueEtaVec = unique(Eta);
uniqueZetaVec = unique(Zeta);

% number of distinct knot values

noUniqueXiKnots = length(uniqueXiVec);
noUniqueEtaKnots = length(uniqueEtaVec);
noUniqueZetaKnots = length(uniqueZetaVec);

nodes  = zeros(noUniqueXiKnots*noUniqueEtaKnots*noUniqueZetaKnots,4);
count = 1;

for wk=1:noUniqueZetaKnots
  zeta = uniqueZetaVec(wk);
  for vk=1:noUniqueEtaKnots
      eta = uniqueEtaVec(vk);
      for uk=1:noUniqueXiKnots
          xi = uniqueXiVec(uk);

          v = evaluateNURB(solid, [xi, eta, zeta]);

          nodes(count,1) = v(1);
          nodes(count,2) = v(2);
          nodes(count,3) = v(3);
          nodes(count,4) = count;

          count = count + 1;
      end
  end
end

% build H8 elements

chan  = zeros(noUniqueXiKnots,noUniqueEtaKnots,noUniqueZetaKnots);

count = 1;

for i=1:noUniqueZetaKnots
    for j=1:noUniqueEtaKnots
        for k=1:noUniqueXiKnots
            chan(i,j,k) = count;
            count       = count + 1;
        end
    end
end

connecU = zeros(noUniqueXiKnots-1,2);
connecV = zeros(noUniqueEtaKnots-1,2);
connecW = zeros(noUniqueZetaKnots-1,2);

for i=1:size(connecU,1)
   connecU(i,:) = [i i+1];	
end

for i=1:size(connecV,1)
   connecV(i,:) = [i i+1];	
end

for i=1:size(connecW,1)
   connecW(i,:) = [i i+1];	
end

noElems  = (noUniqueXiKnots-1) * (noUniqueEtaKnots-1) * (noUniqueZetaKnots-1);
elementV = zeros(noElems,8);

e = 1;
for w=1:noUniqueZetaKnots-1
    wConn = connecW(w,:);
    for v=1:noUniqueEtaKnots-1
        vConn = connecV(v,:);
        for u=1:noUniqueXiKnots-1
            c = 1;
            uConn = connecU(u,:);
            
            for i=1:length(wConn)
                for j=1:length(vConn)
                    for k=1:length(uConn)
                        elementV(e,c) = chan(wConn(i),vConn(j),uConn(k));
                        c = c + 1;
                    end
                end
            end
            e = e + 1;
        end
    end
end

% renumbering nodes according to Jack Chessa's code

col3 = elementV(:,3);
col4 = elementV(:,4);
col7 = elementV(:,7);
col8 = elementV(:,8);

elementV(:,3) = col4;
elementV(:,4) = col3;
elementV(:,7) = col8;
elementV(:,8) = col7;





