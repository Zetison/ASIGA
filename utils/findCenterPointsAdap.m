function [centerPts,subElementMap] = findCenterPointsAdap(patches,pIndex,noElems,index,elRangeXi,elRangeEta,level,Eps)

maxSubElements = (4^(level+1)-1)/3;
centerPts = cell(noElems,1);
for e = 1:noElems
    centerPts{e} = struct('elRangeXi_sub',repmat({NaN(1,2)}, maxSubElements, 1),'elRangeEta_sub',repmat({NaN(1,2)}, maxSubElements, 1),...
                                        'h',repmat({NaN(1,1)}, maxSubElements, 1),'x_5',repmat({NaN(1,3)}, maxSubElements, 1));
end
elRange_t_sub = struct('elRangeXi_t_sub',repmat({NaN(1,2)}, maxSubElements, 1),'elRangeEta_t_sub',repmat({NaN(1,2)}, maxSubElements, 1));
if ~iscell(patches)
    patches = {patches};
end

subElementMap = zeros(maxSubElements,4);
counter = 1;
e_sub = 2;
for i = 0:level
    noSubElems = 2^i;
    for i_eta = 0:noSubElems-1
        for i_xi = 0:noSubElems-1
            subElementMap(counter,:) = e_sub:e_sub+3;
            counter = counter + 1;
            e_sub = e_sub+4;
        end
    end
end
temp = zeros(maxSubElements,5,2);
[temp,elRange_t_sub] = nestedSubElements(temp,elRange_t_sub,subElementMap,1,0,level,[-1,1],[-1,1]);
temp = reshape(temp,maxSubElements*5,2);
[temp, ~, IC] = uniquetol(temp,Eps,'ByRows',true);

for e = 1:noElems  
    patch_y = pIndex(e); % New

    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    
    Xi_e = [Xi_e(1)+eps,Xi_e(2)-eps];
    Eta_e = [Eta_e(1)+eps,Eta_e(2)-eps];
    
%     centerPts{e}(1) = struct('elRangeXi_sub',Xi_e,'elRangeEta_sub',Eta_e,'h',[],'x_5',[]);
%     temp(1,1,:) = [Xi_e(1),Eta_e(1)];
%     temp(1,2,:) = [Xi_e(2),Eta_e(1)];
%     temp(1,3,:) = [Xi_e(1),Eta_e(2)];
%     temp(1,4,:) = [Xi_e(2),Eta_e(2)];
%     temp(1,5,:) = [mean(Xi_e),mean(Eta_e)];
%     e_sub = 2;
%     for i = 0:level-1
%         noSubElems = 2^i;
%         for i_eta = 0:noSubElems-1
%             for i_xi = 0:noSubElems-1
%                 counter = e_sub;
%                 for i_eta2 = 0:1
%                     Eta_e_sub = parent2ParametricSpace(Eta_e, -1+2*((2*i_eta+i_eta2):(2*i_eta+i_eta2+1))/(2*noSubElems));
%                     for i_xi2 = 0:1
%                         Xi_e_sub = parent2ParametricSpace(Xi_e, -1+2*((2*i_xi+i_xi2):(2*i_xi+i_xi2+1))/(2*noSubElems));
%                         centerPts{e}(counter).elRangeXi_sub = Xi_e_sub;
%                         centerPts{e}(counter).elRangeEta_sub = Eta_e_sub;
%                         temp(counter,1,:) = [Xi_e_sub(1),Eta_e_sub(1)];
%                         temp(counter,2,:) = [Xi_e_sub(2),Eta_e_sub(1)];
%                         temp(counter,3,:) = [Xi_e_sub(1),Eta_e_sub(2)];
%                         temp(counter,4,:) = [Xi_e_sub(2),Eta_e_sub(2)];
%                         temp(counter,5,:) = [mean(Xi_e_sub),mean(Eta_e_sub)];
%                         counter = counter + 1;
%                     end
%                 end
%                 e_sub = e_sub+4;
%             end
%         end
%     end
    for e_sub = 1:maxSubElements
        centerPts{e}(e_sub).elRangeXi_sub = parent2ParametricSpace(Xi_e, elRange_t_sub(e_sub).elRangeXi_t_sub);
        centerPts{e}(e_sub).elRangeEta_sub = parent2ParametricSpace(Eta_e, elRange_t_sub(e_sub).elRangeEta_t_sub);
    end
    xi = parent2ParametricSpace(Xi_e, temp(:,1));
    eta = parent2ParametricSpace(Eta_e, temp(:,2));
    yy = evaluateNURBS_2ndDeriv(patches{patch_y}.nurbs, [xi,eta]);
    yy = permute(reshape(yy(IC,:),maxSubElements,5,3),[1,3,2]);
    for e_sub = 1:maxSubElements
        centerPts{e}(e_sub).h = max([norm(yy(e_sub,:,1)-yy(e_sub,:,4)), norm(yy(e_sub,:,2)-yy(e_sub,:,3))]);
        centerPts{e}(e_sub).x_5 = yy(e_sub,:,5);
    end
end


function [temp,elRange_t_sub] = nestedSubElements(temp,elRange_t_sub,subElementMap,e_sub,level,maxLevel,Xi_e_t_sub,Eta_e_t_sub)

elRange_t_sub(e_sub).elRangeXi_t_sub = Xi_e_t_sub;
elRange_t_sub(e_sub).elRangeEta_t_sub = Eta_e_t_sub;
temp(e_sub,1,:) = [Xi_e_t_sub(1),Eta_e_t_sub(1)];
temp(e_sub,2,:) = [Xi_e_t_sub(2),Eta_e_t_sub(1)];
temp(e_sub,3,:) = [Xi_e_t_sub(1),Eta_e_t_sub(2)];
temp(e_sub,4,:) = [Xi_e_t_sub(2),Eta_e_t_sub(2)];
temp(e_sub,5,:) = [mean(Xi_e_t_sub),mean(Eta_e_t_sub)];
if level < maxLevel
    Xi_e_tArr = [Xi_e_t_sub(1),mean(Xi_e_t_sub),Xi_e_t_sub(2)];
    Eta_e_tArr = [Eta_e_t_sub(1),mean(Eta_e_t_sub),Eta_e_t_sub(2)];
    counter = 1;
    for i_eta = 1:2
        for i_xi = 1:2
            [temp,elRange_t_sub] = nestedSubElements(temp,elRange_t_sub,subElementMap,subElementMap(e_sub,counter),...
                level+1,maxLevel,Xi_e_tArr(i_xi:i_xi+1),Eta_e_tArr(i_eta:i_eta+1));
            counter = counter + 1;
        end
    end
end