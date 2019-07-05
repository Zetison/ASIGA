function p_h = calculateScatteredPressureRT(varCol, P_far, plotFarField)

if ~strcmp(varCol.BC, 'SHBC')
    error('This is not implemented')
end
k = varCol.k;
Eps = 1e6*eps;

beamEnergy = varCol.beamEnergy;

p_h = zeros(size(P_far,1),numel(k));
beams = varCol.beams;
beami = beams(:,1);
O = varCol.O;
P_inc = varCol.P_inc;
d_vec = varCol.d_vec;
p_inc = @(v) P_inc*exp(1i*dot3(v,d_vec)*k);
% A = varCol.A;
d = varCol.d;
d_n = d(beami,:);
x_n = O(beami,:,end-1);

if ~plotFarField
    if plotFarField
        P_far = P_far*1e15;
    end
%     for i = 1:size(P_far,1)
    parfor i = 1:size(P_far,1)
        x = P_far(i,:);
        t = dot(repmat(x,size(x_n,1),1)-x_n,d_n,2);
        indices = t > 0;

        triangleFound = false;
        if any(indices)
            d1 = d_n(indices,:);
            x1 = x_n(indices,:);
            P1 = x1 + t(indices,[1,1,1]).*d1;
            xRep = repmat(x,size(P1,1),1);

            for j = 1:6 % loop over triangles in beam
                I2 = beams(indices,j+1);
                x2 = O(I2,:,end-1);
                d2 = d(I2,:);
                tm = dot(P1-x2,d1,2)./dot(d2,d1,2);
                P2 = x2 + tm(:,[1,1,1]).*d2;
                if j == 6
                    I3 = beams(indices,2);
                else
                    I3 = beams(indices,j+2);
                end
                x3 = O(I3,:,end-1);
                d3 = d(I3,:);
                tm = dot(P1-x3,d1,2)./dot(d3,d1,2);
                P3 = x3 + tm(:,[1,1,1]).*d3;

                v1 = P2-P1;
                v1s = sum(v1.^2,2);
                v2 = P3-P1;
                v2s = sum(v2.^2,2);
                xmP1 = xRep-P1;
                xmP1dv1 = dot(xmP1,v1,2);
                xmP1dv2 = dot(xmP1,v2,2);
                v1dv2 = dot(v1,v2,2);
                D = v1s.*v2s-v1dv2.^2;
                xi2 = (v2s.*xmP1dv1-v1dv2.*xmP1dv2)./D;
                xi3 = (v1s.*xmP1dv2-v1dv2.*xmP1dv1)./D;
                hit = and(and(xi2 >= -Eps, xi3 >= -Eps), xi2+xi3 < 1-Eps);
                if any(hit)
                    triangleFound = true;
                    break
                end
            end
            if triangleFound
                temp = find(indices);
                indicesm = temp(hit);
                triIndices = [beams(indicesm,1), I2(hit), I3(hit)];
                indicesm2 = find(beams(:,1) == I2(hit));
                indicesm3 = find(beams(:,1) == I3(hit));
                beamIndices = [indicesm, indicesm2, indicesm3];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             close all
    %             hold on
    %             axis equal
    %             [X,Y,Z] = sphere(1000);
    %             surf(X,Y,Z, 'FaceColor', 1.5*[44 77 32]/255,'LineStyle','none')
    %             view(125,14)
    %             camlight
    %             [X,Y,Z] = sphere(1000);
    %             surf(0.02*X+x(1),0.02*Y+x(2),0.02*Z+x(3), 'FaceColor', [1,0,0],'LineStyle','none')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j = 1:numel(beamIndices)
                    indicesm = beamIndices(j);

                    d1m = d(triIndices(1),:);
                    x1m = O(triIndices(1),:,end-1);
                    P1m = x1m + t(indicesm,[1,1,1]).*d1m;

                    I1m = beams(indicesm,1);
                    x1m = O(I1m,:,end-1);
                    beamArea = 0; 
                    for jj = 1:6 % loop over triangles in beam
                        I2m = beams(indicesm,jj+1);
                        x2m = O(I2m,:,end-1);
                        d2m = d(I2m,:);
                        tm = dot(P1m-x2m,d1m)/dot(d2m,d1m);

                        P2m = x2m + tm*d2m;
                        if jj == 6
                            I3m = beams(indicesm,2);
                        else
                            I3m = beams(indicesm,jj+2);
                        end
                        x3m = O(I3m,:,end-1);
                        d3m = d(I3m,:);
                        tm = dot(P1m-x3m,d1m)/dot(d3m,d1m);
                        P3m = x3m + tm*d3m;
                        triArea = norm(cross(P2m-P1m,P3m-P1m)); % (multiplied by 2)
                        beamArea = beamArea + triArea;
    %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     X = [P1m(1), P1m(1);
    %                          P2m(1), P3m(1)];
    %                     Y = [P1m(2), P1m(2);
    %                          P2m(2), P3m(2)];
    %                     Z = [P1m(3), P1m(3);
    %                          P2m(3), P3m(3)];
    % 
    %                     surf(X,Y,Z, 'FaceColor', [0,70,147]/255,'FaceAlpha',1, 'FaceLighting', 'none')
    %                     X = [x1m(1), P1m(1)];
    %                     Y = [x1m(2), P1m(2)];
    %                     Z = [x1m(3), P1m(3)];
    %                     plot3(X,Y,Z,'black')
    %                     X = [x2m(1), P2m(1)];
    %                     Y = [x2m(2), P2m(2)];
    %                     Z = [x2m(3), P2m(3)];
    %                     plot3(X,Y,Z,'black')
    %                     X = [x3m(1), P3m(1)];
    %                     Y = [x3m(2), P3m(2)];
    %                     Z = [x3m(3), P3m(3)];
    %                     plot3(X,Y,Z,'black')
    % 
    % %                     X = [varCol.o(I1m,1); reshape(O(I1m,1,2:end-1),size(O,3)-2,1)];
    % %                     Y = [varCol.o(I1m,2); reshape(O(I1m,2,2:end-1),size(O,3)-2,1)];
    % %                     Z = [varCol.o(I1m,3); reshape(O(I1m,3,2:end-1),size(O,3)-2,1)];
    % %                     plot3(X(:),Y(:),Z(:),'red')
    %                     axis off
    %                     axis equal
    % %                     export_fig(['../graphics/S1/beam2'], '-png', '-transparent', '-r300')
    %                     keyboard
    %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    beamArea = beamArea/2;
                    d2m = d(triIndices(2),:);
                    x2m = O(triIndices(2),:,end-1);
                    tm = dot(P1m-x2m,d1m)/dot(d2m,d1m);
                    P2m = x2m + tm*d2m;
                    d3m = d(triIndices(3),:);
                    x3m = O(triIndices(3),:,end-1);
                    tm = dot(P1m-x3m,d1m)/dot(d3m,d1m);
                    P3m = x3m + tm*d3m;
                    triArea = norm(cross(P2m-P1m,P3m-P1m)); % (multiplied by 2)
                    xi_1 = norm(cross(P2m-x,P3m-x))/triArea;
                    if plotFarField
                        nrmx = norm(x);
                        expFact = nrmx*exp(-1i*k*dot(x/nrmx,x1m));
                    else
                        R = norm(x-x1m);
                        expFact = exp(1i*k*R);
                    end
        %                 p_h(i,:) = p_h(i,:) + sqrt(beamEnergy/beamArea)*xi_1*A(I1(hit),:).*expFact;
                    p_h(i,:) = p_h(i,:) + sqrt(beamEnergy/beamArea)*xi_1*p_inc(x1m).*expFact;
                    triIndices = [triIndices(2:end), triIndices(1)];
                    if isempty(indicesm2)
                        triIndices = [triIndices(2:end), triIndices(1)];
                    end
                end
            end

        end
    end
else
    temp = norm2(P_far);
    P_far = P_far./temp(:,[1,1,1]);
%     for i = 1:size(P_far,1)
    parfor i = 1:size(P_far,1)
        x = P_far(i,:);
        t = dot(repmat(x,size(x_n,1),1),d_n,2);
        indices = t > 0;

        if any(indices)
            d1 = d_n(indices,:);
            xRep = repmat(x,size(d1,1),1);
            hit = zeros(size(d1,1),6);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             close all
%             hold on
%             axis equal
%             [X,Y,Z] = sphere(1000);
%             surf(X,Y,Z, 'FaceColor', 1.5*[44 77 32]/255,'LineStyle','none')
%             view(125,14)
%             camlight
%             [X,Y,Z] = sphere(1000);
% %             surf(0.02*X+x(1),0.02*Y+x(2),0.02*Z+x(3), 'FaceColor', [1,0,0],'LineStyle','none')
%             plot3(d(:,1),d(:,2),d(:,3),'*','color','blue')
%             plot3(x(1),x(2),x(3),'*','color','red')
%             for j = 1:size(d,1)
%                 text(d(j,1),d(j,2),d(j,3),num2str(j))
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1:6 % loop over triangles in beam
                I2 = beams(indices,j+1);
                d2 = d(I2,:);
                if j == 6
                    I3 = beams(indices,2);
                else
                    I3 = beams(indices,j+2);
                end
                d3 = d(I3,:);

                d2d1 = dot(d2,d1,2);
                d3d1 = dot(d3,d1,2);
                d2d3 = dot(d2,d3,2);
                v1dv2 = d2d3./d2d1./d3d1-1;
                xd1 = dot(xRep,d1,2);
                xd2 = dot(xRep,d2,2);
                xd3 = dot(xRep,d3,2);
                v1s = 1./(d2d1).^2 - 1;
                v2s = 1./(d3d1).^2 - 1;
                D = v1s.*v2s - v1dv2.^2;
                xi2 = (v2s.*(xd2./xd1./d2d1-1) - v1dv2.*(xd3./xd1./d3d1-1))./D;
                xi3 = (v1s.*(xd3./xd1./d3d1-1) - v1dv2.*(xd2./xd1./d2d1-1))./D;
                hit(:,j) = and(and(xi2 >= -Eps, xi3 >= -Eps), xi2+xi3 < 1+Eps)*(j+1);
                
%                 
%                 triArea = norm(cross(d2m/d2d1,d3m/d3d1) + cross(d1m,d2m/d1d2-d3m/d1d3));
%                 xd1 = dot(x,d1m);
%                 xi1m = norm(xd1^2*cross(d2m/dot(d2m,d1m),d3m/dot(d3m,d1m)) - xd1/dot(d3m,d1m)*cross(x,d3m) - xd1/dot(d2m,d1m)*cross(d2m,x))/triArea;
            end
            hits = any(hit,2);
            if any(hits)
                temp = find(indices);
                beamIndices = temp(hits);
                hitsi = find(hits);
                beamIndices(:,2) = 0;
                for j = 1:size(beamIndices,1)
                    beamIndices(j,2) = hit(hitsi(j),find(hit(hitsi(j),:),1));
                end                
                for j = 1:size(beamIndices,1)
                    indicesm = beamIndices(j,1);
                    idx = beamIndices(j,2);
                    if idx == 7
                        idx2 = 2;
                    else
                        idx2 = idx+1;
                    end
                    triIndices = [beams(indicesm,1), beams(indicesm,idx), beams(indicesm,idx2)];
                    d1m = d(triIndices(1),:);
                    x1m = O(triIndices(1),:,end-1);
                    
                    beamArea = 0; 
                    for jj = 1:6 % loop over triangles in beam
                        I2m = beams(indicesm,jj+1);
                        d2m = d(I2m,:);

                        if jj == 6
                            I3m = beams(indicesm,2);
                        else
                            I3m = beams(indicesm,jj+2);
                        end
                        d3m = d(I3m,:);
                        triArea = norm(cross(d2m/dot(d1m,d2m),d3m/dot(d1m,d3m)) + cross(d1m,d2m/dot(d1m,d2m)-d3m/dot(d1m,d3m)));
                        beamArea = beamArea + triArea;
                    end
                    xd1 = dot(x,d1m);
                    beamArea = xd1^2*beamArea/2;
                    d2m = d(triIndices(2),:);
                    d3m = d(triIndices(3),:);
                    xi1m = norm(cross(d2m/dot(d1m,d2m),d3m/dot(d1m,d3m)) + cross(x/xd1,d2m/dot(d1m,d2m)-d3m/dot(d1m,d3m))) ...
                          /norm(cross(d2m/dot(d1m,d2m),d3m/dot(d1m,d3m)) + cross(d1m,  d2m/dot(d1m,d2m)-d3m/dot(d1m,d3m)));
                    
        %                 p_h(i,:) = p_h(i,:) + sqrt(beamEnergy/beamArea)*xi_1*A(I1(hit),:).*expFact;
                    p_h(i,:) = p_h(i,:) + sqrt(beamEnergy/beamArea)*xi1m*p_inc(x1m).*exp(-1i*k*dot(x,x1m));
                end
            end
        end
    end
end