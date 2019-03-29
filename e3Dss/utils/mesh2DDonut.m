function [tri, X] = mesh2DDonut(R_o,R_i,N)


dR = R_o-R_i;
dt = 2*pi/N;
h = sqrt(3)/2*dt*R_o;
Nr = ceil(dR/h)+1;
dr0 = dR/(Nr-1);

X = zeros(Nr*N,2);
counter = 1;
for i = 1:Nr
    r = R_o - (i-1)*dr0;
    N = round(2*pi*r/dr0);
    dt = 2*pi/N;
    theta = linspace(0,2*pi-dt,N).';
    if mod(i,2)
        theta = theta + dt/2;
    end
    X(counter:counter+N-1,:) = r*[cos(theta), sin(theta)];
    counter = counter + N;
end
X = X(1:counter-1,:);
tri = delaunay(X(:,1),X(:,2));
insideIndices = zeros(size(tri,1),1);
fixedIndices = zeros(size(tri,1),1);
counterInside = 1;
counterInside2 = 1;

for i = 1:size(tri,1)
    triInd = tri(i,:);
    x = X(triInd,:);
    isOutsideCircle = norm2(x) - R_i > 10*eps;
    if sum(isOutsideCircle) == 2
        Xind = triInd(~isOutsideCircle);
        Xinside = X(Xind,:);
        theta = atan2(Xinside(2),Xinside(1));
        X(Xind,:) = R_i*[cos(theta), sin(theta)];
        fixedIndices(counterInside2) = Xind;
        counterInside2 = counterInside2 + 1;
    elseif sum(isOutsideCircle) == 1
        Xind = triInd(~isOutsideCircle);
        Xinside = X(Xind,:);
        theta1 = atan2(Xinside(1,2),Xinside(1,1));
        theta2 = atan2(Xinside(2,2),Xinside(2,1));
        X(Xind,:) = R_i*[cos(theta1), sin(theta1);
                       cos(theta2), sin(theta2)];
        fixedIndices(counterInside2:counterInside2+1) = Xind;
        counterInside2 = counterInside2 + 2;
    elseif ~isOutsideCircle
        insideIndices(counterInside) = i;
        counterInside = counterInside + 1;
    end
end
insideIndices = insideIndices(1:counterInside-1);
fixedIndices = fixedIndices(1:counterInside2-1);
removeIndices = sort(setdiff(unique(tri(insideIndices,:)),unique(fixedIndices)));
tri(insideIndices,:) = [];
for i = length(removeIndices):-1:1
    tri(tri > removeIndices(i)) = tri(tri > removeIndices(i)) - 1;
end
X(removeIndices,:) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% triplot(tri,X(:,1),X(:,2))
% for i = 1:size(X,1)
%     text(X(i,1),X(i,2),num2str(i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% end
% for i = 1:size(tri,1)
%     x = mean(X(tri(i,:),:));
%     text(x(1),x(2),num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','color','red')
% end
% axis equal
% hold on
% theta = linspace(0,2*pi,1000);
% plot(R_o*cos(theta),R_o*sin(theta),'red')
% plot(R_i*cos(theta),R_i*sin(theta),'red')
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






