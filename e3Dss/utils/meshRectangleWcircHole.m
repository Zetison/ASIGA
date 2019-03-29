function [tri, X] = meshRectangleWcircHole(xl,xu,R,N)

Lx = xu(1)-xl(1);
Ly = xu(2)-xl(2);
d = Lx/(N-1);
h = sqrt(3)/2*d;
Ny = ceil(Ly/h)+1;
x1 = linspace(xl(1),xu(1),N);
x2 = [xl(1), x1(1:end-1)+d/2, xu(1)];
x3 = [x1,x2];

X = zeros(floor(Ny/2)*length(x3)+mod(Ny,2)*N,2);
X(1:floor(Ny/2)*length(x3),1) = repmat(x3,1,floor(Ny/2));
if mod(Ny,2)
    X(end-N+1:end,1) = x1;
end

N1 = length(x1);
N2 = length(x2);
N = N1+N2;
hy = Ly/(Ny-1);
for i = 1:floor(Ny/2)
    X(N*(i-1)+1:N*(i-1)+N1,2) = (xl(2)+2*hy*(i-1))*ones(1,N1);
    X(N*(i-1)+N1+1:N*(i-1)+N1+N2,2) = (xl(2)+(2*i-1)*hy)*ones(1,N2);
end
if mod(Ny,2)
    X(end-N1+1:end,2) = xu(2)*ones(N1,1);
end
tri = delaunay(X(:,1),X(:,2));
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
% plot(R*cos(theta),R*sin(theta),'red')
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
insideIndices = zeros(size(tri,1),1);
fixedIndices = zeros(size(tri,1),1);
counterInside = 1;
counterInside2 = 1;
% figure(2)
for i = 1:size(tri,1)
    triInd = tri(i,:);
    x = X(triInd,:);
    isOutsideCircle = norm2(x) - R > 10*eps;
    if sum(isOutsideCircle) == 2
        Xind = triInd(~isOutsideCircle);
        Xinside = X(Xind,:);
        theta = atan2(Xinside(2),Xinside(1));
        X(Xind,:) = R*[cos(theta), sin(theta)];
        fixedIndices(counterInside2) = Xind;
        counterInside2 = counterInside2 + 1;
    elseif sum(isOutsideCircle) == 1
        Xind = triInd(~isOutsideCircle);
        Xinside = X(Xind,:);
        theta1 = atan2(Xinside(1,2),Xinside(1,1));
        theta2 = atan2(Xinside(2,2),Xinside(2,1));
        X(Xind,:) = R*[cos(theta1), sin(theta1);
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

% triplot(tri,X(:,1),X(:,2))
% hold on
% theta = linspace(0,2*pi,1000);
% plot(R*cos(theta),R*sin(theta),'red')
% axis equal
% hold off
% 
% 
% 



