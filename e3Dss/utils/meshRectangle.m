function [tri, X] = meshRectangle(xl,xu,N)

Lx = xu(1)-xl(1);
Ly = xu(2)-xl(2);
d = Lx/(N-1);
h = min(sqrt(3)/2*d,Ly/2);
noLayers = floor(Ly/(2*h));
Ny = 2*noLayers+1;
x1 = linspace(xl(1),xu(1),N);
x2 = [xl(1), x1(1:end-1)+d/2, xu(1)];
x3 = [x1,x2];
X = zeros(noLayers*length(x3)+N,2);
X(1:noLayers*length(x3),1) = repmat(x3,1,noLayers);
X(end-N+1:end,1) = x1;

N1 = length(x1);
N2 = length(x2);
N = N1+N2;
hy = Ly/(Ny-1);
for i = 1:noLayers
    X(N*(i-1)+1:N*(i-1)+N1,2) = (xl(2)+2*hy*(i-1))*ones(1,N1);
    X(N*(i-1)+N1+1:N*(i-1)+N1+N2,2) = (xl(2)+(2*i-1)*hy)*ones(1,N2);
end
X(end-N1+1:end,2) = xu(2)*ones(N1,1);
% keyboard
tri = delaunay(X(:,1),X(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% triplot(tri,X(:,1),X(:,2))
% % for i = 1:size(X,1)
% %     text(X(i,1),X(i,2),num2str(i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% % end
% % for i = 1:size(tri,1)
% %     x = mean(X(tri(i,:),:));
% %     text(x(1),x(2),num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','color','red')
% % end
% axis equal
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



