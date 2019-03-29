function [tri, X] = mesh2DDisk(R,N)


dt = 2*pi/N;
h = sqrt(3)/2*dt*R;
Nr = ceil(R/h)+1;
dr0 = R/(Nr-1);

X = zeros(Nr*N,2);
counter = 1;
for i = 1:Nr-1
    r = R - (i-1)*dr0;
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
% r = 0.0001*R;
% X = [X; r*[cos(theta), sin(theta)]];
X(end+1,:) = [0,0];
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
% hold on
% theta = linspace(0,2*pi,1000);
% plot(R*cos(theta),R*sin(theta),'red')
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






