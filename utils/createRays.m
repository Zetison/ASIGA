function varCol = createRays(varCol)


switch varCol.applyLoad
    case 'pointSource'
        error('not implemented')
    case 'planeWave'
        X = varCol.controlPts;
        d_vec = varCol.d_vec;
        [X_m, A] = orthogonalTransform(X, d_vec);

        convexHull = boundary(X_m(:,1),X_m(:,2),0);
        n = round(10^(varCol.N/2));
%         n = round(2^parm)-1;
        minX_m_x = min(X_m(:,1));
        minX_m_y = min(X_m(:,2));
        maxX_m_x = max(X_m(:,1));
        maxX_m_y = max(X_m(:,2));
        delta_x = (maxX_m_x-minX_m_x)/(n-1);
        delta_y = sqrt(3)/2*delta_x;
        maxX_m_x = maxX_m_x + delta_x/2;
        maxX_m_y = maxX_m_y + delta_y/2;
        minX_m_x = minX_m_x - delta_x/2;
        minX_m_y = minX_m_y - delta_y/2;
        X = createBeamTriangulation([minX_m_x,minX_m_y],[maxX_m_x,maxX_m_y],n);
        
        varCol.beamEnergy = abs(varCol.P_inc)^2*3*sqrt(3)/2*delta_x^2;
        XX_m = X(:,1).';
        YY_m = X(:,2).';
        
        indices = (1:numel(XX_m))';
        beams = [indices,indices-n,indices-n+1,indices+1,indices+n,indices+n-1,indices-1];
        IN = inpolygon(XX_m,YY_m,X_m(convexHull,1),X_m(convexHull,2));
        temp = zeros(numel(XX_m),1);
        temp(IN) = 1:sum(IN);
        bndrRays = zeros(1,numel(XX_m),'logical');
        bndrRays([1:(n-1),end-n+2:end,n:2*n-1:end,2*n-1:2*n-1:end,1:2*n-1:end,n-1:2*n-1:end]) = true;
        beams = beams(and(IN, ~bndrRays),:);
        beams = reshape(temp(beams(:)),size(beams,1),size(beams,2));
        beams(any(~beams,2),:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        varCol.oXX_m = XX_m;
        varCol.oYY_m = YY_m;
        varCol.X_m = X_m;
        varCol.convexHull = convexHull;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        XX_m = XX_m(IN).';
        YY_m = YY_m(IN).';
        
        o = (A*[XX_m,YY_m,2*ones(size(XX_m))*min(X_m(:,3))].').';
        varCol.o = o;
        varCol.beams = beams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         close all
%         figure(1)
%         hold on
%         axis equal
%         plot(X_m(convexHull,1),X_m(convexHull,2),'o-','color','red','MarkerFaceColor','red')
%         set(gcf,'color','w');
% %         export_fig('../graphics/S1/createRay1', '-pdf', '-transparent')
%         plot(varCol.oXX_m,varCol.oYY_m,'o','color','green','MarkerFaceColor','green')
% %         ylim([-1,1])
% %         export_fig('../graphics/S1/createRay2', '-pdf', '-transparent')
%         plot(XX_m,YY_m,'o','color','blue','MarkerFaceColor','blue')
% %         export_fig('../graphics/S1/createRay3', '-pdf', '-transparent')
% %         for i = 1:size(XX_m,1)
% %             text(XX_m(i),YY_m(i),num2str(i))
% %         end
%         for i = 1:size(beams,1)
%             indices = beams(i,:);
%             h = plot(XX_m(indices([2,3,4,5,6,7,2,1,3,4,1,5,6,1,7])),YY_m(indices([2,3,4,5,6,7,2,1,3,4,1,5,6,1,7])),'black');
% %             export_fig('../graphics/S1/createRay4', '-pdf', '-transparent') 
%             pause(0.2)
%             set(h,'Visible','off')
%         end
%         keyboard
% 
%         varCol.convexHull = convexHull;
%         varCol.X_m = X_m;
% %         
end


function X = createBeamTriangulation(xl,xu,Nx)
xu2 = xu;
Lx = xu2(1)-xl(1);
Ly = xu2(2)-xl(2);
d = Lx/(Nx-1);
h = sqrt(3)/2*d;
Ny = 2*ceil(Ly/h/2)+1;
xu2(2) = xl(2) + (Ny-1)*h;
x2 = linspace(xl(1),xu2(1),Nx);
x1 = x2(1:end-1)+d/2;
x3 = [x1,x2];
noBeamRows = floor(Ny/2);
X = zeros(Nx*Ny-floor(Ny/2)-1,2);
X(1:noBeamRows*length(x3),1) = repmat(x3,1,noBeamRows);

X(end-Nx+2:end,1) = x1;

N1 = length(x2);
N2 = length(x1);
N = N1+N2;
hy = (xu2(2)-xl(2))/(Ny-1);
for i = 1:noBeamRows
    X(N*(i-1)+1:N*(i-1)+N2,2) = (xl(2)+2*hy*(i-1))*ones(1,N2);
    X(N*(i-1)+N2+1:N*i,2) = (xl(2)+(2*i-1)*hy)*ones(1,N1);
end
X(end-N2+1:end,2) = xu2(2)*ones(N2,1);
% [~,I] = min(norm2(X));
% X(:,1) = X(:,1) - X(I,1);
% X(:,2) = X(:,2) - X(I,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tri = delaunay(X(:,1),X(:,2));
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
% % for i = 1:size(X,1)
% %     text(X(i,1),X(i,2),num2str(i))
% % end
% 
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















