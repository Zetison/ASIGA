function quivers2(x,y,u,v,xLims,yLims,color)

dx = xLims(2)-xLims(1);
dy = yLims(2)-yLims(1);
            
theta = atan2(v,u);
s = 1/sqrt(cos(theta)^2/dx^2+sin(theta)^2/dy^2);
s(theta == pi/2) = dy;
X = [x; x+0.1*u*s];
Y = [y;y+0.1*v*s];
plot(X, Y,'color',color,'LineWidth',2);

phi = 20*pi/180;
Xs = X/dx;
Ys = Y/dy;
scale = 0.3;
rs = norm([Xs(2)-Xs(1),Ys(2)-Ys(1)]);
thetas = atan2(Ys(2)-Ys(1),Xs(2)-Xs(1));
X3s = Xs(2) + rs*scale*cos(pi+thetas+phi);
Y3s = Ys(2) + rs*scale*sin(pi+thetas+phi);
X3 = [X(2); X3s*dx];
Y3 = [Y(2); Y3s*dy];
plot(X3, Y3,'color',color);
X4s = Xs(2) + rs*scale*cos(pi+thetas-phi);
Y4s = Ys(2) + rs*scale*sin(pi+thetas-phi);
X4 = [X(2); X4s*dx];
Y4 = [Y(2); Y4s*dy];
plot(X4, Y4,'color',color);
fill([X(2),X3(2),X4(2)],[Y(2),Y3(2),Y4(2)],color,'EdgeColor',color)
