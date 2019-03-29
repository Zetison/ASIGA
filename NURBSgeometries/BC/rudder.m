function [X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta] = rudder(l_l,l_u,h,delta,R_x,T,sgn,f_l,dfxi_l,dfdxi2_l,f_u,dfxi_u,dfdxi2_u)


X = @(t) repmat(T,1,size(t,2)) + R_x*[-(l_l*t(1,:).^2+t(2,:).*(delta+(l_u-l_l)*t(1,:).^2));
                                      sgn*(l_l*f_l(t(1,:)) + t(2,:).*(l_u*f_u(t(1,:))-l_l*f_l(t(1,:))));
                                      t(2,:)*h];
dXdxi = @(t) R_x*[-(2*l_l*t(1,:)+t(2,:).*(2*(l_u-l_l)*t(1,:)));
                  sgn*(l_l*dfxi_l(t(1,:)) + t(2,:).*(l_u*dfxi_u(t(1,:))-l_l*dfxi_l(t(1,:))));
                  zeros(1,size(t,2))];
dXdeta = @(t) R_x*[-(delta+(l_u-l_l)*t(1,:).^2);
                  l_u*f_u(t(1,:))-l_l*f_l(t(1,:));
                  ones(1,size(t,2))*h];
d2Xdxi2 = @(t) R_x*[-(2*l_l+t(2,:).*(2*(l_u-l_l)));
                    sgn*(l_l*dfdxi2_l(t(1,:)) + t(2,:).*(l_u*dfdxi2_u(t(1,:))-l_l*dfdxi2_l(t(1,:))));
                      zeros(1,size(t,2))];
d2Xdeta2 = @(t) zeros(3,size(t,2));
d2Xdxideta = @(t) R_x*[-(delta+2*(l_u-l_l)*t(1,:));
                  sgn*(l_u*dfxi_u(t(1,:))-l_l*dfxi_l(t(1,:)));
                  zeros(1,size(t,2))];
