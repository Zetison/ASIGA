function X = triangulateRectangle(xl,xu,N)

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