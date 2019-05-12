function [Q, W] = getQuadFromFile(N,type)

if nargin < 2
    type = 'double';
end
switch type
    case 'mp'
        M = mp.read(['integration/legendreData/GL_N' num2str(N) '.csv']);
    otherwise
        M = readmatrix(['integration/legendreData/GL_N' num2str(N) '.csv'],'OutputType',type);
end
Q = zeros(N,1,class(M));
W = zeros(N,1,class(M));
indices = floor(N/2)+1:N;
Q(indices) = M(:,1);
W(indices) = M(:,2);
indices = 1:floor(N/2);
Q(indices) = -Q(N-indices+1);
W(indices) = W(N-indices+1);