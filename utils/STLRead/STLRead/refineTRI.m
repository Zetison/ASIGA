function [P_new,TRI_new] = refineTRI(P,TRI)
temp = [P; meanTRI(P,TRI,[1,2]); meanTRI(P,TRI,[1,3]); meanTRI(P,TRI,[2,3])];
P_new = uniquetol(temp,'ByRows',true);
TRI_new = delaunay(P_new(:,1),P_new(:,2));

function edgeNode = meanTRI(P,TRI,idx)

edgeNode = zeros(size(TRI,1),3);
for i = 1:3
    edgeNode(:,i) = mean([P(TRI(:,idx(1)),i),P(TRI(:,idx(2)),i)],2);
end