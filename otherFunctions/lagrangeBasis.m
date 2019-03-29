function [p,dpdx] = lagrangeBasis(x,i,Xi,j)
Xi = Xi.'; % avoid matlab 'bug' for p = 1
p = prod((x-Xi(j))./(Xi(i)-Xi(j)),2).';
Xi = Xi.';
N = size(i,1);
i_arr = i(:,1).';
dpdx = zeros(1,N);
for i = i_arr
    for j = 1:N
        if j ~= i
            if i < j
                l = [1:i-1,i+1:j-1,j+1:N];
            else
                l = [1:j-1,j+1:i-1,i+1:N];
            end
            dpdx(i) = dpdx(i) + prod((x-Xi(l))./(Xi(i)-Xi(l)),2)/(Xi(i)-Xi(j));
        end
    end
end