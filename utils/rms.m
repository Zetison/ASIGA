function A = rms(B)
n = numel(B);
A = sqrt(1/n*sum(B(:).^2));