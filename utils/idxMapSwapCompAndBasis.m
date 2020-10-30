function I = idxMapSwapCompAndBasis(d_f,n_en)
% Create index map with the first index looping over the d_f components
% instead of the n_en basis functions

I = zeros(d_f*n_en,d_f*n_en);
for i = 1:d_f
    for j = 1:d_f
        I(i:d_f:end, j:d_f:end) = (1+(i-1)*n_en:i*n_en).'  + ((1+(j-1)*n_en:j*n_en)-1)*d_f*n_en;
    end
end
I = reshape(I,(d_f*n_en)^2,1);