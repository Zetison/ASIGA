function layer = setM4Parameters(noDomains)
if nargin < 1
    noDomains = 1;
end

layer{1} = struct('media', 'fluid', ...
                  'R', 3, ...
                  't', 0.02, ...
                  'c_f', 1500, ...
                  'rho', 1000);
layer{2} = struct('media', 'solid', ...
                  'E', 210e9, ...
                  'nu', 0.3, ...
                  'rho', 7850);

layer = layer(1:noDomains);