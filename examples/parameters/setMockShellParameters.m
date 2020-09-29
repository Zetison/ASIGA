function layer = setMockShellParameters(noDomains)

layer{1} = struct('media', 'fluid', ...
                  't', 0.008, ...
                  'R', 1, ...
                  'L',  10, ...
                  'c_f', 1500, ...
                  'rho', 1000);
layer{2} = struct('media', 'solid', ...
                  'E', 210e9, ...
                  'nu', 0.3, ...
                  'rho', 7850);
layer{3} = struct('media', 'fluid', ...
                  'c_f', 1500, ...
                  'rho', 1000);

layer = layer(1:noDomains);