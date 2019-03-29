function p = farField_(v,options)

options.calc_farField = true;
d_vec = options.d_vec;
if size(d_vec,2) > 1
    nPts = size(v,1);
    v = v(1,:); % Assume monostatic scattering
    options.d_vec = d_vec(:,1);
    data = e3Dss(v,options);
    p = data(1).p;
    p = p*ones(nPts,1);
else
    data = e3Dss(v,options);
    p = data(1).p;
end