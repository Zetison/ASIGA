function nurbs = GordonHall(faces)
% Algorithm based on that by William Gordan and Charles Hall available in
% their article: https://doi.org/10.1002/nme.1620070405
% Another openly available reference: https://wiki.math.ntnu.no/_media/ma8502/2014h/deformed2.pdf
% It is assumed that the faces has the ordering xi=0, xi=1, eta=0, eta=1,
% and so on.

d_p = faces{1}.d_p+1;
d = faces{1}.d;
switch d_p
    case 2
        degree = [faces{3}.degree,faces{1}.degree];
        knots  = [faces{3}.knots, faces{1}.knots];
        number = [faces{3}.number, faces{1}.number];
        g = cell(1,d_p);
        for i = 1:numel(knots)
            g{i} = aveknt(knots{i}, degree(i)+1);
        end

        g{1} = reshape(g{1},1,[]);
        g{2} = reshape(g{2},1,1,[]);
        P = zeros([d+1,number]);
        P(:,1,:)    = reshape(faces{1}.coeffs,d+1,1,[]);
        P(:,end,:)  = reshape(faces{2}.coeffs,d+1,1,[]);
        P(:,:,1)    = faces{3}.coeffs;
        P(:,:,end)  = faces{4}.coeffs;
        P =  P(:,1,:).*(1-g{1}) + P(:,end,:).*g{1} + P(:,:,1).*(1-g{2}) + P(:,:,end).*g{2} ...
           - P(:,1,1).*(1-g{1}).*(1-g{2}) - P(:,end,1).*g{1}.*(1-g{2}) - P(:,1,end).*(1-g{1}).*g{2} - P(:,end,end).*g{1}.*g{2};
        nurbs = createNURBSobject(P,knots);   
    case 3
        degree = [faces{3}.degree(2),faces{1}.degree];
        knots = [faces{3}.knots{2},faces{1}.knots];
        number = [faces{3}.number(2),faces{1}.number];
        g = cell(1,d_p);
        for i = 1:numel(knots)
            g{i} = aveknt(knots{i}, degree(i)+1);
        end

        g{1} = reshape(g{1},1,[]);
        g{2} = reshape(g{2},1,1,[]);
        g{3} = reshape(g{3},1,1,1,[]);
        P = zeros([d+1,number]);
        P(:,1,:,:)   = faces{1}.coeffs(1:d+1,:,:);
        P(:,end,:,:) = faces{2}.coeffs(1:d+1,:,:);
        P(:,:,1,:)   = permute(faces{3}.coeffs(1:d+1,:,:),[1,3,2]);
        P(:,:,end,:) = permute(faces{4}.coeffs(1:d+1,:,:),[1,3,2]);
        P(:,:,:,1)   = faces{5}.coeffs(1:d+1,:,:);
        P(:,:,:,end) = faces{6}.coeffs(1:d+1,:,:);
        P =    P(:,1,:,:).*(1-g{1}) + P(:,end,:,:).*g{1} + P(:,:,1,:).*(1-g{2}) + P(:,:,end,:).*g{2} + P(:,:,:,1).*(1-g{3}) + P(:,:,:,end).*g{3} ...
             - P(:,1,1,:).*(1-g{1}).*(1-g{2}) - P(:,1,end,:).*(1-g{1}).*g{2} - P(:,1,:,1).*(1-g{1}).*(1-g{3}) ...
             - P(:,1,:,end).*(1-g{1}).*g{3} - P(:,end,1,:).*g{1}.*(1-g{2}) - P(:,end,end,:).*g{1}.*g{2} - P(:,end,:,1).*g{1}.*(1-g{3}) ...
             - P(:,end,:,end).*g{1}.*g{3} - P(:,:,1,1).*(1-g{2}).*(1-g{3}) - P(:,:,1,end).*(1-g{2}).*g{3} ...
             + P(:,1,1,1).*(1-g{1}).*(1-g{2}).*(1-g{3}) + P(:,1,1,end).*(1-g{1}).*(1-g{2}).*g{3} ...
             + P(:,end,1,1).*g{1}.*(1-g{2}).*(1-g{3}) + P(:,end,1,end).*g{1}.*(1-g{2}).*g{3} - P(:,:,end,1).*g{2}.*(1-g{3}) ...
             - P(:,:,end,end).*g{2}.*g{3} + P(:,1,end,1).*(1-g{1}).*g{2}.*(1-g{3}) + P(:,1,end,end).*(1-g{1}).*g{2}.*g{3} ...
             + P(:,end,end,1).*g{1}.*g{2}.*(1-g{3}) + P(:,end,end,end).*g{1}.*g{2}.*g{3};
        nurbs = createNURBSobject(P,knots);
    otherwise
        error('not implemented')
end