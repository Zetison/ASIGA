function XYZ_mirror = mirror3D(XYZ, plane)

dim = length(size(XYZ));

if dim == 3
    n = plane(1:3);
    c = plane(4);
    proj_XYZ = zeros(size(XYZ));

    proj_XYZ(:,:,1) = 2*(dot2(XYZ,n) - c)/norm(n)^2 * n(1);
    proj_XYZ(:,:,2) = 2*(dot2(XYZ,n) - c)/norm(n)^2 * n(2);
    proj_XYZ(:,:,3) = 2*(dot2(XYZ,n) - c)/norm(n)^2 * n(3);

    XYZ_mirror = XYZ - proj_XYZ;
elseif dim == 2
    n = plane(1:3);
    c = plane(4);
    proj_XYZ = zeros(size(XYZ));

    proj_XYZ(1,:) = 2*(dot2(XYZ,n) - c)/norm(n)^2 * n(1);
    proj_XYZ(2,:) = 2*(dot2(XYZ,n) - c)/norm(n)^2 * n(2);
    proj_XYZ(3,:) = 2*(dot2(XYZ,n) - c)/norm(n)^2 * n(3);

    XYZ_mirror = XYZ - proj_XYZ;
end