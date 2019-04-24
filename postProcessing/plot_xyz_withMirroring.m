
plot3(xyz(1,:), xyz(2,:), xyz(3,:), 'color',lineColor)

hold on
if size(mirrorMatrix,1) > 0
    xyz_mirror = mirror3D(xyz, mirrorMatrix(1,:));

    plot3(xyz_mirror(1,:), xyz_mirror(2,:), xyz_mirror(3,:), 'color',lineColor)

    if size(mirrorMatrix,1) > 1

        xyz_mirror2 = mirror3D(xyz, mirrorMatrix(2,:));

        plot3(xyz_mirror2(1,:), xyz_mirror2(2,:), xyz_mirror2(3,:), 'color',lineColor)

        xyz_mirror3 = mirror3D(xyz_mirror, mirrorMatrix(2,:));

        plot3(xyz_mirror3(1,:), xyz_mirror3(2,:), xyz_mirror3(3,:), 'color',lineColor)

        if size(mirrorMatrix,1) > 2

            xyz_mirror4 = mirror3D(xyz, mirrorMatrix(3,:));

            plot3(xyz_mirror4(1,:), xyz_mirror4(2,:), xyz_mirror4(3,:), 'color',lineColor)

            xyz_mirror5 = mirror3D(xyz_mirror, mirrorMatrix(3,:));

            plot3(xyz_mirror5(1,:), xyz_mirror5(2,:), xyz_mirror5(3,:), 'color',lineColor)

            xyz_mirror6 = mirror3D(xyz_mirror2, mirrorMatrix(3,:));

            plot3(xyz_mirror6(1,:), xyz_mirror6(2,:), xyz_mirror6(3,:), 'color',lineColor)

            xyz_mirror7 = mirror3D(xyz_mirror3, mirrorMatrix(3,:));

            plot3(xyz_mirror7(1,:), xyz_mirror7(2,:), xyz_mirror7(3,:), 'color',lineColor)
        end
    end
end