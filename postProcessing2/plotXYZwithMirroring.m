
if plotElementEdges
    surf(XYZ(:,:,1), XYZ(:,:,2), XYZ(:,:,3),W,'EdgeColor','none','LineStyle','none')
    hold on
    if size(mirrorMatrix,1) > 0
        XYZ_mirror = mirror3D(XYZ, mirrorMatrix(1,:));

        surf(XYZ_mirror(:,:,1), XYZ_mirror(:,:,2), XYZ_mirror(:,:,3),W,'EdgeColor','none','LineStyle','none')

        if size(mirrorMatrix,1) > 1

            XYZ_mirror2 = mirror3D(XYZ, mirrorMatrix(2,:));

            surf(XYZ_mirror2(:,:,1), XYZ_mirror2(:,:,2), XYZ_mirror2(:,:,3),W,'EdgeColor','none','LineStyle','none')

            XYZ_mirror3 = mirror3D(XYZ_mirror, mirrorMatrix(2,:));

            surf(XYZ_mirror3(:,:,1), XYZ_mirror3(:,:,2), XYZ_mirror3(:,:,3),W,'EdgeColor','none','LineStyle','none')

            if size(mirrorMatrix,1) > 2

                XYZ_mirror4 = mirror3D(XYZ, mirrorMatrix(3,:));

                surf(XYZ_mirror4(:,:,1), XYZ_mirror4(:,:,2), XYZ_mirror4(:,:,3),W,'EdgeColor','none','LineStyle','none')

                XYZ_mirror5 = mirror3D(XYZ_mirror, mirrorMatrix(3,:));

                surf(XYZ_mirror5(:,:,1), XYZ_mirror5(:,:,2), XYZ_mirror5(:,:,3),W,'EdgeColor','none','LineStyle','none')

                XYZ_mirror6 = mirror3D(XYZ_mirror2, mirrorMatrix(3,:));

                surf(XYZ_mirror6(:,:,1), XYZ_mirror6(:,:,2), XYZ_mirror6(:,:,3),W,'EdgeColor','none','LineStyle','none')

                XYZ_mirror7 = mirror3D(XYZ_mirror3, mirrorMatrix(3,:));

                surf(XYZ_mirror7(:,:,1), XYZ_mirror7(:,:,2), XYZ_mirror7(:,:,3),W,'EdgeColor','none','LineStyle','none')
            end
        end
    end
else
    surf(XYZ(:,:,1), XYZ(:,:,2), XYZ(:,:,3),W)
    hold on
    if size(mirrorMatrix,1) > 0
        XYZ_mirror = mirror3D(XYZ, mirrorMatrix(1,:));

        surf(XYZ_mirror(:,:,1), XYZ_mirror(:,:,2), XYZ_mirror(:,:,3),W)

        if size(mirrorMatrix,1) > 1

            XYZ_mirror2 = mirror3D(XYZ, mirrorMatrix(2,:));

            surf(XYZ_mirror2(:,:,1), XYZ_mirror2(:,:,2), XYZ_mirror2(:,:,3),W)

            XYZ_mirror3 = mirror3D(XYZ_mirror, mirrorMatrix(2,:));

            surf(XYZ_mirror3(:,:,1), XYZ_mirror3(:,:,2), XYZ_mirror3(:,:,3),W)

            if size(mirrorMatrix,1) > 2

                XYZ_mirror4 = mirror3D(XYZ, mirrorMatrix(3,:));

                surf(XYZ_mirror4(:,:,1), XYZ_mirror4(:,:,2), XYZ_mirror4(:,:,3),W)

                XYZ_mirror5 = mirror3D(XYZ_mirror, mirrorMatrix(3,:));

                surf(XYZ_mirror5(:,:,1), XYZ_mirror5(:,:,2), XYZ_mirror5(:,:,3),W)

                XYZ_mirror6 = mirror3D(XYZ_mirror2, mirrorMatrix(3,:));

                surf(XYZ_mirror6(:,:,1), XYZ_mirror6(:,:,2), XYZ_mirror6(:,:,3),W)

                XYZ_mirror7 = mirror3D(XYZ_mirror3, mirrorMatrix(3,:));

                surf(XYZ_mirror7(:,:,1), XYZ_mirror7(:,:,2), XYZ_mirror7(:,:,3),W)
            end
        end
    end
end