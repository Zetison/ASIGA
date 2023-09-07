function plotControlPts(nurbsPatches,colorXi,colorEta,colorZeta,markerColor)
error('Deprecated. Use plotNURBS() instead')
if nargin < 2
    colorXi = 'red';
    colorEta = 'red';
    colorZeta = 'red';
end
if nargin < 5
    markerColor = 'black';
end
markerEdgeColor = markerColor;
if ~iscell(nurbsPatches)
    nurbsPatches = {nurbsPatches};
end
for patch = 1:numel(nurbsPatches)
    nurbs = nurbsPatches{patch};
    switch nurbs.type
        case '3Dvolume'
            for j = 1:size(nurbs.coeffs,3)
                for l = 1:size(nurbs.coeffs,4)
                    v = nurbs.coeffs(1:3,:,j,l);
                    v = reshape(v,3,size(nurbs.coeffs,2));
                    plot3(v(1,:),v(2,:),v(3,:),'o-','color',colorXi,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
                end
            end
            for i = 1:size(nurbs.coeffs,2)
                for j = 1:size(nurbs.coeffs,3)
                    v = nurbs.coeffs(1:3,i,j,:);
                    v = reshape(v,3,size(nurbs.coeffs,4));
                    plot3(v(1,:),v(2,:),v(3,:),'o-','color',colorEta,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
                end
            end
            for i = 1:size(nurbs.coeffs,2)
                for l = 1:size(nurbs.coeffs,4)
                    v = nurbs.coeffs(1:3,i,:,l);
                    v = reshape(v,3,size(nurbs.coeffs,3));
                    plot3(v(1,:),v(2,:),v(3,:),'o-','color',colorZeta,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
                end
            end
        case '3Dsurface'
            for j = 1:size(nurbs.coeffs,3)
                v = nurbs.coeffs(1:3,:,j);
                v = reshape(v,3,size(nurbs.coeffs,2));
                plot3(v(1,:),v(2,:),v(3,:),'o-','color',colorXi,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
            end
            for i = 1:size(nurbs.coeffs,2)
                v = nurbs.coeffs(1:3,i,:);
                v = reshape(v,3,size(nurbs.coeffs,3));
                plot3(v(1,:),v(2,:),v(3,:),'o-','color',colorEta,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
            end

        case '2Dsurface'
            for j = 1:size(nurbs.coeffs,3)
                v = nurbs.coeffs(1:2,:,j);
                v = reshape(v,2,size(nurbs.coeffs,2));
                plot(v(1,:),v(2,:),'o-','color',colorXi,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
            end
            for i = 1:size(nurbs.coeffs,2)
                v = nurbs.coeffs(1:2,i,:);
                v = reshape(v,2,size(nurbs.coeffs,3));
                plot(v(1,:),v(2,:),'o-','color',colorEta,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
            end
        case '1Dnurbs'
            grev = aveknt(nurbs.knots, nurbs.degree+1);
            plot(grev,nurbs.coeffs(1,:),'o-','color',colorXi,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
        case '2Dcurve'
            plot(nurbs.coeffs(1,:),nurbs.coeffs(2,:),'o-','color',colorXi,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
        case '3Dcurve'
            plot3(nurbs.coeffs(1,:),nurbs.coeffs(2,:),nurbs.coeffs(3,:),'o-','color',colorXi,'MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)

    end
end