function invA = inv3(A,J)
% Formula on http://mathworld.wolfram.com/MatrixInverse.html
invA = zeros(size(A));
invA(1,1,:,:,:) = det2(A(2:3,2:3,:,:,:))./J;
invA(2,1,:,:,:) = det2(A(2:3,[3,1],:,:,:))./J;
invA(3,1,:,:,:) = det2(A(2:3,1:2,:,:,:))./J;
invA(1,2,:,:,:) = det2(A([1,3],[3,2],:,:,:))./J;
invA(2,2,:,:,:) = det2(A([1,3],[1,3],:,:,:))./J;
invA(3,2,:,:,:) = det2(A([1,3],[2,1],:,:,:))./J;
invA(1,3,:,:,:) = det2(A(1:2,2:3,:,:,:))./J;
invA(2,3,:,:,:) = det2(A(1:2,[3,1],:,:,:))./J;
invA(3,3,:,:,:) = det2(A(1:2,1:2,:,:,:))./J;