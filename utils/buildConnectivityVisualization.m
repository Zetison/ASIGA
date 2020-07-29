function conn = buildConnectivityVisualization(extraKnots, noElems)

conn = zeros(noElems,extraKnots+2);

counter = 1;
counter2 = 1;
for e = 1:noElems
    conn(e,:) = counter2:(counter2+extraKnots+1);
    counter2 = counter2 + extraKnots +1;
    counter = counter + 1;
end