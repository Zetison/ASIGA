function output = objFun(prm)
output = norm(evaluateNURBS(scatterer,[prm(2), prm(3)]) - [sqrt(prm(1)^2-Upsilon^2)*sin(theta)*cos(phi);
                                                                          sqrt(prm(1)^2-Upsilon^2)*sin(theta)*sin(phi);
                                                                          prm(1)*cos(theta)]);