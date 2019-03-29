
Kinft = zeros(Nxi,Neta,N,Nxi,Neta,N);
for i = 1:Nxi
    for j = 1:Neta
        for l = 1:N
            for it = 1:Nxi
                for jt = 1:Neta
                    for lt = 1:N
                        X = c(1:3,it,jt,end);
                        [~, theta] = evaluateProlateCoords2(X(1),X(2),X(3),Upsilon);
                        
                        % Add A1*B1
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(i,it)*delta(j,jt)*rho{1}(it)*rho{2}(jt)*sin(theta)*J3(1,1,it,jt)*BB(lt,l,1);
                        
                        % Add A2*B2
                        temp = 0;
                        for alpha = 1:Nxi
                            temp = temp + J5(1,1,alpha,j)*Dxi(i,alpha)*Dxi(it,alpha);
                        end
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(j,jt)*temp*BB(lt,l,2);
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + J5(1,2,it,j)*Dxi(i,it)*Deta(jt,j)*BB(lt,l,2);
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + J5(1,2,i,jt)*Dxi(it,i)*Deta(j,jt)*BB(lt,l,2);
                        temp = 0;
                        for beta = 1:Neta
                            temp = temp + J5(1,3,i,beta)*Deta(j,beta)*Deta(jt,beta);
                        end
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(i,it)*temp*BB(lt,l,2);
                        
                        % Add A3*B3
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(i,it)*delta(j,jt)*rho{1}(it)*rho{2}(jt)*cos(theta)^2*sin(theta)*J3(1,1,it,jt)*BB(lt,l,3);
                        
                        % Add A4*B4
                        temp = 0;
                        for alpha = 1:Nxi
                            temp = temp + J5(2,1,alpha,j)*Dxi(i,alpha)*Dxi(it,alpha);
                        end
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(j,jt)*temp*BB(lt,l,4);
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + J5(2,2,it,j)*Dxi(i,it)*Deta(jt,j)*BB(lt,l,4);
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + J5(2,2,i,jt)*Dxi(it,i)*Deta(j,jt)*BB(lt,l,4);
                        temp = 0;
                        for beta = 1:Neta
                            temp = temp + J5(2,3,i,beta)*Deta(j,beta)*Deta(jt,beta);
                        end
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(i,it)*temp*BB(lt,l,4);
                        
                        % Add A5*B5
                        temp = 0;
                        for alpha = 1:Nxi
                            temp = temp + J5(3,1,alpha,j)*Dxi(i,alpha)*Dxi(it,alpha);
                        end
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(j,jt)*temp*BB(lt,l,5);
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + J5(3,2,it,j)*Dxi(i,it)*Deta(jt,j)*BB(lt,l,5);
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + J5(3,2,i,jt)*Dxi(it,i)*Deta(j,jt)*BB(lt,l,5);
                        temp = 0;
                        for beta = 1:Neta
                            temp = temp + J5(3,3,i,beta)*Deta(j,beta)*Deta(jt,beta);
                        end
                        Kinft(it,jt,lt,i,j,l) = Kinft(it,jt,lt,i,j,l) + delta(i,it)*temp*BB(lt,l,5);
                    end
                end
            end
        end
    end
end
Kinft = reshape(Kinft,dofsInInfElements,dofsInInfElements);
Kinf = reshape(Kinf,dofsInInfElements,dofsInInfElements);
sparsityIdx = nnz(Kinft)/dofsInInfElements^2
ttt = Kinf-Kinft;
rms(rms(ttt))
% spy(ttt)
tttt = 1;
% fid = fopen('sparsity.txt', 'a+');
% fprintf(fid, '%f\n', nnz(K)/noDofsPatch^2);
% fclose(fid);

function d = delta(i,j)

if i == j
    d = 1;
else
    d = 0;
end
end
