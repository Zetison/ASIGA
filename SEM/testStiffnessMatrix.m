
% 
% Kt = zeros(noDofsPatch,noDofsPatch);
% for i = 1:Nxi
%     for j = 1:Neta
%         for l = 1:Nzeta
%             for it = 1:Nxi
%                 for jt = 1:Neta
%                     for lt = 1:Nzeta
%                         tempXi = 0;
%                         for alpha = 1:Nxi
%                             tempXi = tempXi + G(1,1,alpha,j,l)*Dxi(i,alpha)*Dxi(it,alpha);
%                         end
%                         tempEta = 0;
%                         for beta = 1:Neta
%                             tempEta = tempEta + G(2,2,i,beta,l)*Deta(j,beta)*Deta(jt,beta);
%                         end
%                         tempZeta = 0;
%                         for gamma = 1:Nzeta
%                             tempZeta = tempZeta + G(3,3,i,j,gamma)*Deta(l,gamma)*Deta(lt,gamma);
%                         end
%                         Kt(it+(jt-1)*Nxi+(lt-1)*Nxi*Neta,i+(j-1)*Nxi+(l-1)*Nxi*Neta) = ...
%                               delta(lt,l)*G(1,2,i,jt,l)*Dxi(it,i)*Deta(j,jt) ...
%                             + delta(lt,l)*G(1,2,it,j,l)*Dxi(i,it)*Deta(jt,j) ...
%                             + delta(jt,j)*G(1,3,i,j,lt)*Dxi(it,i)*Dzeta(l,lt) ...
%                             + delta(jt,j)*G(1,3,it,j,l)*Dxi(i,it)*Dzeta(lt,l) ...
%                             + delta(it,i)*G(2,3,i,j,lt)*Deta(jt,j)*Dzeta(l,lt) ...
%                             + delta(it,i)*G(2,3,i,jt,l)*Deta(j,jt)*Dzeta(lt,l) ...
%                             + delta(jt,j)*delta(lt,l)*tempXi ...
%                             + delta(it,i)*delta(lt,l)*tempEta ...
%                             + delta(it,i)*delta(jt,j)*tempZeta;
%                     end
%                 end
%             end
%         end
%     end
% end

Kt = zeros(noDofsPatch,noDofsPatch);
for i = 1:Nxi
    for j = 1:Neta
        for l = 1:Nzeta
            for it = 1:Nxi
                for jt = 1:Neta
                    for lt = 1:Nzeta
                        w = 0;
                        for alpha = 1:Nxi
                            for beta = 1:Neta
                                for gamma = 1:Nzeta
                                    w = w + [  Dxi(it,alpha)*delta(jt,beta)*delta(lt,gamma), ...
                                             delta(it,alpha) *Deta(jt,beta)*delta(lt,gamma), ...
                                             delta(it,alpha)*delta(jt,beta)*Dzeta(lt,gamma)]*G(:,:,alpha,beta,gamma) ...
                                           *[  Dxi(i,alpha)*delta(j,beta)*delta(l,gamma);
                                             delta(i,alpha)* Deta(j,beta)*delta(l,gamma);
                                             delta(i,alpha)*delta(j,beta)*Dzeta(l,gamma)];
                                end
                            end
                        end
                        Kt(it+(jt-1)*Nxi+(lt-1)*Nxi*Neta,i+(j-1)*Nxi+(l-1)*Nxi*Neta) = w;
                    end
                end
            end
        end
    end
end

K = reshape(Kloc,noDofsPatch,noDofsPatch);
sparsityIdx = nnz(Kt)/noDofsPatch^2
ttt = K-Kt;
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