function dirs = findDirsNSD_2D(cps,dgdx,dgdy,dim,zeta)

if dim == 1
    dg = @(xi,r) dgdx{r}(xi,zeta);
else
    dg = @(xi,r) dgdy{r}(zeta,xi);
end
    

dirs = NaN(2*size(cps,2)-2,1);
for i = [1,size(cps,2)]
    xi = cps(1,i); 
    r = cps(2,i);
    if r == 0
        if i == 1
            dirs(1:2) = sign(imag(1i/sign(dg(xi,1))));
        elseif i == size(cps,2)
            dirs(end-1:end) = sign(imag(1i/sign(dg(xi,1))));
        end
    end
end
for i = 1:size(cps,2)
    xi = cps(1,i); 
    r = cps(2,i);
    if r < 2
        sgnimg = sign(imag(exp(1i*pi/2/(r+1))/sign(dg(xi,r+1))));
        if i == 1
            dirs(1:2) = sgnimg;
        elseif i == size(cps,2)
            if r == 0
                dirs(end-1:end) = sgnimg;
            else % r == 1
                dirs(end-1:end) = -sgnimg;
            end
        else
            dirs(2*i-3:2*i) = sgnimg*[-1,-1,1,1];
        end
    end
end
while any(isnan(dirs))
    for i = 1:size(cps,2)
        r = cps(2,i);
        if r >= 2
            if i == 1
                if ~isnan(dirs(2))
                    dirs(1) = dirs(2);
                end
            elseif i == size(cps,2)
                if ~isnan(dirs(end-1))
                    dirs(end) = dirs(end-1);
                end
            else
                prev = dirs(2*i-3);
                next = dirs(2*i);
                if ~isnan(prev) && ~isnan(next)
                    dirs(2*i-2:2*i-1) = [prev,next];
                else
                    if mod(r,2) % r is odd
                        if ~isnan(next)
                            dirs(2*i-2:2*i-1) = prev*[1,-1];
                        elseif ~isnan(next)
                            dirs(2*i-2:2*i-1) = next*[-1,1];
                        end
                    else
                        if ~isnan(prev)
                            if sign(imag(exp(1i*pi/2/(r+1))/sign(dg(xi,r+1)))) == prev
                                dirs(2*i-2:2*i-1) = prev*[1,1];
                            else
                                dirs(2*i-2:2*i-1) = prev*[1,-1];
                            end
                        elseif ~isnan(next)
                            if sign(imag(exp(1i*pi/2/(r+1))/sign(dg(xi,r+1)))) == next
                                dirs(2*i-2:2*i-1) = next*[1,1];
                            else
                                dirs(2*i-2:2*i-1) = next*[-1,1];
                            end
                        end
                    end     
                end
            end
        end
    end
end
for i = 1:size(cps,2)
    xi = cps(1,i); 
    r = cps(2,i);
    turn = sign(imag(exp(1i*pi/2/(r+1))/sign(dg(xi,r+1)))) == -1;
    if turn
        j = [1:r,0];
    else
        j = 0:r;
    end
    if i == 1
        if dirs(1) == 1
            dirs(1) = j(1);
        else
            dirs(1) = j(end);
        end
    elseif i == size(cps,2)
        if dirs(end) == 1
            dirs(end) = j(1);
        else
            dirs(end) = j(end);
        end
    else
        if turn
            idx = floor((r+1)/2);
        else
            idx = ceil((r+1)/2);
        end
        if dirs(2*i-2) == 1
            dirs(2*i-2) = j(idx);
        else
            dirs(2*i-2) = j(idx+1);
        end
        if dirs(2*i-1) == 1
            dirs(2*i-1) = j(1);
        else
            dirs(2*i-1) = j(end);
        end
    end
end