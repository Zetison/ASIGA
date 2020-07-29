function f = objFunc(k,options,v,c_f,total)

    options.omega = k*c_f;
    data = e3Dss(v, options);
    if total
        d_vec = -[0,0,1].';
        p_inc = @(v) 1*exp(1i*dot3(v,d_vec)*k.');
        f = -20*log10(abs(data(1).p + p_inc(v).'));
    else
        f = -20*log10(abs(data(1).p));
%         f = -abs(data(1).p);
    end
end
