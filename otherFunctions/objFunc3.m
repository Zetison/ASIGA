function f = objFunc3(options,alpha_f_arr,beta_f)

    v = options.R_o(1)*[cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]';
    data = e3Dss(v, options);
    
    f = 20*log10(abs(data(1).p)).';
end
