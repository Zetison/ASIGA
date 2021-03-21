function task = computeDerivedFFPquantities(task,p_h)
v = getFarFieldPoints(task);

task.results.p = p_h;
task.results.p_Re = real(p_h);
task.results.p_Im = imag(p_h);
task.results.abs_p = abs(p_h);
task.results.TS = 20*log10(abs(p_h/task.misc.P_inc));
if task.analyticSolutionExist
    if task.ffp.plotFarField
        p_ref = task.p_0_(v);
%         p_ref = exactKDT(task.misc.omega/task.varCol{1}.c_f,task.misc.P_inc,task.varCol{1}.R_i);
    else
        p_ref = task.varCol{1}.p_(v);
    end
    task.results.p_ref = p_ref;
    task.results.p_Re_ref = real(p_ref);
    task.results.p_Im_ref = imag(p_ref);
    task.results.abs_p_ref = abs(p_ref);
    task.results.TS_ref = 20*log10(abs(p_ref/task.misc.P_inc));

    task.results.error_pAbs = 100*abs(abs(p_ref)-abs(p_h))./abs(p_ref);
    task.results.error_p = 100*abs(p_ref-p_h)./abs(p_ref);
end