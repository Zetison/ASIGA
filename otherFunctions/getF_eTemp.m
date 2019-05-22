function F_eTemp = getF_eTemp(F_eTemp,useNeumanProj,SHBC,psiType,useCBIE,useHBIE,useRegul,R_x,sctr_x,x,nx,...
                U,dU,p_inc,dp_inc,dpdn,alpha,integrals,k,constants,sgn)

if useNeumanProj
    if SHBC
        if useCBIE
            p_inc_x = R_x*U(sctr_x,:);
        end
        if useHBIE
            dp_inc_x = R_x*dU(sctr_x,:);
        end
    else
        if useHBIE || useRegul
            dpdn_x = R_x*U(sctr_x,:);
        end
    end
else
    if SHBC
        if useCBIE
            p_inc_x = p_inc(x);
        end
        if useHBIE
            dp_inc_x = dp_inc(x,nx);
        end
    else
        if useHBIE || useRegul
            dpdn_x = dpdn(x,nx);
        end
    end
end
if SHBC
    if useCBIE
        F_eTemp = F_eTemp - p_inc_x;
    end
    if useHBIE
        F_eTemp = F_eTemp - alpha*dp_inc_x;
    end
else 
    if useRegul
        temp2 = dpdn_x*(integrals{2} - integrals{4});
        if psiType == 1
            C1 = constants{2};
            C2 = constants{3};
            temp2 = temp2 + dpdn_x*1i*C1/(k*C2)*(1-exp(2*1i*k*C1))/2;
        end
        F_eTemp = F_eTemp + temp2;
    end
    if useHBIE
        F_eTemp = F_eTemp + alpha*dpdn_x*(integrals{1} + 0.5*(1+sgn) - nx*integrals{3});
    end
end