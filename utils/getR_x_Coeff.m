function R_xScaled = getR_x_Coeff(R_x,useEnrichedBfuns,k,d_vec,x,useRegul,integrals,sgn,constants,...
                psiType,useCBIE,useHBIE,dXIdv,dR_xdxi,dR_xdeta,v_1,v_2,alpha)
if useEnrichedBfuns
    R_x = R_x*exp(1i*k*dot(d_vec, x));
end
if useRegul
    temp2 = integrals{3} - integrals{1};
    switch psiType
        case 1
            C1 = constants{2};
            temp2 = temp2 + 0.5*(1-sgn - (1+1i/(k*C1))*(1-exp(2*1i*k*C1)));
        case 2
            temp2 = temp2 + 0.5*(1-sgn);
        case 3
            temp2 = temp2 - 0.5*(1+sgn);
    end
    R_xScaled = R_x*temp2;
else
    R_xScaled = complex(zeros(size(R_x)));
    if useCBIE
        R_xScaled = R_x*(-0.5*(1+sgn) - integrals{1});
    end
    if useHBIE
        dphidv = dXIdv*[dR_xdxi; dR_xdeta];
        gR_x = dphidv(1,:).'*v_1 + dphidv(2,:).'*v_2;
        if useEnrichedBfuns
            gR_x = exp(1i*k*dot(d_vec, x))*gR_x + 1i*k*(d_vec*R_x).';
        end
        R_xScaled = R_xScaled + alpha*(-R_x*integrals{2} + (gR_x*integrals{3}).');
    end
end