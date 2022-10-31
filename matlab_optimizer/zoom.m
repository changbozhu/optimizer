function [ alpha_star ] = zoom(alpha_lo,alpha_hi,zoomPara)
%zoom function is compatible with Line Search.
    c1 = zoomPara.c1;
    c2 = zoomPara.c2;
    xk = zoomPara.xk;
    pk = zoomPara.pk;
    interp_method = zoomPara.interp_method;
    obj_func = zoomPara.obj_func;
    obj_grad = zoomPara.obj_grad;
    DPhi_0 = zoomPara.dphi0;
    Phi_0  = phi(obj_func, xk, pk, 0);
    Phi_lo = phi(obj_func, xk, pk, alpha_lo);
    iterINF = 5;
    zoomiter = 0;
    interp_func = str2func(interp_method);
    while zoomiter<iterINF
        %disp(['alpha_lo = ' num2str(alpha_lo) ';' 'alpha_hi = ' num2str(alpha_hi) ])
        %Interpolate (using quadratic,cubic or bisection) 
        %to find a trial steplength alpha_j 
        %between alpha_lo and alpha_hi;
        alpha_j = interp_func(alpha_lo,alpha_hi);    
        
        %Evaluate Phi(alpha_j,phipara) and  Phi(alpha_lo,phipara)
        Phi_j  = phi(obj_func, xk, pk, alpha_j);
        
        
        %Wolfe Condition
        Gline = Phi_0 + c1*alpha_j*DPhi_0;
        if  ( Phi_j > Gline || Phi_j >= phi(obj_func, xk, pk, alpha_lo) )
            alpha_hi = alpha_j ;           
        else
            %Evaluate DPhi(alpha_j)
            DPhi_j = dphi(obj_grad, xk, pk, alpha_j);
            
            if (abs(DPhi_j)<=-c2*DPhi_0)
                alpha_star=alpha_j;
                break;
            end
            ctemp = DPhi_j*(alpha_hi-alpha_lo);
            if ctemp>=0
                alpha_hi=alpha_lo;
            end                
            
            alpha_lo=alpha_j ;         
        end           
        zoomiter=zoomiter+1;
    end
    if zoomiter >= iterINF
        disp('  ')
        disp('Please notice:')
        disp('Warning: zoom terminated because the maximum number of iterations is reached.')
        alpha_star = 1;
    end
    return;
end

