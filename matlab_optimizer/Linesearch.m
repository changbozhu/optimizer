function [ alpha_star ] = Linesearch( obj_func, xk, pk, obj_grad, varargin)
    % Line Search Algorithm With the Wolfe Conditions
    %
    if nargin > 4
        options = varargin{1};
        options = Linesearch_options(options);
    else
        options = Linesearch_options();
    end
    
    % Initialize some corresponding parameters
    c1  = options.c1;
    c2  = options.c2;
    % Initialize step length alpha0  and alpha_max
    % Define the left point of step length interval
    alpha0          = options.alpha0;
    alpha_i         = options.alpha_i;  % Initial step length
    alpha_max       = options.alpha_max;
    interp_method   = options.interp_method;
    % get handle of specific interpolation function
    interp_func = str2func(interp_method);
     
    %Evaluate phi(alpha0) and dphi0
    phi0=phi(obj_func, xk, pk, alpha0);
    alpha_i_old = alpha0;
    phi_i_old = phi0;
    dphi0=dphi(obj_grad, xk, pk, alpha0);
    
    zoomPara = struct;
    zoomPara.dphi0 = dphi0;
    zoomPara.c1 = options.c1;
    zoomPara.c2 = options.c2;
    zoomPara.interp_method = options.interp_method;
    zoomPara.xk = xk;
    zoomPara.pk = pk;
    
    zoomPara.obj_func  = obj_func;
    zoomPara.obj_grad = obj_grad;
    
    Lineiter=1;
    while(1)
        %Evaluate phi(alpha_i)
        phi1=phi(obj_func, xk, pk, alpha_i);
        if (  (phi1>phi0+c1*alpha_i*dphi0) || (phi1>=phi_i_old && Lineiter>1) )
            alpha_star=zoom(alpha_i_old,alpha_i,zoomPara);
            break;
        end
        
        %Evaluate  dphi(alpha_i)
        dphi1= dphi(obj_grad, xk, pk, alpha_i);
        if (abs(dphi1)<=-c2*dphi0)
            alpha_star=alpha_i;
            break;
        end
        if (  dphi1>=0   )
            alpha_star=zoom(alpha_i,alpha_i_old,zoomPara) ;
            break;
        end
        %Get Name of the interpolating function
        alpha_i_old = alpha_i;
        alpha_i=interp_func(alpha_i,alpha_max);
        phi_i_old = phi1;
        Lineiter=Lineiter+1;
    end
    return;
end

