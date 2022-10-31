    format bank;
    clear;
    clc;
    % configure parameters of Newton Method    
    options = quasi_Newton_BFGS_options_rosenbrock();
    objpara.cofficient=[1 ;10];
    objpara.rosenbrock = '2D';
    obj_func = @(opt_var) rosenbrock( opt_var, objpara ); 
    obj = struct;
    obj.obj_func = obj_func;
    obj_grad = @(opt_var) rosenbrock_gradient(opt_var, objpara );    
    %obj.obj_grad = obj_grad;    
    opt_var0=[0; 0];
    [argMinVal, optims]= quasi_Newton_BFGS(obj, opt_var0, options);