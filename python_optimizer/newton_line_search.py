# -*- coding: utf-8 -*-
"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AUTHER:  Changbo Zhu
  E-mail:  changbozhu@outlook.com
  DATE:    Dec. 5th 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
import numpy as np
import numpy.linalg as lng
import time
from .numerdiff import numerDiff
from .line_search_model import LineSearchModel

class NewtonMethods(object):
    def __init__(self, options=None):
        self.interp_methods_avail = ['bisection']
        ####################################################################\
        #
        #    Default value of all parameters
        #    
        ####################################################################
        # norm precision or minimal norm of the gradient, 
        # it should be close to 0.
        
        self.tolgrad = 1e-5   
        # max iteration steps
        self.max_iters = 20000000
        self.verbose = True
        
        #for BFGS
        self.init_grad_norm = 1
        
        
        # Configure parameters for Linesearch and Wolfe Conidtion
        # Please reference to the file Linesearch.m for More details 
        self.c1 = 0.0001 # c1 is chosen form (0, 1) and close to 0. 
        self.c2 = 0.9    # c2 is chosen from (c1, 1).
        # Initialize minimal step length alpha0
        # or define the left point of the step length interval
        self.alpha0   = 0 # a lower boundary to step length alpha 
        # Initialize step length alpha_i which chose from (alpha0, alpha_max)
        self.alpha_i  = 1 
        # Initialize maximal step length alpha_max
        # or define the right point of the step length interval
        self.alpha_max = 120 # an upper boundary to step length alpha 
        # interpolation method in a search interval 
        # split two partitios of step length interval
        self.interp_method = 'bisection' 
        
        if options is not None:
            self.tolgrad        = options.get('tolgrad', self.tolgrad)
            self.max_iters      = options.get('max_iters', self.max_iters)
            self.verbose        = options.get('verbose', self.verbose)
            self.c1             = options.get('c1', self.c1)
            self.c2             = options.get('c2', self.c2)
            self.alpha0         = options.get('alpha0', self.alpha0)
            self.alpha_i        = options.get('alpha_i', self.alpha_i)
            self.alpha_max      = options.get('alpha_max', self.alpha_max)
            self.interp_method  = options.get('interp_method', self.interp_method)
    
    def __call__(self, obj_func, opt_vars0, obj_deriv=None, obj_hessian=None, 
                 optimizer = 'quasi-newton-bfgs'):
        if optimizer.lower() == 'quasi-newton-bfgs' or optimizer.lower() =='bfgs':
            optims = self.quasiNewtonBFGS(obj_func, opt_vars0, obj_deriv=obj_deriv)
        elif optimizer.lower() == 'newton':
            optims = self.Newton_method(obj_func, opt_vars0, obj_deriv=obj_deriv, obj_hessian=obj_hessian)
        return optims
    
    def phi_LSM(self, obj_func, xk, pk):
        return lambda alpha: obj_func(xk + alpha * pk)
    
    def dphi_LSM(self, obj_derive, xk, pk):
        return lambda alpha: np.dot(obj_derive(xk+alpha*pk),  pk)
    
    def Newton_method(self, obj_func, opt_vars0, obj_deriv=None, obj_hessian=None):
        if obj_deriv is None:
            ndiff = numerDiff()
            obj_deriv = lambda x: ndiff(obj_func, x, mode = 'gradient')
        
        if obj_hessian is None:
            ndiff = numerDiff()
            obj_hessian = lambda x: ndiff(obj_func, x, mode = 'Hessian')
            
        tstart = time.time()
        # initialize starting point
        xk   = opt_vars0
        dims = np.size(xk)
        xk = np.reshape(xk, (dims, )).astype(np.float64)
        # initialize gradient and Hessian at the starting point opt_vars0
        gradk = obj_deriv(xk)
        Hk    = obj_hessian(xk)
        
       
        #Do loop 
        k = 0
        
        while ( lng.norm(gradk) > self.tolgrad ) and ( k <= self.max_iters ):
            #print('\ngradient norm:',np.sqrt(np.dot(grad,grad)),'\n')
            if k == self.max_iters:
                print('  \n')
                print('Please notice:\n')
                print('Warning: optimization terminated because the maximum number of iterations is reached.\n')
                print('Warning: the norm of gradient at terminated point is ' + str(lng.norm(gradk) )+ ' > ' + str(self.tolgrad)  + '!')
                argMinVal = xk
                minGrad   = gradk
            pk = - lng.solve(Hk, gradk)       # Compute the descent 
            
            LSM   = LineSearchModel()
            phi   = self.phi_LSM(obj_func, xk, pk)
            dphi  = self.dphi_LSM(obj_deriv, xk, pk)
            #print('pk=%s, grad=%s, Hk=%s.'%(pk, grad, Hk))
            #print('phi(0):%s;dphi(0):%s.'%(phi(0), dphi(0)))
            alpha = LSM(phi,dphi)
              
            xk = xk + alpha*pk
            
            grad_old = gradk
            H_old    = Hk
            
            gradk = obj_deriv(xk) 
            Hk    = obj_hessian(xk)
            
            k = k + 1
        
        if np.sqrt(np.dot(gradk,gradk)) < self.tolgrad:
            argMinVal = xk
            minGrad   = gradk
            if self.verbose:
                print('\n')
                print('--Newton Method: \n')
                print('   Converged at the %sth iteration.\n' % (k))
                print('   The norm of gradient at terminated point is ' + str(lng.norm(gradk)) + ' < ' + str(self.tolgrad)  + '.\n')
          

        telapse = time.time() - tstart
        if self.verbose:
            print('\n')
            print('---Optimization Over: \n')
            print('---Elapsed time:' + str(telapse) +
                  '.\n')
   
        # Collect results
        optims = {}
        optims['argMinVal'] =  argMinVal
        optims['minVal']    =  obj_func(argMinVal)
        optims['minGrad']   =  minGrad
        
        return optims
    
    def quasiNewtonBFGS(self, obj_func, opt_vars0, obj_deriv=None):
        '''
        %quasi-Newton Method       
        % check gradient function and determin the parameter obj_deriv
        % BFGS Method Algorithm
        %     Given starting point x0, convergence tolerance epsilon >0,
        %           inverse Hessian approximation H0;
        %     iter = 0;
        %     While ||gradient(f(xk))|| > epsilon
        %           Compute search direction 
        %                                pk = - Hk * gradient(f(xk));
        %           Set xk+1 = xk + alpha_k * pk where alpha_k is computed
        %           from a line search procedure to satisfy the (strong) 
        %           Wolfe conditions ;
        %
        %           Define sk = xk+1 - xk 
        %                 and yk = gradient(f(xk+1)) - gradient(f(xk));
        %           Compute Hk+1 by means of  BFGS update equation;
        %           iter = iter + 1;
        %     end(while)
        %--------------------------------------------------------------------
        %   AUTHER:  Changbo Zhu
        %   E-mail:  changbozhu@outlook.com
        %   DATE:    Dec. 5th 2020
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        '''
        #-------------Start the quasi-Newton Method with BFGS-----------------
        # get function handle of computing gradient of the objective function 
        if obj_deriv is None:
            ndiff = numerDiff()
            obj_deriv = lambda x: ndiff(obj_func, x, mode = 'gradient')
            
        tstart = time.time()
        xk   = opt_vars0
        dims = np.size(xk)
        xk = np.reshape(xk, (dims, )).astype(np.float64)
        
        gradk = obj_deriv(xk)
        gradk = np.reshape(gradk, (dims, )).astype(np.float64)
        I    = np.eye( dims, dtype=int)
        '''
            %---------Initialize inverse Hessian approximation H0-----------        
            % the initial matrix H0 often is set to some mutiple beta*I of the
            % identity,but there is no good general strategy for choosing the 
            % multiple beta.We ask the user to prescribe  a value delta for 
            % the norm of the first step,and then set 
            %                 H0 = delta / || g0 || * I 
            % to achieve this norm.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        '''
        beta  = self.init_grad_norm/lng.norm(gradk)
        Hk = beta * I          # Generate the initial inverse Hessian 
                               # approximation
        
        #Do loop 
        k = 0
        
        while ( lng.norm(gradk) > self.tolgrad ) and ( k <= self.max_iters ):
            #print('\ngradient norm:',np.sqrt(np.dot(gradk,gradk)),'\n')
            if k == self.max_iters:
                print('  \n')
                print('Please notice:\n')
                print('Warning: optimization terminated because the maximum number of iterations is reached.\n')
                print('Warning: the norm of gradient at terminated point is ' + str(lng.norm(gradk) )+ ' > ' + str(self.tolgrad)  + '!')
                argMinVal = xk
                minGrad   = gradk
            pk = - Hk.dot(gradk)        # Compute the descent 
            
            LSM   = LineSearchModel()
            phi   = self.phi_LSM(obj_func, xk, pk)
            dphi  = self.dphi_LSM(obj_deriv, xk, pk)
            #print('pk=%s, gradk=%s, Hk=%s.'%(pk, gradk, Hk))
            #print('phi(0):%s;dphi(0):%s.'%(phi(0), dphi(0)))
            alpha = LSM(phi,dphi)
              
            xk = xk + alpha*pk
            grad_old = gradk
            gradk = obj_deriv(xk)            
                  
            # Define 
            #         sk = xk+1 - xk = alpha*pk 
            # and 
            #         yk = gradient(xk+1) - gradient(xk)
            # 
            sk   = alpha*pk 
            yk   = gradk - grad_old
                      
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # A heuristic that is often quite effective  is to scale the
            # starting matrix after the first step has been computed but
            # before the first BFGS update is performed. We change the 
            # provisional value H0 = I by setting
            #                    yk'sk
            #              H0 = ------- * I
            #                    yk'yk
            # before applying the update BFGS equation to obtain H1.
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if k ==0:
                coff = np.dot(yk, sk)/np.dot( yk, yk)
                Hk = coff*I
            '''
            % BFGS update equation:
            %                         1
            %                 rho = ------
            %                       yk'sk
            % and
            %          Hk+1 = (I - rho*sk*yk')*Hk*(I - rho*yk*sk') + rho*sk*sk'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            '''
            #print('\nbfgs:', np.dot(sk, yk))
            rho = 1.0 / (np.dot(yk, sk))            
            A1 = I - rho * sk[:, np.newaxis] * yk[np.newaxis, :]
            A2 = I - rho * yk[:, np.newaxis] * sk[np.newaxis, :]
            Hk = np.dot(A1, np.dot(Hk, A2)) + (rho * sk[:, np.newaxis] * sk[np.newaxis, :])
            k = k + 1
        
        if np.sqrt(np.dot(gradk,gradk)) < self.tolgrad:
            argMinVal = xk
            minGrad   = gradk
            if self.verbose:
                print('\n')
                print('--quasi_Newton_BFGS: \n')
                print('   Converged at the %sth iteration.\n' % (k))
                print('   The norm of gradient at terminated point is ' + str(lng.norm(gradk)) + ' < ' + str(self.tolgrad)  + '.\n')
          

        telapse = time.time() - tstart
        if self.verbose:
            print('\n')
            print('---Optimization Over: \n')
            print('---Elapsed time:' + str(telapse) +
                  '.\n')
   
        # Collect results
        optims = {}
        optims['argMinVal'] =  argMinVal
        optims['minVal']    =  obj_func(argMinVal)
        optims['minGrad']   =  minGrad
        
        return optims

