# -*- coding: utf-8 -*-
"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AUTHER:  Changbo Zhu
  E-mail:  changbozhu@outlook.com
  DATE:    Dec. 3th 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
import numpy as np
import sys

class numerDiff(object):
    def __init__(self, options = None):
        '''
          options:
        %    Optionally, the third argument of the function can be a structure containing further
        %    settings for Ridder's method.
        %
        %    options['init_eps']      Initial finite difference (default: 1)
        %    options['ratio']         Divisor used to reduce h on each step (default: 1.2)
        %    options['min_steps']     Minimum number of steps in h (default: 3)
        %    options['max_steps']     Maximum number of steps in h (default: 100)
        %    options['threshold']     Terminate if last step worse than preceding by a factor of threshold
        %                             (default: 1.5)
        '''
        self.init_eps        = 1
        self.ratio           = 0.8
        self.min_steps       = 5
        self.max_steps       = 100
        
        self.threshold       = 1.5
        self.min_err         = sys.maxsize
        self.tolerance       = 1e-20
        
        if options is not None:            
            self.init_eps     = options.get('init_eps',  self.init_eps)
            self.ratio        = options.get('ratio',     self.ratio)
            self.min_steps    = options.get('min_steps', self.min_steps)
            self.max_steps    = options.get('max_steps', self.max_steps)
            self.threshold    = options.get('threshold', self.threshold)
            self.tolerance    = options.get('tolerance', self.tolerance)
    
    def __call__(self, obj_func, x, mode='gradient'):
        if mode.lower() in ['gradient', 'deriv1']:
            result, err = self.evalgradient(obj_func, x)
        elif mode.lower() in ['hessian', 'deriv2']:
            result, err = self.evalHessian(obj_func, x)
        return result
    
    def evalscarlarderiv(self, phi):
        ''' 
          Evaluates derivative of the function phi at point 0 according to Ridders' method:
        
          Ridders, CJF. (1982). Accurate computation of F'(x) and F'(x) F''(x). 
          Advances in Engineering Software, 4(2), 75-6.
         
          INPUT:
            phi           Function handle of a scalar real function of one real variable
        
          OUTPUT:
             deriv         First derivative of phi at 0
             min_err       Error estimate
        '''
        #---------------------------------------------------------------------
        # Initialize recording matrix A by max_steps times max_steps
        A = np.full([self.max_steps, self.max_steps], np.nan)
        # |------|--------------|---------------|---------------|--------------
        # |      |O(eps^(2*0+2))|O(eps^(2*1+2)) |O(eps^(2*2+2)) |O(eps^(2*n+2))
        # |------|--------------|---------------|---------------|--------------
        # | m=0: |  A0(eps)     |               |               |
        # | m=1: |  A0(r*eps)   |   A1(eps)     |               |
        # | m=2: |  A0(r^2*eps) |   A1(r*eps)   |   A2(eps)     |
        # | m=3: |  A0(r^3*eps) |   A1(r^2*eps) |   A2(r*eps)   |
        # |  :   |      :       |          :    |       :       |
        # |  :   |      :       |          :    |       :       |
        # | m=M; | A0(r^M*eps)  |A1(r^(M-1)*eps)|A2(r^(M-2)*eps)|
        # |------|--------------|---------------|---------------|--------------
        #  N = M = max_steps - 1  
        #  A(m,0) = (phi(r^m *eps) - phi(- r^m * eps))/(2*r^m *eps)
        #  A(m,n) = (A(m,n-1) -  r^2*A(m-1 ,n-1))/(1 - r^2), 0 < n < m
        #  absolutely_error(i,j) 
        #  = max(abs(A(m,n) -  A( m,n-1)), abs(A(m,n) -  A(m-1,n-1)))/(1 - r^2)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize finite difference step epsilon (eps)
        eps = self.init_eps
        min_err = self.min_err
        M   = self.max_steps - 1
        N   = self.max_steps - 1
        deriv = np.nan
        # Compute central difference formula at initial step
        # A[0,0] = ( phi(eps) - phi(-eps) )/(2*eps)
        A[0,0] = phi(eps)
        # Loop to reduce eps 
        for m in range(1, M+1):
            # Generate a new step size eps_m = ratio^m * init_eps
            eps = self.ratio*eps
            # Compute 2nd order approximation by central difference formula at 
            # this step 
            # A[m,0] = ( phi(eps) - phi(-eps) )/(2*eps)
            A[m,0] = phi(eps)
        
            # Use square of ratio for extrapolation because errors increase
            # quadratically with eps (here, of course, they decrease quadratically
            # because we're reducing eps ...)
            ratio_sq = self.ratio**2
            coeff = ratio_sq
            # Fill the current row using Richardson extrapolation
            for n in range(1, N+1):
            # Richardson extrapolation
                A[m,n] = ( A[m, n-1] - coeff*A[m-1, n-1] )/(1 - coeff)
                
                # Increment extrapolation factor
                coeff = coeff*ratio_sq
            
                # Error on this trial is defined as the maximum absolute difference
                # to the extrapolation parents
                err_mn = np.amax([np.abs(A[m,n] - A[m,n-1]), np.abs(A[m,n] - A[m-1,n-1])]);
                
                if err_mn < min_err:
                    min_err = err_mn
                    deriv = A[m,n]
            # Stop if errors start increasing (to be expected for very small
            # values of eps
            if m > (self.min_steps - 1) and abs( A[m,n] - A[m-1,n-1] ) > (self.threshold*min_err) and min_err < self.tolerance:
                return  (deriv,  min_err)
        return (deriv,  min_err)
    
    def phi_1(self, obj_func, x, p, eps):
        x = np.reshape(x, (np.size(x)  ,))
        p = np.reshape(p, (np.size(p)  ,))
        return obj_func(x + eps*p )
    
    def phi_1_central_diff_formula(self, obj_func, x, p, eps):
        x = np.reshape(x, (np.size(x), ))
        p = np.reshape(p, (np.size(p), ))
        ydiff = self.phi_1(obj_func, x, p, eps) - self.phi_1(obj_func, x, p, -eps)
        diff_phi_1 = ydiff/(2*eps)
        return diff_phi_1
    
    def evalgradient(self, obj_func, x):
        d = np.size(x)
        I = np.eye(d)
        grad = np.full([d, ], np.nan)
        min_err =  sys.maxsize * np.ones([d, ])
        # Loop through each component of variable x
        for i in range(d):
            ei = I[:,i]
            # Construct filehandle to be passed to eval
            phi = lambda eps: self.phi_1_central_diff_formula(obj_func, x, ei, eps)
            
            # Calculate derivative
            [grad[i], min_err[i]] = self.evalscarlarderiv(phi)
        return grad, min_err
    
    def phi2(self, eps, obj_func, x, ei, ej):
        d  = np.size(x)
        x  = np.reshape(x,  (d, ) )
        ei = np.reshape(ei, (d, ) )
        ej = np.reshape(ej, (d, ) )
        
        diff2 = obj_func(x + eps*ei + eps*ej ) - obj_func(x - eps*ei + eps*ej ) - obj_func(x + eps*ei - eps*ej ) + obj_func(x - eps*ei - eps*ej )
        diff2 = diff2/(4*eps**2)
        return diff2
    
    def evalderiv2(self, phi2):
        deriv2 = np.nan
        
        # Initialize recording matrix A by max_steps times max_steps
        A = np.full([self.max_steps, self.max_steps], np.nan)
        '''
        |------|--------------|---------------|---------------|--------------
        |      |O(eps^(2*0+2))|O(eps^(2*1+2)) |O(eps^(2*2+2)) |O(eps^(2*n+2))
        |------|--------------|---------------|---------------|--------------
        | m=0: |  A0(eps)     |               |               |
        | m=1: |  A0(r*eps)   |   A1(eps)     |               |
        | m=2: |  A0(r^2*eps) |   A1(r*eps)   |   A2(eps)     |
        | m=3: |  A0(r^3*eps) |   A1(r^2*eps) |   A2(r*eps)   |
        |  :   |      :       |          :    |       :       |
        |  :   |      :       |          :    |       :       |
        | m=M; | A0(r^M*eps)  |A1(r^(M-1)*eps)|A2(r^(M-2)*eps)|
        |------|--------------|---------------|---------------|--------------
        N = M = max_steps - 1  
        A(m,0) = (phi(r^m *eps) - phi(- r^m * eps))/(2*r^m *eps)
        A(m,n) = (A(m,n-1) -  coeff(k)*A(m-1 ,n-1))/(1 - coeff(k)), 0 < n < m
        coeff(k) = ratio^2 * coeff(k-1), k = m*(m-1)/2 + n
        absolutely_error(i,j) 
        = max(abs(A(m,n) -  A( m,n-1)), abs(A(m,n) -  A(m-1,n-1)))/(1 - r^2)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        '''
        # Initialize finite difference step epsilon (eps)
        eps = self.init_eps
        min_err = self.min_err
        M = self.max_steps - 1
        N = self.max_steps - 1
        # Compute central difference formula at initial step
        A[0,0] = phi2(eps)
    
        # Loop to reduce eps 
        for m in range(1, M+1):
            n=0
            # Generate a new step size eps_m = ratio^m * init_eps
            eps = self.ratio*eps
        
            # Compute 2nd order approximation by central difference formula at 
            # this step 
            A[m,n] = phi2(eps)
        
            # Use square of ratio for extrapolation because errors increase
            # quadratically with eps (here, of course, they decrease quadratically
            # because we're reducing eps ...)
            ratio_sq = self.ratio**2
            coeff = ratio_sq
            # Fill the current row using Richardson extrapolation
            for n in range(1,N+1):
                # Richardson extrapolation
                A[m,n] = ( A[m,n-1] - coeff*A[m-1,n-1] )/(1 - coeff)
            
                # Increment extrapolation factor
                coeff = coeff*ratio_sq
            
                # Error on this trial is defined as the maximum absolute difference
                # to the extrapolation parents
                err_mn = max(abs(A[m,n]-A[m,n-1]),abs(A[m,n] - A[m-1,n-1]))
            
                if err_mn < min_err:
                    min_err = err_mn
                    deriv2 = A[m,n]

            # Stop if errors start increasing (to be expected for very small
            # values of eps
            if m > (self.min_steps - 1) and abs( A[m,n] - A[m-1,n-1] ) > (self.threshold*min_err)\
                and min_err < self.tolerance:
                return  deriv2,  min_err
        return  deriv2,  min_err
    def evalHessian(self, obj_func, x):
        '''
        % Calculates the Hessian (i.e., d^2f/(dx_idx_j)) of the function f at point x
        % according to Ridders' method:
        %
        % Ridders, CJF. (1982). Accurate computation of F'(x) and F'(x) F''(x). 
        % Advances in Engineering Software, 4(2), 75-6.
        %
        % INPUT:
        %    obj_func      Function handle of a scalar real function of a real vector variable
        %                  which are passed as *one* vector with *d* elements
        %      x           Point at which to differentiate f
        %
        % OUTPUT:
        %    Hessian       a set of 2nd derivatives of a given function f with
        %                  respect to a real vector x.
        %    min_err       Error estimate, a matrix with the same shape of Hessian
        %
        % 
        '''
    
        d = np.size(x)
        I = np.eye(d, dtype=np.int)
        # results
        Hessian  = np.full((d,d), np.nan)
        min_err  = self.min_err*np.ones([d,d])
    
        # Biuld struct of configuration for evaluating Hessian of function obj_func.    
        for i in range(d):
            ei = I[i,:]
            for j  in range(i+1):
                ej = I[:,j]
                phi2 = lambda eps: self.phi2(eps, obj_func, x, ei, ej)
                [deriv2_ei_ej , min_err_ei_ej] = self.evalderiv2(phi2)
                Hessian[i,j] = deriv2_ei_ej
                min_err[i,j] = min_err_ei_ej
                if j<i:
                    Hessian[j,i] = deriv2_ei_ej
                    min_err[j,i] = min_err_ei_ej
        return Hessian, min_err
    


if __name__== '__main__':
    def rosenbrock(x):
        return (1-x[0])**2 + 100*(x[1]-x[0]**2)**2

    def test_func(x):
        x = np.reshape(x, (np.size(x), ))
        P = np.identity(np.size(x))
        P[0,0] = 4
        return np.matmul(np.matmul(np.reshape(x, (1, np.size(x)) ), P), x) 
    
    #print('test_fun:',test_fun([1,1]))
    numDiff = numerDiff()
    #print('[1,1]grad:',numDiff.phi_1_central_diff_formula(test_fun, [1,1], [1,0], 0.0001))
    #print('[1,1]grad:',numDiff.phi_1_central_diff_formula(test_fun, [1,1], [0,1], 0.0001))
    print('[1,1]grad:',numDiff.evalgradient(test_func, [-100,100]))
    #phi2 =lambda  eps:numDiff.phi2(eps, test_fun, [1,1],[0,1],[0,1])
    print('[1,1]grad:',numDiff(test_func, [-100,100]))
    #print(nDiff.evalscarlarderiv(test_fun))
    