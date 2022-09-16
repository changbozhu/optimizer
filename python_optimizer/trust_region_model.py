# -*- coding: utf-8 -*-
"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AUTHER:  Changbo Zhu
  E-mail:  changbozhu@outlook.com
  DATE:    Dec. 5th 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
import numpy as np

class TrustRegionModel(object):
    '''
    Algorithm 4.1 (Trust Region)：
        Given hat_Delta >= 0, Delta0 in (0, hat_Delta), and eta in [0, 1/4):
        for k = 0,1,2,3,...
            Obtain pk by (approximately) solving 
                    min_{p in R^n} m_k(p) = f_k + g_k^T p + 1/2 p^T H_k p
                    s.t. ||p|| <= Delta_k
            
            Evaluate rho_k:
            rho_k = [f(x_k) - f(x_k + p_k)]/[ m_k(0) - m_k(p_k) ]
            
            if rho_k < 1/4
                Delta_{k+1}  = 1/4 Delta_k
            else
                if rho_k > 3/4 and ||p_k|| = Delta_k
                    Delta_{k+1} = min{2Delta_k, hat_Delta }
                else
                    Delta_{k+1} = Delta_k
            
            if rho_k > eta
                x_{k+1} = x_k + p_k
            else
                x_{k+1} = x_{k}
        end(for)
    '''
    def __init__(self, options = None):
        # Configure parameters for Trust Region
        # Initialize maximum radius of trust region ball
        # an overall bound on the step lengths
        self.hat_Delta  = 12 
        # Initialize radius of trust region ball
        self.Delta0     = 1/2*self.hat_Delta
        # maximum iterations
        self.max_iters  = 200
        self.eta        = None
        self.verbose    = True
        
        if options is not None:
            self.verbose        = options.get('verbose', self.verbose)
            self.hat_Delta      = options.get('hat_Delta', self.hat_Delta)
            self.Delta0         = options.get('Delta0', self.Delta0)
            self.max_iters      = options.get('max_iters', self.max_iters)
            self.eta            = options.get('eta',self.eta)
    
    def trust_region(self, m):
        '''
            Trust Region Algorithm
            f(x_k + p) ~ m_k(p)
            B_k ~ Hessian(x_k)
            m_k(p) = f(x_k)  +  nabla f(x_k)^T * p  + 1/2 p^T*B_k*p
        '''
        max_iters = self.max_iters
        # Initialize some corresponding parameters
        if self.eta is None:
            eta = np.random.rand()/4
        else:
            eta = self.eta
        
        
        Delta_k = self.Delta0
        for k in range(max_iters):
            pk = self.solving_subproblem(Delta_k)
            rho_k = self.comput_ratio(pk)
            
        #end(for)
            
        return alpha
    def quadra_func_m(self, obj_func, g, B, x):
        d = np.size(p)
        m = lambda p: obj_func(x) + np.dot(g, p) + 1/2*np.matmul(np.reshape(p, [d,]), np.matmul(B, p) )        
    def comput_ratio(self, pk):
        obj_func(x) + 
    def solving_subproblem(self, Delta_k):
        pass