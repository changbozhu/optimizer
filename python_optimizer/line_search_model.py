# -*- coding: utf-8 -*-
"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AUTHER:  Changbo Zhu
  E-mail:  changbozhu@outlook.com
  DATE:    Dec. 5th 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
# information
version = 'OPT-LSM-V2020.12'
author  = 'Changbo Zhu'
address = 'Beijing, China'
email   = 'changbozhu@outlook.com'
# information 

from . import utilities as utils 

class LineSearchModel(object):
    '''
        We now consider techniques for finding a minimum of the one-dimensional function
                phi(a) = f(xk + a*p_k)
        or for simply finding a step length a_k satisfying Wolfe Condition.
        
    '''
    def __init__(self, options = None):
        self.mode = 'Wolfe'
        self.interp_methods_avail = [ 'bisection' ]
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
            self.verbose        = options.get('verbose', self.verbose)
            self.c1             = options.get('c1', self.c1)
            self.c2             = options.get('c2', self.c2)
            self.alpha0         = options.get('alpha0', self.alpha0)
            self.alpha_i        = options.get('alpha_i', self.alpha_i)
            self.alpha_max      = options.get('alpha_max', self.alpha_max)
            self.interp_method  = options.get('interp_method', self.interp_method)
        
        if self.interp_method in self.interp_methods_avail:
            if self.interp_method=='bisection':
                self.interp_func = utils.bisection
    
    
    def __call__(self, phi, dphi):
        return self.linesearch_Wolfe(phi, dphi)
    
    
    def linesearch_Wolfe(self, phi, dphi):
        '''
        Algorithm: Line Search Algorithm With the Wolfe Conditions
            set alpha_0 <--0, 
            choose alpha_max > 0 ,
            and 0 < alpha_1 < alpha_max;
            i <--1
            phi_0 <--phi(alpha_0) 
            d_phi_0 <-- phi'(alpha_0)
            repeat
                phi_i <-- phi(alpha_i)
                if phi_i> phi_0 + c1*alpha_i*d_phi_0 or [ phi_i >= phi_{i-1}  and i >1 ]
                    alpha_star <--zoom(alpha_i, alpha_{i-1})  and stop
                
                d_phi_i <-- phi'(alpha_i)
                
                if |d_phi_i| <= -c2*d_phi_0
                    set alpha_star <---alpha_i and stop
                
                Choose alpha_{i+1} in (alpha_i, alpha_{i-1})
                i <-- i+1           
            end(repeat)
        '''
        
        # Initialize some corresponding parameters
        max_iters = 5000
        c1 = self.c1
        c2 = self.c2
        
        # Initialize step length alpha0  and alpha_max
        # Define the left point of step length interval
        alpha0 =0 # self.alpha0
        alpha_i = self.alpha_i
        alpha_max = self.alpha_max        
        
        phi_0 = phi(0)
        dphi_0 = dphi(0)
        #print('phi0=%s, dphi0=%s'%(phi_0, dphi_0) )
        
        alpha_i_old = alpha0
        phi_i_old = phi_0     
        
        # Do loop to search a step length 
        termin = False
        iters  = 1
        while (not termin) and (iters < max_iters):
            phi_i = phi(alpha_i)
            Gline = phi_0 + c1* alpha_i * dphi_0
            #print('linesearch:', iters, 'Gline:',Gline, 'phi_i:', phi_i, 'phi_i_old:', phi_i_old, 'dphi_i:', dphi_i )
            if (phi_i > Gline) or ( phi_i >= phi_i_old  and iters > 1):
                alpha = self.zoom(alpha_i_old, alpha_i, phi, dphi)
                termin = True
                break
            #endif
            
            dphi_i = dphi(alpha_i)
            if abs(dphi_i) <= -c2*dphi_0:
                alpha = alpha_i
                termin = True
                break
            #endif
            
            if dphi_i >= 0:
                alpha = self.zoom(alpha_i, alpha_i_old, phi, dphi)
                termin = True
                break
            #endif
            alpha_i_old = alpha_i
            phi_i_old = phi_i
            alpha_i = self.interp_func(alpha_i, alpha_max)
            
            iters += 1
        # end loop 
            
        if termin:
            return alpha
        else:
            assert 0, 'No step length alpha satisfying Wolfe Condition!'
    
    
    
    def zoom(self, alpha_lo, alpha_hi, phi, dphi):
        ''' zoom function is compatible with Line Search.
        repeat    
            Interpolate (using quadratic, cubic or bisection) to find a trial 
            step length alpha_j between alpha_lo and alpha_hi;
            
            phi_j <-- phi(alpha_j)
            if phi_j > phi_0 + c1*alpha_j*d_phi_0  or phi_j >= phi(alpha_lo)
                alpha_hi <-- alpha_j;
            else
                d_phi_j = phi'(alpha_j)
                if |d_phi_j| <= -c2*d_phi_0
                    set alpha_star  <--- alpha_j and stop;
                
                if d_phi_j*(alpha_hi - alpha_lo) >= 0
                    alpha_hi <---alpha_lo;
                alpha_lo <-- alpha_j;
        end(repeat)    
        '''
        max_iters = 10
        c1 = self.c1
        c2 = self.c2
        alpha0 = self.alpha0
        phi_0 = phi(alpha0)
        dphi_0 = dphi(alpha0)
        #print('zoom1:','alpha_lo:',alpha_lo, 'alpha_hi', alpha_hi)
        
        termin = False
        iters = 1
        while  iters <= max_iters:
            #print('zoom:', iters)
            #Interpolate (using quadratic,cubic or bisection) 
            #to find a trial steplength alpha_j 
            #between alpha_lo and alpha_hi;
            alpha_j = self.interp_func(alpha_lo, alpha_hi)
            
            #print(alpha_j)
            phi_j = phi(alpha_j)
            Gline = phi_0 + c1*alpha_j*dphi_0
            if (phi_j > Gline) or (phi_j >= phi(alpha_lo)):
                alpha_hi = alpha_j
            else:
                dphi_j = dphi(alpha_j)
                
                if abs(dphi_j) <= -c2*dphi_0:
                    alpha = alpha_j
                    termin = True
                    break
                if dphi_j*(alpha_hi - alpha_lo) >= 0:
                    alpha_hi = alpha_lo
                alpha_lo = alpha_j
                
            iters = iters + 1
        
        if termin:
            return alpha
        else:
            assert 0, 'No step length alpha satisfied to Wolfe Condition!'