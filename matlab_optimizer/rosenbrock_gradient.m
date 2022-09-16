function [ out_g ] = rosenbrock_gradient( xk , objpara )
% In mathematical optimization, the Rosenbrock function is a non-convex function 
% used as a performance test problem for optimization algorithms introduced by 
% Howard H. Rosenbrock in 1960.[1] It is also known as Rosenbrock's valley or 
% Rosenbrock's banana function.
% The global minimum is inside a long, narrow, parabolic shaped flat valley. 
% To find the valley is trivial. To converge to the global minimum, however,
% is difficult.
% The function is defined by
%             f(x,y) = (a - x)^2 + b*(y - x^2)^2
% It has a global minimum at  (x,y)=(a,a^2), where  f(x,y)=0. 
% Usually a=1 and b=100. 
% Only in the trivial case of a=0 is the function symmetric and the minimum at origo.
% Two variants are commonly encountered. One is the sum of N/2 uncoupled 2D Rosenbrock problems,
%  f(X) = f(x(1),x(2),x(3),...,x(N)) = sum(i=1,N/2){a*( x(2i-1) - x(2i) )^2 + ( x(2i-1) - 1 )^2 }
% This variant has been shown to have exactly one minimum for  N=3 (at [1, 1, 1]) 
% and exactly two minima for  4<= N <= 7 ¡ªthe global minimum of all ones and 
% a local minimum near [x(1),x(2),x(3),...,x(N)] = [-1,1,1,...,1]. 
% This result is obtained by setting the gradient of the function equal to zero,
% noticing that the resulting equation is a rational function of x. 
% For small  N the polynomials can be determined exactly and Sturm's theorem 
% can be used to determine the number of real roots, while the roots can 
% be bounded in the region of |x(i)|<2.4. For larger  N this method breaks down
% due to the size of the coefficients involved.
rosenModel = objpara.rosenbrock;
coff = objpara.cofficient;
a = coff(1);
b = coff(2);
if  strcmp(rosenModel, '2D')
    x = xk(1);
    y = xk(2);
    out_g = [-2*(a - x) - 4*b*x*(y - x^2);2*b*(y - x^2)];
    
elseif strcmp(rosenModel, 'uncoupled_2D')
   N = length(xk);
    if rem(N ,2) ==0
        
         m = N/2;
         out_g = zeros(N,1);
         
         for i=1:1:m
             out_g(2*i-1) = 4*a*xk(2*i-1) * ( xk(2*i-1)^2 - xk(2*i) ) + 2*b*( xk(2*i-1) - 1 );
             out_g(2*i  ) = -2*a*( xk(2*i-1)^2 - xk(2*i) );
         end   
         
    else
        display('Error:Dimensions are error!');
    end
    
elseif  strcmp(rosenModel,'coupled_2D')
     N = length(xk);
     out_g = zeros(N,1);
     out_g(1) = -4*a*xk(1)*(xk(2) - xk(1)^2 ) - 2*b*(1 - xk(1));
     for i=2:1:N-1
         out_g(i) = -4*a*xk(i)*( xk(i+1) - xk(i)^2 )  - 2*b*( 1 - xk(i) ) + 2*a*(xk(i) - xk(i-1)^2 );
     end  
     out_g(N) = 2*a*(xk(N) - xk(N-1)^2);
else
     display('Error:Parameter is incorrect!');
end
return;  
end