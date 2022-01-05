
function U = fem1d_ND(x,source,fexact,pex,cex,qex)        
                                                
%  Simple Matlab code for 1D FEM for             
%                                                  
%     -(pu')' + c u' + qu= f(x),   a <= x <= b,  
% 
%   with Neumann (left) and nonhomegeneous Dirichlet BC (right)
%
%  Input:  x,      Nodal points      
%          source, Source function
%          fexact, Exact solution
%          pex,    function p
%          cex,    function c
%          qex,    function q
%
%  Output: U, FEM solution at nodal points                         
%                                                                      
                                    
M = length(x);

% Compute step sizes
h = zeros(M,1);
for i=1:M-1,
  h(i) = x(i+1)-x(i);
end
 
A = sparse(M,M);   % Initialize stiffness matrix
F = zeros(M,1);    % Initialize load vector

for i=1:M-1
  
  [p11,p12,p21,p22] = int_hat_p(x(i),x(i+1),pex);  % contribution from p part
  [c11,c12,c21,c22] = int_hat_c(x(i),x(i+1),cex);  % contribution from c part
  [q11,q12,q21,q22] = int_hat_q(x(i),x(i+1),qex);  % contribution from q part   
    
  A(i,i)     = A(i,i)     +  p11 + c11 + q11;
  A(i,i+1)   = A(i,i+1)   +  p12 + c12 + q12;                                   
  A(i+1,i)   = A(i+1,i)   +  p21 + c21 + q21;                  
  A(i+1,i+1) = A(i+1,i+1) +  p22 + c22 + q22; 
  
  [f1,f2] = int_hat_f(x(i),x(i+1),source);         % right-hand side
               
  F(i)       = F(i)   + f1 ; 
  F(i+1)     = F(i+1) + f2 ;
  
end

% Add Neumann BC 
[~,ux] = feval(fexact,-1);
F(1)   = F(1) - ux*feval(pex,-1);

% Dirichet Solutions
UD    = zeros(M,1);
UD(M) = feval(fexact,2);    % right node

% Extract Dirichlet solutions from rhs
F = F - A*UD;

% Evaluate at Dirichlet points
A(:,M) = 0 ; A(M,:) = 0;  A(M,M) = 1; 

F(M) = feval(fexact,2);    % right node

% Solve linear system
U = A\F;

end


% Contributions to stiffnes matrix coming from p part 
function [p11,p12,p21,p22] = int_hat_p(x1,x2,pex)

xm = (x1+x2)*0.5;

p1  = feval(pex,x1);
pxm = feval(pex,xm);
p2  = feval(pex,x2);

p11 = (x2-x1)*( p1*hatd1(x1,x1,x2)*hatd1(x1,x1,x2) + 4*pxm*hatd1(xm,x1,x2)*hatd1(xm,x1,x2) + p2*hatd1(x2,x1,x2)*hatd1(x2,x1,x2) )/6;
p12 = (x2-x1)*( p1*hatd1(x1,x1,x2)*hatd2(x1,x1,x2) + 4*pxm*hatd1(xm,x1,x2)*hatd2(xm,x1,x2) + p2*hatd1(x2,x1,x2)*hatd2(x2,x1,x2) )/6;
p21 = (x2-x1)*( p1*hatd2(x1,x1,x2)*hatd1(x1,x1,x2) + 4*pxm*hatd2(xm,x1,x2)*hatd1(xm,x1,x2) + p2*hatd2(x2,x1,x2)*hatd1(x2,x1,x2) )/6;
p22 = (x2-x1)*( p1*hatd2(x1,x1,x2)*hatd2(x1,x1,x2) + 4*pxm*hatd2(xm,x1,x2)*hatd2(xm,x1,x2) + p2*hatd2(x2,x1,x2)*hatd2(x2,x1,x2) )/6;

end 

% Contributions to stiffnes matrix coming from c part 
function [c11,c12,c21,c22] = int_hat_c(x1,x2,cex)

xm = (x1+x2)*0.5;

c1  = feval(cex,x1);
cxm = feval(cex,xm);
c2  = feval(cex,x2);

c11 = (x2-x1)*( c1*hat1(x1,x1,x2)*hatd1(x1,x1,x2) + 4*cxm*hat1(xm,x1,x2)*hatd1(xm,x1,x2) + c2*hat1(x2,x1,x2)*hatd1(x2,x1,x2) )/6;
c12 = (x2-x1)*( c1*hat1(x1,x1,x2)*hatd2(x1,x1,x2) + 4*cxm*hat1(xm,x1,x2)*hatd2(xm,x1,x2) + c2*hat1(x2,x1,x2)*hatd2(x2,x1,x2) )/6;
c21 = (x2-x1)*( c1*hat2(x1,x1,x2)*hatd1(x1,x1,x2) + 4*cxm*hat2(xm,x1,x2)*hatd1(xm,x1,x2) + c2*hat2(x2,x1,x2)*hatd1(x2,x1,x2) )/6;
c22 = (x2-x1)*( c1*hat2(x1,x1,x2)*hatd2(x1,x1,x2) + 4*cxm*hat2(xm,x1,x2)*hatd2(xm,x1,x2) + c2*hat2(x2,x1,x2)*hatd2(x2,x1,x2) )/6;

end 

% Contributions to stiffnes matrix coming from q part 
function [q11,q12,q21,q22] = int_hat_q(x1,x2,qex)

xm = (x1+x2)*0.5;

q1  = feval(qex,x1);
qxm = feval(qex,xm);
q2  = feval(qex,x2);

q11 = (x2-x1)*( q1*hat1(x1,x1,x2)*hat1(x1,x1,x2) + 4*qxm*hat1(xm,x1,x2)*hat1(xm,x1,x2) + q2*hat1(x2,x1,x2)*hat1(x2,x1,x2) )/6;
q12 = (x2-x1)*( q1*hat1(x1,x1,x2)*hat2(x1,x1,x2) + 4*qxm*hat1(xm,x1,x2)*hat2(xm,x1,x2) + q2*hat1(x2,x1,x2)*hat2(x2,x1,x2) )/6;
q21 = (x2-x1)*( q1*hat2(x1,x1,x2)*hat1(x1,x1,x2) + 4*qxm*hat2(xm,x1,x2)*hat1(xm,x1,x2) + q2*hat2(x2,x1,x2)*hat1(x2,x1,x2) )/6;
q22 = (x2-x1)*( q1*hat2(x1,x1,x2)*hat2(x1,x1,x2) + 4*qxm*hat2(xm,x1,x2)*hat2(xm,x1,x2) + q2*hat2(x2,x1,x2)*hat2(x2,x1,x2) )/6;

end 

% Contributions to load vector
function [f1,f2] = int_hat_f(x1,x2,source)

xm = (x1+x2)*0.5;

f1  = feval(source,x1);
fxm = feval(source,xm);
f2  = feval(source,x2);

f1  = (x2-x1)*( f1*hat1(x1,x1,x2) + 4*fxm*hat1(xm,x1,x2) + f2*hat1(x2,x1,x2) )/6;
f2  = (x2-x1)*( f1*hat2(x1,x1,x2) + 4*fxm*hat2(xm,x1,x2) + f2*hat2(x2,x1,x2) )/6;

end
