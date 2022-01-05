
% Purpose: Solve one dimensional Poisson equation with
%          a) nonhomegeneous DBC both sides
%          b) NBC at left and nonhomegeneous DBC at right
%
% Author : Hamdullah Yuecel
%
% Input  : ---
%
% Output : --
%
% Usage : drive()


function drive()

% Clear every thing so it won't mess up with other existing variables.
close all;  
clc

% Initilize 

errND_L2    = zeros(7,1);
errND_inf   = zeros(7,1);
errND_H1    = zeros(7,1);
ratioND_L2  = zeros(6,1);
ratioND_inf = zeros(6,1);
ratioND_H1  = zeros(6,1);

N = zeros(7,1);

for i=1:7

  % Generate a mesh
  x   = linspace(-1,2,2^(i+1)+1);
  N(i)= 2^(i+1);
  
  
  % Get the FEM solution at nodal points for the ND case
  UND = fem1d_ND(x,@source,@u_exact,@p_exact,@c_exact,@q_exact);
  
  
  % Using the FEM coefficients evaluate at different points
  x2 = linspace(-1,2,1000);
  k  = length(x2);
  
  % Initialize  
  uND  = zeros(k,1);
  uNDx = zeros(k,1);
  
  for j=1:k
   [uND(j),uNDx(j)] = fem_soln(x,UND,x2(j)); 
  end
     
  % Compute the analytic solution
  [uex,uex_x] = u_exact(x2);
  
  % Compute the errors
  errND_L2(i)  = norm(uex'-uND);
  errND_H1(i)  = sqrt(errND_L2(i)^2 + norm(uex_x'-uNDx)^2);  
  errND_inf(i) = max(abs(uex'-uND));
  
end


% Display the results
fprintf('Neumann-Dirichlet case:\n\n')
fprintf('   N      L2-error   rate   H1-error   rate    Inf-error   rate\n')

for i=1:7
   if(i>1)
    % Compute ratios  
    ratioND_L2(i-1)  = log2(errND_L2(i-1)/errND_L2(i))/log2(2);
    ratioND_inf(i-1) = log2(errND_inf(i-1)/errND_inf(i))/log2(2);
    ratioND_H1(i-1)  = log2(errND_H1(i-1)/errND_H1(i))/log2(2);
    
    fprintf('%4d      %3.2e   %4.2f   %3.2e   %4.2f    %3.2e    %4.2f\n', ...
            N(i),errND_L2(i),ratioND_L2(i-1),errND_H1(i),ratioND_H1(i-1), errND_inf(i), ratioND_inf(i-1));
  
   else     
     r='-';    
     fprintf('%4d      %3.2e    %c     %3.2e    %c      %3.2e     %c\n', ...
            N(i),errND_L2(i),r,errND_H1(i),r, errND_inf(i),r); 
  
  end 
end

% Draw the figures
figure(1);
loglog(N,errND_L2,  '-ro','LineWidth',2);  hold on;
loglog(N,errND_H1,  '--s','LineWidth',2);  hold on;
loglog(N,errND_inf, '-.kh','LineWidth',2);
legend('L^2','H^1','L^\infty')
title('Neumann-Dirichlet','FontSize',18)
xlabel('N','FontSize',14);
ylabel('Error','FontSize',14);

end

%% Define functions
% Exact solution 
function [u,ux] = u_exact(x)
  u  = sin(x).*exp(x);
  ux = exp(x).*(cos(x) + sin(x));
end

function p = p_exact(x)
  p = 1 + x.^2;
end

function c = c_exact(~)
  c = 1;
end

function q = q_exact(x)
  q = exp(-x);
end

% Source function
function y = source(x)
  y = exp(x).*(-2*x.^2 -2*x -1).* cos(x) + ( exp(x).*(1-2*x) + 1).*sin(x); 
end


