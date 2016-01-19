% This function returns the fitness of the solution for a combination of
% parameters. It should not modify any global variables, but reads E from
% the optimization wrapper and dimensions from the recursive algorithm
% call. SA is the solution matrix in the form: [a11 a12 b11 b12].

% SA: [a1/a2 b1-b2]

function [ SSR ] = optimization_overlaps( SA )
       
    global E1 E2;
    
    SSR = (E1 * SA(1) - E2 * SA(2) - SA(3) + SA(4)).^2;
    SSR = sum(SSR(:));
    
end