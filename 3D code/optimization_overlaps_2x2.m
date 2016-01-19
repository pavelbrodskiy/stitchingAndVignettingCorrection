% This function returns the fitness of the solution for a combination of
% parameters. It should not modify any global variables, but reads E from
% the optimization wrapper and dimensions from the recursive algorithm
% call. SA is the solution matrix in the form: [a11 a12 b11 b12].

% SA: a, then b in form [11 21 12 22]

function [ SSR ] = optimization_overlaps_2x2( SA )
       
    global E11 E12 E13 E14 E21 E22 E23 E24;
    
    SSR1 = (E11 * SA(1) - E21 * SA(2) - SA(5) + SA(6)).^2;
    SSR2 = (E12 * SA(2) - E22 * SA(4) - SA(6) + SA(8)).^2;
    SSR3 = (E13 * SA(4) - E23 * SA(3) - SA(8) + SA(7)).^2;
    SSR4 = (E14 * SA(3) - E24 * SA(1) - SA(7) + SA(5)).^2;
    
    SSR = sum(SSR1(:))+sum(SSR2(:))+sum(SSR3(:))+sum(SSR4(:));
    
end