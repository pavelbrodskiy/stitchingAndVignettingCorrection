% Hand this function an M * N * x * y * 4 set of edge intensities (E) and
% the function will make several decisions: if the matrix is sufficiently
% small (2x1, 1x2, 1x1) it will directly call the optimization wrapper and
% return a and b parameters in the form E* = aE - b that minimizes the
% difference between corresponding E matricies. If the matrix is larger, it
% recursively calls itself for each half of the matrix, obtains a and b
% parameters, obtains E1* and E2* from the a, b and E values and returns
% new parameters. This algorithm should not modify E1, E2, x and y global
% variables, but instead call the wrapper to do so.

% Unlike the optimization wrapper, this function returns an M x N matrix of
% parameters for a and b. This function does not handle binning of edge
% matricies, that is done in the main function.

function [ a, b ] = a_b_recursive( E )

%% Obtain the size parameters from input matrix

[M, N, ~, y, ~] = size(E); %Determine tile layout

xlength = round((M + 0.5) / 2); %Number of tiles if split in x direction
ylength = round((N + 0.5) / 2); %Number of tiles if split in y direction

if xlength > 2
    xlength = xlength - mod(xlength, 2); 
end

if ylength > 2
    ylength = ylength - mod(ylength, 2);    
end

%% Determine which algorithm to run depending on M x N dimensions of input matrix

if M == 2 && N == 2
    
    qq = E(1,1,:,:,2);
    ww = E(2,1,:,:,3);
    ee = E(2,2,:,:,4);
    rr = E(1,2,:,:,1);
    tt = E(2,1,:,:,4);
    yy = E(2,2,:,:,1);
    uu = E(1,2,:,:,2);
    ii = E(1,1,:,:,3);
    
    [ a, b ] = optimize_a_b_2x2(qq,ww,ee,rr,tt,yy,uu,ii);
    
elseif M == 2 && N == 1 % Run optimization for 2x1
    
    [ a, b ] = optimize_a_b(E(1,1,:,:,2), E(2,1,:,:,4));
    a = a';
    b = b';
    
elseif M == 1 && N == 2 % Run optimization for 1x2
    
    [ a, b ] = optimize_a_b(E(1,1,:,:,3), E(1,2,:,:,1));
    
elseif M == 1 && N == 1 % Optimization is finished for this tile
    
    a = 1;
    b = 0;
    
else
    
    [temp_a11, temp_b11] = a_b_recursive(E(1:xlength,1:ylength,:,:,:)); % Obtain parameters for left side tiles
    [temp_a21, temp_b21] = a_b_recursive(E((xlength+1):M,1:ylength,:,:,:)); % Obtain parameters for right side tiles
    [temp_a12, temp_b12] = a_b_recursive(E(1:xlength,(ylength+1):N,:,:,:)); % Obtain parameters for left side tiles
    [temp_a22, temp_b22] = a_b_recursive(E((xlength+1):M,(ylength+1):N,:,:,:)); % Obtain parameters for right side tiles
    
    % Overlap 1 (vertical top)  
    for j = 1:ylength
        left_temp(:,:) = E(xlength,j,:,:,2) * temp_a11(xlength,j) - temp_b11(xlength,j);
        right_temp(:,:) = E(xlength+1,j,:,:,4) * temp_a21(1,j) - temp_b21(1,j);
        qq (:,((j-1)*y+1):(j*y)) = left_temp;
        tt (:,((j-1)*y+1):(j*y)) = right_temp; 
    end
    
    % Overlap 2 (horizontal right)
    for i = 1:(M-xlength)
        top_temp(:,:) = E(i+xlength,ylength,:,:,3) * temp_a21(i,ylength) - temp_b21(i,ylength);
        bottom_temp(:,:) = E(i+xlength,ylength+1,:,:,1) * temp_a22(i,1) - temp_b22(i,1);
        ww (:,((i-1)*y+1):(i*y)) = top_temp;
        yy (:,((i-1)*y+1):(i*y)) = bottom_temp; 
    end

    % Overlap 3 (vertical bottom)
    for j = 1:(N-ylength)
        left_temp(:,:) = E(xlength,j+ylength,:,:,2) * temp_a12(xlength,j) - temp_b12(xlength,j);
        right_temp(:,:) = E(xlength+1,j+ylength,:,:,4) * temp_a22(1,j) - temp_b22(1,j);
        uu (:,((j-1)*y+1):(j*y)) = left_temp;
        ee (:,((j-1)*y+1):(j*y)) = right_temp; 
    end
    
    % Overlap 4 (horizontal left)
    for i = 1:xlength
        top_temp(:,:) = E(i,ylength,:,:,3) * temp_a11(i,ylength) - temp_b11(i,ylength);
        bottom_temp(:,:) = E(i,ylength+1,:,:,1) * temp_a12(i,1) - temp_b12(i,1);
        rr (:,((i-1)*y+1):(i*y)) = bottom_temp;
        ii (:,((i-1)*y+1):(i*y)) = top_temp; 
    end
    
    [super_a, super_b] = optimize_a_b_2x2(qq,ww,ee,rr,tt,yy,uu,ii);

    a(1:xlength,1:ylength) = temp_a11 * super_a(1,1);
    a((xlength+1):M,1:ylength) = temp_a21 * super_a(2,1);
    a(1:xlength,(ylength+1):N) = temp_a12 * super_a(1,2);
    a((xlength+1):M,(ylength+1):N) = temp_a22 * super_a(2,2);  
    
    b(1:xlength,1:ylength) = temp_b11 * super_a(1,1) + super_b(1,1);
    b((xlength+1):M,1:ylength) = temp_b21 * super_a(2,1) + super_b(2,1);
    b(1:xlength,(ylength+1):N) = temp_b12 * super_a(1,2) + super_b(1,2);
    b((xlength+1):M,(ylength+1):N) = temp_b22 * super_a(2,2) + super_b(2,2);    

end

