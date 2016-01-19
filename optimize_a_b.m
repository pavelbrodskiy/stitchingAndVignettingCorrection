% The purpose of this function is to serve as a wrapper for the
% optimization problem of two overlapping regions (E1 and E2). The function
% takes in two regions and returns the a and b values which will result in
% the most seamless fit of the two images. No other function should modify 
% E1 and E2. This function returns a and b as arrays of length 2.

function [ a, b ] = optimize_a_b(matrix1, matrix2)

    global E1 E2 output;
    global lb ub x0 MFE MI TolCon TolFun TolX time ST1 Trial
    
    [row, ~] = size(output);
    
    E1 = matrix1;
    E2 = matrix2;

    problem.objective = @optimization_overlaps;
    problem.nvars = 4;
    problem.lb = [ones(1,2)*lb(1) ones(1,2)*lb(2)];
    problem.ub = [ones(1,2)*ub(1) ones(1,2)*ub(2)];
    problem.solver = 'fmincon';
    problem.x0 = [ones(1,2)*x0(1) ones(1,2)*x0(2)];
    problem.options = optimoptions('fmincon','Display','none');
    problem.options.MaxFunEvals = MFE;
    problem.options.MaxIter = MI;
    problem.TolCon = TolCon;
    problem.TolFun = [ones(1,2)*TolFun(1) ones(1,2)*TolFun(2)];
    problem.TolX = TolX;
    
    gs = GlobalSearch('MaxTime',time,'NumStageOnePoints',ST1,'NumTrialPoints',Trial,'Display','off');
    
    [solution,fval,flag,~,solutions] = run(gs, problem);
    
    output{row+1,1} = fval;
    output{row+1,2} = solutions;
    
    disp(['2x1 GS finished with flag: ' num2str(flag) ' GS run: ' num2str(row+1)]);
    
    a = reshape(solution(1:2),[1, 2]);
    b = reshape(solution(3:4),[1, 2]);
    
    %if solution(1) > 1
    %    a = [solution(1), 1];
    %else
    %    a = [1 1/solution(1)];
    %end
    
    %if solution(2) > 0
    %    b = [solution(2), 0];
    %else
    %    b = [0, -solution(2)];        
    %end
    
end

