% The purpose of this function is to serve as a wrapper for the
% optimization problem of two overlapping regions (E1 and E2). The function
% takes in two regions and returns the a and b values which will result in
% the most seamless fit of the two images. No other function should modify 
% E1 and E2. This function returns a and b as arrays of length 2.

function [ a, b ] = optimize_a_b_2x2(q, w, e, r, t, y, u, i)

    global E11 E12 E13 E14 E21 E22 E23 E24 output;
    global lb ub x0 MFE MI TolCon TolFun TolX time ST1 Trial
    
    E11 = q;
    E12 = w;
    E13 = e;
    E14 = r;
    E21 = t;
    E22 = y;
    E23 = u;
    E24 = i;
    
    [row, ~] = size(output);

    problem.objective = @optimization_overlaps_2x2;
    problem.nvars = 8;
    problem.lb = [ones(1,4)*lb(1) ones(1,4)*lb(2)];
    problem.ub = [ones(1,4)*ub(1) ones(1,4)*ub(2)];
    problem.solver = 'fmincon';
    problem.x0 = [ones(1,4)*x0(1) ones(1,4)*x0(2)];
    problem.options = optimoptions('fmincon','Display','none');
    problem.options.MaxFunEvals = MFE;
    problem.options.MaxIter = MI;
    problem.TolCon = TolCon;
    problem.TolFun = [ones(1,4)*TolFun(1) ones(1,4)*TolFun(2)];
    problem.TolX = TolX;
    
    gs = GlobalSearch('MaxTime',time,'NumStageOnePoints',ST1,'NumTrialPoints',Trial,'Display','off');

    [solution,fval,flag,~,solutions] = run(gs, problem);
    
    output{row+1,1} = fval;
    output{row+1,2} = solutions;
    
    disp(['2x2 GS finished with flag: ' num2str(flag) ' GS run: ' num2str(row+1)]);
    
    a = reshape(solution(1:4),[2 2]);
    b = reshape(solution(5:8),[2 2]);
    
end

