function [x,lambda]=EqualityQPSolver(H,g,A,b,solver)
switch solver
    case 'LUdense'
        [x,lambda]=EqualityQPSolverLUdense(H,g,A,b);
    case 'LUsparse'
        [x,lambda]=EqualityQPSolverLUsparse(H,g,A,b);
    case 'LDLdense'
        [x,lambda]=EqualityQPSolverLDLdense(H,g,A,b);
    case 'LDLsparse'
        [x,lambda]=EqualityQPSolverLDLsparse(H,g,A,b);
    case 'RangeSpace'
        [x,lambda]=EqualityQPSolverRangeSpace(H,g,A,b);
    case 'NullSpace'
        [x,lambda]=EqualityQPSolverNullSpace(H,g,A,b);
end
end