function [A,lambda] = convert_parameters(P,fitType)
% Converts parameters used for fitting into parameters reported in paper
if(strcmp(fitType,'eq5'))
    A = 10.^P(3);
    lambda = 10.^P(3)*exp(P(2));
    fprintf('A = %f\n',A);
    fprintf('l = %f\n',lambda);
    fprintf('tau1 = %f\n',P(1));
elseif(strcmp(fitType,'eq6'))
    A = 10.^P(4);
    lambda = 10.^P(4)*exp(P(3));
    fprintf('A = %f\n',A);
    fprintf('l = %f\n',lambda);
    fprintf('tau1 = %f\n',P(1));
    fprintf('tau2 = %f\n',P(2));
else
    error('fitType must be ''eq5'' or ''eq6''');
end