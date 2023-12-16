function result = multi_monexp_model_nte (beta,x)
%Input:
%	beta   --  parameter of the function, must have >2 element
%	x      --  input vector, must be a column vector, the last element
%	should be the number of TEs (nTE)
%Output:
%	result --  corresponding function value, length(result) == length(x)-1
%Peiying Liu, 5/29/2012
% Dengrong Jiang, August 2022, remove the use of global variable
% Example: [beta_VOI, resid_VOI, jacob_VOI] = nlinfit_hlu([eTE_All;length(eTE(ind2fit))], sig2fit, 'multi_monexp_model_nte', init_beta);

s0 = beta(1:end-1);
r2_star=beta(end);

nTE = x(end);
x1 = x(1:end-1);
result = zeros(length(x1), 1);

for i=1:length(s0)
    result((i-1)*nTE+1:i*nTE) = s0(i)*exp(-x1((i-1)*nTE+1:i*nTE).*r2_star);
end
