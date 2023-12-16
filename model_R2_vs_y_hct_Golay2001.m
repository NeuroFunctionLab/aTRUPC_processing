function result = model_R2_vs_y_hct_Golay2001 (beta,x)
%Input:
%	beta   --  parameter of the function, must have 2 element
%	x      --  input vector, must have one column
%Output:
%	result --  corresponding function value, 
%             has the same dimension as x
%Hanzhang Lu
%Date: 12/08/2009

ii=size(x,1);
result=zeros(ii,1);
A=beta(1);
B=beta(2);
C=beta(3);
D=beta(4);
E=beta(5);
F=beta(6);
for i=1:ii
result(i) = (A*(1-x(i,2))+B*x(i,2)+C*x(i,2)*(1-x(i,2)))+(D*x(i,2)+E*x(i,2)^2)*(1-x(i,1))+F*x(i,2)*(1-x(i,2))*(1-x(i,1))^2;
end
