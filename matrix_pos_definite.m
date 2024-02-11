function [returnString] = matrix_pos_definite( A )

%Purpose: check matrix for positive definiteness using eigenvalues and rank

% Inputs
% A - symmetric matrix

eig_A = eig(A);
flag = 0;
for i=1:rank(A)
   if(eig_A(i) <= 0)
       flag = 1;
   end
end

if(flag == 1)
	returnString='matrix is NOT positive definite';
	else
	returnString='matrix is positive definite';
end

end