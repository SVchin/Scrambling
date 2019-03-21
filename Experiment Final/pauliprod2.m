% Requires:		nothing
% Author:       Xiao-Ming Lu (luxiaoming@gmail.com)
% Date:         2010/4/17
%
% Description: Generate the tensor products of Pauli matrices.

global pauli pauli2;
pauli = cell(4,1);
pauli2 = cell(4);
pauli{1} = eye(2);
pauli{2} = [0 1; 1 0];
pauli{3} = [0 -i; i 0];
pauli{4} = [1 0; 0 -1];

for k1=1:4
    for k2=1:4
        pauli2{k1,k2} = kron(pauli{k1},pauli{k2});
    end%for
end%for