% Requires:     pauliprod2.m 
% Author:       Xiao-Ming Lu (luxiaoming@gmail.com)
% Date:         2010/4/17

%
% Description: Generate a matrix whose elements are the expected value of pauli tensor products on a 4X4 density matrix.
%
% Usage: cm = correlation_matrix (rho)
%     rho: The density matrix (Require 4by4).

function cm = dm2cm (rho)
global pauli2;
if(isempty(pauli2))
    pauliprod2;
end%if

cm = zeros(4,4);

for k1=1:4
    for k2=1:4
        cm(k1,k2) = real(trace(rho * pauli2{k1,k2}));
    end%for
end%for