function [ke,fe] = domain(xi, eta, coord, k, Q)

% This function calculate element Stiffness Matrix

% INPUT:
% ======
% k = Thermal_conductivity
% Q = Heat generation per unit volume
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% ke = Element stiffness matrix
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);


J = [x(1)-x(3), y(1)-y(3); x(2)-x(3), y(2)-y(3)];
N = [xi, eta, 1-xi-eta];
B = inv(J)*[1, 0, -1; 0, 1, -1];


ke = 0.5*k*det(J)*B'*B;

fe = N'*0.5*Q*det(J);
