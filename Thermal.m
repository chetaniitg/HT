k = 1.5;  h = 50; Tinf = 25; Q = 0;

% Gauss Points and weights for one point gauss quadrature
% for triangular area integration
% -------------------------------------------------------
xi1 = 1/3;  eta1 = 1/3;  w1 = 1;

% Gausspoints for line integration
% --------------------------------
beta1 = -0.774597;   ew1 = 5/9;
beta2 = 0.774597;    ew2 = 5/9;
beta3 = 0;           ew3 = 8/9;

% Element 1
coord = [0.0, 0.0;
         0.4, 0.0;
         0.4, 0.15];
[ke1d,fe1d] = domain(xi1, eta1, coord, k, Q);

[ke1hg1,fe1hg1] = gamah(beta1, coord, h, Tinf, 2);
[ke1hg2,fe1hg2] = gamah(beta2, coord, h, Tinf, 2);
[ke1hg3,fe1hg3] = gamah(beta3, coord, h, Tinf, 2);
ke1h = ke1hg1*ew1 + ke1hg2*ew2 + ke1hg3*ew3;
fe1h = fe1hg1*ew1 + fe1hg2*ew2 + fe1hg3*ew3;

ke1 = ke1d + ke1h;
fe1 = fe1d + fe1h;


% Element 2
coord = [0.0, 0.0;
         0.4, 0.15;
         0.0, 0.3];
[ke2,fe2] = domain(xi1, eta1, coord, k, Q);

% Element 3
coord = [0.4, 0.15;
         0.4, 0.3;
         0.0, 0.3];
[ke3d,fe3d] = domain(xi1, eta1, coord, k, Q);

[ke3hg1,fe3hg1] = gamah(beta1, coord, h, Tinf, 1);
[ke3hg2,fe3hg2] = gamah(beta2, coord, h, Tinf, 1);
[ke3hg3,fe3hg3] = gamah(beta3, coord, h, Tinf, 1);
ke3h = ke3hg1*ew1 + ke3hg2*ew2 + ke3hg3*ew3;
fe3h = fe3hg1*ew1 + fe3hg2*ew2 + fe3hg3*ew3;

ke3 = ke3d + ke3h;
fe3 = fe3d + fe3h;


% Assembly
K = zeros(5,5);
F = zeros(5,1);

K(1:3,1:3) = ke1(1:3,1:3);
K([1,3,5],[1,3,5]) = K([1,3,5],[1,3,5]) + ke2(1:3,1:3);
K([3,4,5],[3,4,5]) = K([3,4,5],[3,4,5]) + ke3(1:3,1:3);

F(1:3) = fe1(1:3);
F([1,3,5]) = F([1,3,5]) + fe2(1:3);
F([3,4,5]) = F([3,4,5]) + fe3(1:3);


% Imposition of B.C.
Kreduce = K(1:3,1:3);
Freduce = F(1:3) - (K(4,1:3)*180 + K(5,1:3)*180)';

% Finding Solution
ureduce = inv(Kreduce)*Freduce;
un = [ureduce;180;180]

% % ==============================================
% %% Printing Intermediate Result to The Output File
% % ------------------------------------------------
fid=fopen('Steps','w');
fprintf(fid,'The Element Stiffness matrices are\n');
fprintf(fid,'===================================\n');
fprintf(fid,'k = %12.4e, h = %12.4e, Tinf = %12.4e, Q = %12.4e\n\n',k,h,Tinf,Q);
fprintf(fid,'Ke1 \n');
fprintf(fid,'----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',ke1(i,1:3));
end
fprintf(fid,'Ke2 \n');
fprintf(fid,'----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',ke2(i,1:3));
end
fprintf(fid,'Ke3 \n');
fprintf(fid,'----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',ke3(i,1:3));
end



fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');

fprintf(fid,'fe1 \n');
fprintf(fid,'-----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',fe1(i));
end
fprintf(fid,'fe2 \n');
fprintf(fid,'-----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',fe2(i));
end
fprintf(fid,'fe3 \n');
fprintf(fid,'-----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',fe3(i));
end

fprintf(fid,'\n\nThe Global Stiffness matrix is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K\n');
fprintf(fid,'--\n');
for i = 1:5
   fprintf(fid,'%12.4e\t8%12.4e\t%12.4e\t%12.4e\t%12.4e\n\n',K(i,1:5));
end

fprintf(fid,'\n\nThe Global Load Vector is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'F\n');
fprintf(fid,'--\n');
for i = 1:5
   fprintf(fid,'%12.4e\n\n',F(i));
end

fprintf(fid,'\n\nImposition of Boundary Condition\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K u = F\n');
fprintf(fid,'--------\n');
for i = 1:5
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',K(i,1:5),F(i));
end

fprintf(fid,'\n\nReduced Equations\n');
fprintf(fid,'=======================\n\n');
fprintf(fid,'K_reduce u = F_reduce\n');
fprintf(fid,'--------\n');
for i = 1:3
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',Kreduce(i,1:3),Freduce(i));
end

fprintf(fid,'\n\nThe Final Solution\n');
fprintf(fid,'=========================\n\n');
fprintf(fid,'Tn\n');
fprintf(fid,'--\n');
for i = 1:5
   fprintf(fid,'%12.4e\n\n',un(i));
end

