function [ Ry,Rz1,Ry2 ] = RotMatrix( phi,betaZ,psi )
% Objective is to express the three Euler's angles by three rotation
% matricees in the coordinates of initial tube axis.
%
%   [ Ry,Rz1,Ry2 ] = RotMatrix( phi,betaZ,psi )
%
% VERSION-2: possibility to start with a coordinates basis already rotated
%
% Inputs :
%   - [ex,ey,ez]: A coordinates axis base
%   - phi,betaZ, psi : Euler's rotation angles, corresponding resp. Ry,Rz1,Ry2
% Outputs:
%   - Three rotation matrices
%
% Keerthi Krishna PARVATHANENI  2014

ex=[1 0 0];  ey=[0 1 0];  ez=[0 0 1];

A = [0 ey(3) -1*ey(2); -1*ey(3) 0 ey(1); ey(2) -1*ey(1) 0];
B = ey'*ey;
C = diag([1,1,1]) - B;
Ry = A.*sin(phi) + C.*cos(phi) + B;


ez1=ez*Ry;
  
A = [0 ez1(3) -1*ez1(2); -1*ez1(3) 0 ez1(1); ez1(2) -1*ez1(1) 0];
B = ez1'*ez1;
C = diag([1,1,1]) - B;
Rz1 = A.*sin(betaZ) + C.*cos(betaZ) + B;

ey2 = ey*Rz1;

  
A = [0 ey2(3) -1*ey2(2); -1*ey2(3) 0 ey2(1); ey2(2) -1*ey2(1) 0];
B = ey2'*ey2;
C = diag([1,1,1]) - B;
Ry2 = A.*sin(psi) + C.*cos(psi) + B;

end

