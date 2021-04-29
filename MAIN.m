%main program
clear all;
% first we call function PREPROC to do preprocessing step
% this function reads the input file for node coordinates and
% for element attributes (node numbers)
% then based on the input data it calculates every element length 
% and element direction (Cosine and Sine)
% the function returns the following variables to the main program
% 1- [elc] which is the array of element cosine, its dimension is the
%    same as
% the number of elements
% 2- [els] which is the array of element sin, its dimension is the same as
%    the number of elements
% 3- [elen] which is the array of element length, its dimension is the 
%    same as
% the number of elements
% [nc] is a single variable for number of nodes defined
% [ec] is a single variable for number of elements defined
[elc,els,elen,nc,ec,ex,a,m]=PREPROC;
% Here we assume the Young's modulus of 200E9 (Pa)
% and element cross section of 0.0025 (m^2)
% ***********************************************************************
% now we call function assembstif to calculate every element local
% stiffness matrix and finally construct the global stiffness matrix
ex=30e6;
[kglob]=ASSEMBSTIF(elc,els,elen,nc,ec,ex,a,m);
% now it is time to call function BC to apply boundary conditions 
% and also generate the force vector ([K][U]=[F])
Original_kglob=kglob;
[kglob,fglob]=BC(kglob);
% to calculate displacement vector we make use of {U}=[K]^-1*{F}
u=inv(kglob)*transpose(fglob);
% now time to calculate element forces and stresses
% this procedure will be done thru the function POSTPROC
[axstr,totsts]=POSTPROC(u,nc,ec,elc,els,elen,a);
% Calculating force vector {f}=[K]{U}
Original_fglob=Original_kglob*u;
fid = fopen('C:\Users\Matt\Downloads\Model (extract.me)\Model\truss.dat','w+');
fprintf(fid,'Nodal displacements in global cartesian \n');
fprintf(fid,'Node number         UX                    UY                    UM      \n');
for i=1:nc;
    fprintf(fid,'%8u %20.13f %20.13f %20.13f\n',[i;u(3*i-2);u(3*i-1);u(3*i)]);
end;
fprintf(fid,'********************************************************** \n');
fprintf(fid,'Element calculation \n');
fprintf(fid,'Element number      Total stress     \n');
for i=1:ec;
    fprintf(fid,'%8u    %20.13f\n',[i;totsts(i)]);
end;
fprintf(fid,'********************************************************** \n');
fprintf(fid,'Force Vector \n');
fprintf(fid,'Node number         FX                    FY                    M     \n');
for i=1:nc;
    fprintf(fid,'%8u %20.13f %20.13f %20.13f\n',[i;Original_fglob(3*i-2);Original_fglob(3*i-1);Original_fglob(3*i)]);
end;
fclose(fid);



        