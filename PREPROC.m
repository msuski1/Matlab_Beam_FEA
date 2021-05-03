function[elcos,elsin,l,ncount,ecount,exx,ar,m]=PREPROC;
% We request for data be printed on the screen in long format
format long;
% We read columns in the file 'node.txt' which is located at 'c:\model' 
% in the
% following variables
% 1- The first column of the file is read into the 'numnode' array, this 
% integer array represents
% the number assigned to each node
% 2- The second column of the file is read into the 'x' array, this
% real array
% represents the corresponding node x coordinate 
% 3- The third column of the file is read into the 'y' array, this 
% real array
% represents the corresponding node y coordinate
[numnod,x,y]=textread('C:\Users\Matt\Downloads\Model (extract.me)\Model\node.txt','%u %20.13f %20.13f');
% Here we read columns in the file 'elem.txt' which is located at
% 'c:\model' in the
% following variables
% 1- The first column of the file is read into the 'numelem' array, this
% integer array represents
% the number assigned to each element
% 2- The second column of the file is read into the 'ni' array, this 
% integer array
% represents the corresponding node i of each element.
% 3- The third column of the file is read into the 'nj' array, this array
% represents the corresponding node j of each element.
% 4- The fourth column of the file is read into 'ar' array, which
% represents cross-sectional area of each element.
[numelem,ni,nj,ar,m]=textread('C:\Users\Matt\Downloads\Model (extract.me)\Model\elem.txt','%u %u %u %20.13f %f');
% we have to know the number of nodes and elements defined
% so we make use of the 'length' command. this command returns the size of
% the array. So the number of defined nodes is equal to the length of array
% numnod and the number of defined elements is equal to the length of array
% numelem.
ncount=length(numnod);
ecount=length(numelem);
% Assign the material and geometric properties of elements. A constant
% value is assumed here.
exx=30e6;
% Here the program calculates the length and the element angle 
% (cosine and sine of the angle) with respect 
% to the global X direction (see figure 1)
% Plot what has been inputted. This is the data check.
for i=1:ecount;
    a=ni(i);
    b=nj(i);
    l(i)=sqrt((x(b)-x(a))^2+(y(b)-y(a))^2);
    elcos(i)=(x(b)-x(a))/l(i);
    elsin(i)=(y(b)-y(a))/l(i);
end
% Here the function plots the FE model
i=1;
a=[x(ni(i)) x(nj(i))];
b=[y(ni(i)) y(nj(i))];
plot(a,b,'k+-');
xlim([-5,50]);
ylim([-5, 20]);

hold;
for i=2:ecount;
    a=[x(ni(i)) x(nj(i))];
    b=[y(ni(i)) y(nj(i))];
    plot(a,b,'k+-');

end
hold;



