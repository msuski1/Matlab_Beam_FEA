function[kg,fg]=BC(kg);
% now we have to read the load step file to impose BC and loads
% in the load step file we have the format like below
% F 1 FX 80000
% The label [F] means that this is a concentrated force 
% the number [1] means it is applied on node 1
% the label [FX] means it is applied along the X direction
% the number [80000] is the force magnitude
% similar way is used for displacement loading in the loas step file
% i.e. D 2 UY 0.0
[lid,lnum,ldir,lval]=textread('C:\Users\Matt\Downloads\Model (extract.me)\Model\load.txt','%1c %u %2c %20.13f');
nol=length(lid);
j=1;
k=1;
for i=1:nol;
    if lid(i)=='D' 
        dpdir(j,:)=ldir(i,:);
        dpval(j)=lval(i);
        dpnum(j)=lnum(i);
        j=j+1;
    elseif lid(i)=='F'
        fdir(k,:)=ldir(i,:);
        fval(k)=lval(i);
        fnum(k)=lnum(i);
        k=k+1;
    end
end
% We must know how many loading we have
nu=length(dpval);
nf=length(fval);
% the methodology to apply boundary condition is to zero the row and 
% the column of the correponding node number in the global stiffness matrix
% and make the diagonal member 1.0. It means that if on node number i we 
% have set the ux=0.0 we have to do the following in the global stiffness:
% kg(2*i-1,j)=kg(j,2*i-1)=0.0 (for j=1,2*ncount); kg(2*i-1,2*i-1)=1.0
for i=1:nu
    sa=dpdir(i,:);
    dnum=dpnum(i);
    if sa=='UX'
        for j=1:length(kg)
            kg(3*dnum-2,j)=0.0;
            kg(j,3*dnum-2)=0.0;
        end
        kg(3*dnum-2,3*dnum-2)=1.0;
    elseif sa=='UY'
        for j=1:length(kg)
            kg(3*dnum-1,j)=0.0;
            kg(j,3*dnum-1)=0.0;
        end
        kg(3*dnum-1,3*dnum-1)=1.0;
     elseif sa=='UM'
        for j=1:length(kg)
            kg(3*dnum,j)=0.0;
            kg(j,3*dnum)=0.0;
        end
        kg(3*dnum,3*dnum)=1.0;
    end
end
for i=1:length(kg)
    fg(i)=0.0;
end
for i=1:nf;
    sa=fdir(i,:);
    fno=fnum(i);
    if sa=='FX'
        fg(3*fno-2)=fval(i);
    elseif sa=='FY'
        fg(3*fno-1)=fval(i);
    elseif sa=='FM'
        fg(3*fno)=fval(i);
    end
end
    