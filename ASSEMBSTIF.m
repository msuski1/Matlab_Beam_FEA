function[gstif]=ASSEMBSTIF(elcos,elsin,l,ncount,ecount,e,ae,m);
[numelem,ni,nj,test,moment]=textread('C:\Users\Matt\Downloads\Model (extract.me)\Model\elem.txt','%u %u %u %20.13f %f');
ndim=3*ncount;
for i=1:ndim;
    for j=1:ndim
        gstif(i,j)=0.0;
    end
end
for i=1:ecount;
    landa=elcos(i);
    mu=elsin(i);
    ael=e*ae(i)/l(i);
    kloc(1,1)=ael;
    kloc(1,2)=0;
    kloc(1,3)=0;
    kloc(1,4)=-ael;
    kloc(1,5)=0;
    kloc(1,6)=0;
    kloc(2,2)=12*e*m(i)/(l(i)^3);
    kloc(2,3)=6*e*m(i)/(l(i)^2);
    kloc(2,4)=0;
    kloc(2,5)=-12*e*m(i)/(l(i)^3);
    kloc(2,6)=6*e*m(i)/(l(i)^2);
    kloc(3,3)=4*e*m(i)/l(i);
    kloc(3,4)=0;
    kloc(3,5)=-6*e*m(i)/(l(i)^2);
    kloc(3,6)=2*e*m(i)/l(i);
    kloc(4,4)=ael;
    kloc(4,5)=0;
    kloc(4,6)=0;
    kloc(5,5)=12*e*m(i)/(l(i)^3);
    kloc(5,6)=-6*e*m(i)/(l(i)^2);
    kloc(6,6)=4*e*m(i)/l(i);
    for j=2:6;
        for k=1:(j-1);
            kloc(j,k)=kloc(k,j);
        end
    end
    r(1,:) = [landa, mu, 0,0,0,0];
    r(2,:) = [-mu,landa, 0,0,0,0];
    r(3,:) = [0,0,1,0,0,0];
    r(4,:) = [0,0,0,landa,mu,0];
    r(5,:) = [0,0,0,-mu,landa,0];
    r(6,:) = [0,0,0,0,0,1];
    rt = transpose(r);
    kloc = rt*kloc*r;
    % here to save space we assemble the local stiffness matrix [kloc] to 
    % the global stiffness matrix [gstif]
    na=ni(i);
    nb=nj(i);
    npe=3;
    id=[npe*na-2 npe*na-1 npe*na npe*nb-2 npe*nb-1 npe*nb];
    for j=1:6;
        for k=1:6;
            a=id(j);
            b=id(k);
            val=kloc(j,k);
            gstif(a,b)=gstif(a,b)+val;
        end
    end
    kloc=0.0;
end

    