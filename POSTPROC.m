function[axstrain,totstress]=POSTPROC(us,ncount,ecount,elcos,elsin,l,ae);
[numnod,x,y]=textread('C:\Users\Matt\Downloads\Model (extract.me)\Model\node.txt','%u %20.13f %20.13f');
[numelem,ni,nj,area, m]=textread('C:\Users\Matt\Downloads\Model (extract.me)\Model\elem.txt','%u %u %u %20.13f %f');
e=30e6;
% calculate every element local stiffness matrix
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
    ug=[us(3*na-2); us(3*na-1); us(3*na); us(3*nb-2); us(3*nb-1); us(3*nb);];
    % calculate element forces using the equation {f}=[k]*{u}
    % in the global cartesian 
    fg=kloc*ug;
    % tmat as the transformation matrix
    % transforming the displacements of the nodes of each element 
    % to the local element coordinate system, i.e. axial direction as x. 
    uloc=r*ug;
    fprintf('%d \n', uloc);
    fprintf('\n');
    
    syms Ro Ri
    A = pi()*(Ro^2-Ri^2);
    I = pi()*(Ro^4-Ri^4)/4;
    [Ro, Ri] = solve(A==ae(i), I==m(i),Ro,Ri);
    Ro2 = max(double(Ro));
    fprintf('Ro2 %d \n', Ro2);
    % calculating element axial strains and stresses using the simple
    % equation strain=(deltaU)/(element length) 
    % and      stress=Young's modulus*axial strains
    uxa=uloc(1);
    uxb=uloc(4);
    uya=uloc(2);
    uyb=uloc(5);
    uma = uloc(3);
    umb = uloc(6);
    axstrain(i)=(uxb-uxa)/l(i);
    axstress(i)=axstrain(i)*e;
    bxstress1(i)= Ro2*e*(6*(uyb-uya)/(l(i)^2)-2*(2*uma+umb)/l(i));
    bxstress2(i)= Ro2*e*(6*(uya-uyb)/(l(i)^2)+2*(2*umb+uma)/l(i));
    fprintf('bx1 %d \n', bxstress1(i));
    fprintf('bx2 %d \n', bxstress2(i));
    fprintf('ax %d \n', axstress(i));
    maxstress = max(abs(bxstress1(i)), abs(bxstress2(i)));
    if axstress(i) < 0
       totstress(i) = axstress(i)-maxstress;
    else
       totstress(i) = axstress(i)+maxstress;
    end
   end
   