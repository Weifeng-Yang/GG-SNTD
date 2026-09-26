% Parameter update functionI
function [core,var,LC,L,btc,bt]=L1SGNTDupdate(core,var,num,ngmar,r,coreK,varK,wk,L,LC,LK,LCK,lamda,alpha,Lapk,btmax)
    btc=min(wk,btmax);
    core=core+btc*(core-coreK);
    [V,LC]=gradSGNTDcore(core,var,ngmar,r,num,Lapk,alpha(num));
    core=PROXn1(V);

    for j=1:num
    bt(j)=min(wk,btmax);
    if(j<num)
    var{j}=var{j}+bt(j)*(var{j}-varK{j});
    [V,L(j)]=gradSGNTD(core,var,ngmar,r,j,num,Lapk,alpha(j));
    var{j}= PROXL1(V,lamda(j),1/(L(j)*r));
    else
    var{j}=var{j}+bt(j)*(var{j}-varK{j});
    [V,L(j)]=gradSGNTD(core,var,ngmar,r,j,num,Lapk,alpha(j));
    var{j}=PROXn1(V);
    end
    end
end


function [V,L]=gradSGNTD(core,var,ngmar,r,n,num,Lapk,alpha)
core=tensor(core);
index=1:num;
index(n)=[];
coreg=ttm(core,var,index);
tempB=double(tenmat(coreg,n));
temp=tempB*tempB';
Xn=double(tenmat(ngmar,n));
if(n==num)
corenum=double(tenmat(core,num));
coretemp=corenum*corenum';
U=var{n}*temp-Xn*tempB'+alpha*Lapk{n}*var{n}*coretemp;
L=norm(temp,'fro')+alpha*norm(Lapk{n},'fro')*norm(coretemp,'fro'); %%
else
U=var{n}*temp-Xn*tempB'+alpha*Lapk{n}*var{n};
L=norm(temp,'fro')+alpha*norm(Lapk{n},'fro'); %%
end
V=var{n}-1/(r*L)*U;
end


function [V,LC,U]=gradSGNTDcore(core,var,ngmar,r,num,Lapk,alpha)
ngmar=tensor(ngmar);
core=tensor(core);
Vc=core;
Vx=ngmar;
temp=1;

for i=1:num
    vart=var{i}'*var{i};
    Vc=ttm(Vc,vart,i);
    Vx=ttm(Vx,var{i}',i);
    temp=temp*norm(vart,'fro');
end
U=Vc-Vx+alpha*ttm(core,var{num}'*Lapk{num}*var{num},num);
LC=temp+alpha*norm(var{num}'*Lapk{num}*var{num},'fro');
V=core-1/(r*LC)*U;

end

function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end
