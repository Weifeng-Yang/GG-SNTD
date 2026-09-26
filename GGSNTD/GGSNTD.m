%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [core,var,loss,timerun]=GGSNTD(core,var,ngmar,maxiteropt,stopindex,r,alphat,lamda,btmax)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));
LK=zeros(1,num);
LCK=0;
L=ones(1,num);
LC=1;
tk=1;
rho=1e-10;

Lapk=LLaplace(ngmar);
varK=var;
coreK=core;
wk=(tk-1)/(tk);

returnloss=norm(tensor(ngmar));
for i=1:num
    alpha(i)=0;
end
alpha(num)=alphat;
loss(1)=computeloss(ngmar,core,var,lamda,alpha,Lapk);


rate=[0];
t1=clock;


for i=1:maxiteropt
%% update parameters

fprintf("%d\n",i);

LtempK=LK;
LCtempK=LCK;
vartempK=varK;
coretempK=coreK;    
varK=var;
coreK=core;
LK=L;
LCK=LC;    
[core,var,LC,L,btc,bt]=GGSNTDupdate(core,var,num,ngmar,r,coretempK,vartempK,wk,L,LC,LtempK,LCtempK,lamda,alpha,Lapk,btmax);
loss(i+1)=computeloss(ngmar,core,var,lamda,alpha,Lapk);

check=norm(tensor(core-coreK-btc*(coreK-coretempK)))^2;
for j=1:num
    check=check+norm(var{j}-varK{j}-bt(j)*(varK{j}-vartempK{j}),'fro')^2;
end
%% Judging whether to extrapolate
if(loss(i+1)>loss(i)-rho*check)
    var=varK;
    core=coreK;
    L=LK;
    LC=LCK;
    [core,var,LC,L,btc,bt]=GGSNTDupdate(core,var,num,ngmar,r,coretempK,vartempK,0,L,LC,LtempK,LCtempK,lamda,alpha,Lapk,btmax);
    loss(i+1)=computeloss(ngmar,core,var,lamda,alpha,Lapk);
end


%% Check if termination condition is met

fprintf("GGSNTD\n");

fprintf("nonzero:%d\n",nnz(core));  

check1=norm(tensor(core-coreK))^2;
check2=norm(tensor(coreK))^2;
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));  
    check1=check1+norm(var{j}-varK{j},'fro')^2;
    check2=check2+norm(varK{j},'fro')^2;
end

bts{i}=[btc,bt];
t2=clock;
timerun(i+1)=etime(t2,t1);
% Res=abs(loss(i+1)-loss(i));
Res=sqrt(check1/check2);
fprintf("cri：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end

tknew=(1+sqrt(1+4*tk^2))/2;
wk=(tk-1)/tknew;
tk=tknew;



end
end


function loss=computeloss(ngmar,core,var,lamda,alpha,Lapk)
    loss=compute(core,var,ngmar);
    coret=double(tenmat(core,length(size(ngmar))));
    temp=var{end}*coret;
    loss=loss+alpha(end)/2*trace(temp'*Lapk{end}*temp);
    
    for j=1:length(var)-1
        loss=loss+lamda(j)*sum(sqrt(sum(var{j}.^2, 2)));
    end

end
