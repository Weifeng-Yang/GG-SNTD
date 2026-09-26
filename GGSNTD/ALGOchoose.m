% Input.
% var,core    : initial matrix and core tensor
% ngmar       : decomposed tensor
% The remaining parameters are explained the same as the 'main_Run_me' function

% Output.
% cores, vars : Decomposition matrix and core tensor resulting from the final iterative result
% loss:       : Array of loss functions generated during iteration
% tr:         : Runtime array during iteration

function [data,varss]=ALGOchoose(core,var,ngmar,maxiteropt,Rdims,flag,stopindex,r,alphat,lamda,kiter)

if(flag==1)
[vars,loss,tr]=L1GSNMF(var,ngmar,maxiteropt,stopindex,r,lamda(1));
varss{1}=vars;
lossdata=loss;
trdata=tr;
 
  

elseif(flag==2)
[cores,vars,loss,tr]=SNTD(core,var,ngmar,maxiteropt,stopindex,r,lamda(1));
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  





elseif(flag==3)
[cores,vars,loss,tr]=GSNTD(core,var,ngmar,maxiteropt,Rdims,stopindex,r,1.7,10);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  




elseif(flag==4) 
[cores,vars,loss,tr]=MGNTD(core,var,ngmar,maxiteropt,stopindex,r,kiter);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  



elseif(flag==5) 
[cores,vars,loss,tr]=AMGRNTD(core,var,ngmar,maxiteropt,stopindex,r);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  



elseif(flag==6) 
[cores,vars,loss,tr]=L1SGNTD(core,var,ngmar,maxiteropt,stopindex,r,alphat,lamda,0.9999);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  



elseif(flag==7) 
[cores,vars,loss,tr]=GGSNTD(core,var,ngmar,maxiteropt,stopindex,r,0,lamda,0.9999);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  



elseif(flag==8) 
[cores,vars,loss,tr]=GGSNTD(core,var,ngmar,maxiteropt,stopindex,r,alphat,lamda,0.9999);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  
end


data{1}=lossdata;
data{2}=trdata;







end
