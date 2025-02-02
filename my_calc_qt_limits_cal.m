function [model] = my_calc_qt_limits_cal (Xcal,Ycal,Xnew,nvl);
% Routine that calculates the PLS-DA model and the values of Q and T2 and
% their respective limits for each sample. 
% INPUT
% Xcal = Calibration Samples
% Ycal = Classes Matrix
% Xnew = Validation Samples
% nvl = Number of Latent Variables

% OUTPUT
% model = Arquivo estrutura contendo o modelo PLS-DA

% Rotina desenvolvida por Felipe Bachion
% para maiores informa��es ou informar erros contate por e-mail
% felipebachion@gmail.com

[Xce,mX,Xnew_ce]=centrar(Xcal,1,Xnew);		
	[Yce,mY]=centrar(Ycal,1);
	[B,C,P,T,U,R,R2X,R2Y]=simpls(Xce,Yce,nvl,Xce'*Yce,[]);	        
	Yhat=Xce*B;			
	[Yhat]=centrar_inverso(Yhat,mY);
    
    
    Ynew_hat=Xnew_ce*B;			
	[Ynew_hat]=centrar_inverso(Ynew_hat,mY);
    
    
    Xrecuperado = T*P';
    [Xrecuperado]=centrar_inverso(Xrecuperado,mX);
    Qcont=Xcal-Xrecuperado;
    comp=nvl;
    nobj=size(Xcal,1);
    
  
   %%
   
   
   %    Calc T2
fvar = sqrt(1./(diag(T'*T)/(size(Xcal,1) - 1)));
Thot = sum((T*diag(fvar)).^2,2);
Tcont = (T*diag(fvar)*P');


    % T2 limit
lev_conf = 0.95;
if license('test','statistics_toolbox')
    F = finv(lev_conf,comp,nobj-comp);
    tlim = (comp*(nobj - 1)/(nobj - comp))*F;
else
    disp(['Enter F-Value with (Namostras-N� latent variables -1) degrees of freedom ','. ']);
    tlim=input('');   
end
   

%% Calc Q
    for i=1:size(T,1)
    Qres(i) = Qcont(i,:)*Qcont(i,:)';
    end   

% Q limit

[m,n] = size(Xcal);
s = svd(Qcont');
s = s.^2/(m-1);
Qcont=s;

[m,n] = size(Qcont);
if n>m
  s   = s';
  m   = n;
end
Qcont=s;
t1 = sum(Qcont(0+1:m,1));
t2 = sum(Qcont(0+1:m,1).^2);
t3 = sum(Qcont(0+1:m,1).^3);

ho = 1-2*t1*t3/3/(t2.^2);
if ho<0.001; ho = 0.001; end
ca = sqrt(2)*erfinv(2*lev_conf-1);
term1    = ca*sqrt(2*t2*ho.^2)/t1;
term2    = t2*ho*(ho-1)/(t1.^2);
qlim = t1*(1+term1+term2).^(1/ho);


%% Model
model.T = T;
model.P = P;
model.U = U;
% model.Q = Q;
% model.W = W;
model.b = B;
model.qlim=qlim;
model.tlim=tlim;
model.Qres=Qres;
model.Tcont=Tcont;
model.Thot=Thot
model.Ypre_cal=Yhat;
% model.rescl=rescl;

% Plotting for each class
figure;
for i = 1:size(Ycal,2)
    classes{i}=find(Ycal(:,i)==1)
plot(model.Thot(classes{i}),model.Qres(classes{i}),'o')
hold on
end




% figure;plot(model.Thot,model.Qres,'o')
% hold on
vline(model.tlim)
hline (model.qlim)
legend ('show')
xlabel ('Hotelling T^2')
ylabel ('Residues in Q')
title ('Calibration Samples')

end
