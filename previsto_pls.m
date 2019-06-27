function [Ynew_hat]=previsto_pls(Xcal,Ycal,Xnew,standard,nvl);
% INPUT =
% Xcal = matriz calibra��o (xcal)
% Ycal = Vetor de calibra��o (ycal)
% Xnew =  matriz que ser� prevista Ex: xcal, xval
% standard = 1 para autoescalar os dados e 0 para centrar os dados na m�dia
% nvl = n�mero de variaveis latentes

% OUTPUT = Ynew_hat > valores previstos de Y para o conjunto Xnew

% fun��o para construir o modelo PLS para um dado n�mero de variavies
% latentes e obter os resultados de previs�o de Y para o conjunto Xnew
% Caso queira os resultados de calibra��o utilizo o xcal novamente, exemplo
% para dados centrados na media e n�mero de variaveis latentes = 5
% [ycal_pred]=plspred(xcal,ycal,xcal,0,5);

[n,px]=size(Xcal);   					           
[n,m]=size(Ycal); 

if standard==1						% Caso seja selecionando o autoscalamento
	[Xauto,mX,sX,Xnew_auto]=scale(Xcal,1,Xnew);
	[Yce,mY]=centrar(Ycal,1);				
	[B]=simpls(Xauto,Yce,nvl,Xauto'*Yce,[]);	
	Ynew_hat=Xnew_auto*B;				
	[Ynew_hat]=centrar_inverso(Ynew_hat,mY);
end

if standard==0						
	[Xce,mX,Xnew_ce]=centrar(Xcal,1,Xnew);		
	[Yce,mY]=centrar(Ycal,1);
	[B]=simpls(Xce,Yce,nvl,Xce'*Yce,[]);	        
	Ynew_hat=Xnew_ce*B;			
	[Ynew_hat]=centrar_inverso(Ynew_hat,mY);
end
end
