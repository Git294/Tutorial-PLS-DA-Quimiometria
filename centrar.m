function [cdata,me,cnewdata]=centrar(data,opt,newdata);
% Rotina que centraliza a matriz de dados
% INPUT
%            x        > Matriz de dados x
%            opt      > 1 centra pela coluna, 2 centra pela linha, 3 centra na m�dia duas vezes 						
%            xnew     > Matriz de dados x da valida��o que ser�o centrados na
%            media utilizando a media das amostras de calibra��o
% OUTPUT
%            xcent    > matriz x centrada na m�dia
%            media    > m�dia da matriz x
%            xnewcent > xnew centrado na m�dia
[m,n]=size(data);

if nargin==1;
  opt=[4];
  while opt>3 | opt<=0 
    opt=input('column centering(1), row centering(2), double centering(3):');
  end
end


if opt==1			% column centering 
   me=mean(data);
   cdata=data-ones(m,1)*me;
end

if opt==2			% row centering
   me=mean(data')';
   cdata=data-me*ones(1,n);
end

if opt==3 	% double centering
   me=mean(mean(data));
   mej=mean(data');
   mei=mean(data);
   cdata=data-(ones(m,1)*mei)-(ones(n,1)*mej)'+(ones(m,n)*me);
end

if exist('newdata')==1			% center new data
    [mt,n]=size(newdata);
    
    if opt==1				% column centering 
        me=mean(data);
        cnewdata=newdata-ones(mt,1)*me;
    else
        error('Row centering and double centering are impossible to perform on a test set');
    end
    
end