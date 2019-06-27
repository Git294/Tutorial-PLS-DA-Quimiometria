function [xnew]=centrar_inverso(x,media);
% Rotina que realiza a Centraliza��o inversa (Inverse centering)
% INPUT 
%         X     > matriz X
%         media > m�dia da matriz X
% OUTPUT
%         xnew  > Matriz apos a centraliza��o inversa
[m,n]=size(x);
[mme,nme]=size(media);


if nme==n			
  xnew=x+ones(m,1)*media;
end

if mme==m			
  xnew=x+media*ones(1,n);
end

if size(media)==1		
  xnew=x+ones(m,n)*media;
end


  
