function [Ymatrix,ylabel,classes] = vet_matrix(yvector)
% Rotina para transformar um vetor de caracteristicas em matriz de
% caracteristicas - aplicavel em problemas qualitativos.
%
% input
%       yvector  : vetor y com a caracteristica a ser transformada.
%   Caso: yvector seja de string, a sa�da ylabel ser� um vetor n�merico.
%   Caso: yvector seja uma matriz, a sa�da Ymatrix ser� um vetor e ylabel
%         ser� um vetor de string.
% output:
%       Ymatrix : matriz com classes. Cada coluna da Ymatriz refere-se 
% a uma das classes.
%       ylabel  : vetor de classes transformado em string.
%       classes : identifica��o das classes da entrada yvector.
%
% Exemplo:
% [Ymatrix,ylabel,classes] = vet_matrix(yvector);
%
% +++  Paulo R. Filgueiras   - 19/08/2014
%

if isvector(yvector) && isnumeric(yvector)
        [Ymatrix,ylabel,classes] = vet_matrix_vetor(yvector);  % entrada � vetor
    elseif isvector(yvector) && ~isnumeric(yvector)
        [Ymatrix,ylabel,classes] = vet_matrix_string(yvector); % entrada � string
    elseif ~isvector(yvector) && isnumeric(yvector)
        [Ymatrix,ylabel,classes] = vet_matrix_matriz(yvector); % entrada � matriz
end
    

function [Ymatrix,ylabel,classes] = vet_matrix_vetor(yvector)
% subrotina para entrada yvector = vetor.
classes=unique(yvector);  % classes
Ymatrix=[];
ylabel=[];
for ki=1:length(yvector)
    for kj=1:length(classes)
        if yvector(ki)==classes(kj)
            Ymatrix=[Ymatrix;1];
        else
            Ymatrix=[Ymatrix;0];
        end
    end
    ylabel=[ylabel;num2str(yvector(ki))];
end
Ymatrix=reshape(Ymatrix,length(classes),length(yvector));
Ymatrix=Ymatrix';

function [Ymatrix,ylabel,classes] = vet_matrix_string(yvector)
% subrotina para entrada yvector = string.
classes=unique(yvector,'rows');  % classes
Ymatrix=[];
ylabel=[];
for ki=1:size(yvector,1)
    ylabel=[ylabel;str2num(yvector(ki))];
    for kj=1:size(classes,1)
        if str2num(yvector(ki))==str2num(classes(kj))
            Ymatrix=[Ymatrix;1];
        else
            Ymatrix=[Ymatrix;0];
        end
    end
end
Ymatrix=reshape(Ymatrix,length(classes),length(yvector));
Ymatrix=Ymatrix';

function [Ymatrix,ylabel,classes] = vet_matrix_matriz(yvector)
% subrotina para entrada yvector = mstriz.
classes=(1:size(yvector,2))';
Ymatrix=zeros(size(yvector,1),1);  % yvetor com as classes
for kj=1:length(classes)
    aa1=find(yvector(:,kj)==1);
    Ymatrix(aa1)=kj;
end

ylabel=[];   % ylabel
for ki=1:length(Ymatrix)
    ylabel=[ylabel;num2str(Ymatrix(ki))];
end

