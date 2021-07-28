%% Routine developed according to the article published in Quimica nova
% Download all the necessary functions from:
% https://github.com/felipebachion/Tutorial-PLS-DA-Quimiometria
% The functions and the sample bank can also be requested via e-mail:
% rjpoppi@unicamp.br

%% 1) loading the dataset
load amostras
%% 2) plot complete set of spectra
figure(1)
subplot(2,1,1)
plot(num,X);
axis tight
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A','LineStyle','none');
subplot(2,1,2)
plot(num(:,1300:1750),X(:,1300:1750));
axis tight
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.3,0.1,0.1],'String', 'B','LineStyle','none');
[ax,h1]=suplabel('Wavenumber (cm^-^1)','x',[.11 .11 .84 .84]); 
[ax,h2]=suplabel('Absorbance','y'); 
set(h1,'FontSize',24,'FontName','Times New Roman') 
set(h2,'FontSize',24,'FontName','Times New Roman') 
%% 3) create the class vector to segregate the spectra of each type of oil
y=[ones(108,1);2*ones(54,1);3*ones(54,1);4*ones(54,1)];
%% 4) use the CALTESTDA algorithm to separate the samples from the calibration, test and preprocess sets
% observe that a region slightly larger than the fingerprint region is
% selected.
[~,xcal,xval,ycal,yval]=caltestda(X(:,1300:1750),y,70,'k',[],[]);
numc=num(1300:1750); % Performing the spectral band cutoff.
%% 5) Plot the calibration and validation spectra after selecting the spectral region of interest
figure(2)
subplot(2,1,1);
plot(numc,xcal)
axis tight
set(gca,'xtick',[],'FontSize',18);...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A');
title('Calibration Samples','fontsize',14);
subplot(2,1,2);
plot(numc,xval)
axis tight
set(gca,'xtick',[],'FontSize',18);...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.3,0.1,0.1],'String', 'B');
title('Validation Samples','fontsize',14); hold on
[ax,h1]=suplabel('Wavenumber (cm^-^1)','x',[.11 .11 .84 .84]); 
[ax,h2]=suplabel('Absorbance','y',[.11 .11 .84 .84]); 
set(h1,'FontSize',24) 
set(h2,'FontSize',24) 
%% 6) Use the function vet_matrix to transform a vector y1 into a matrix y2
%   y1[1        y2[1 0 0
%      1           1 0 0
%      2           0 1 0
%      2           0 1 0
%      3           0 0 1
%      3]          0 0 1]

ycal = vet_matrix(ycal);
yval = vet_matrix(yval);
% type ycal and yval in the prompt of the command window to view the
% creation of the class matrix.
%% 7) Use the PRETRAT routine for the spectra plotting; Remember that I did not focus on the mean
[xcal,xval]=pretrat(xcal,xval,{'deriv';[9,2,1]});

%% 8) plot the preprocessed calibration and validation matrices for
% visualisation of the spectral preprocessing employed.
figure(3)
subplot(2,1,1);
plot(numc,xcal)
axis tight
title('Calibration Samples','fontname','Times New Roman','fontsize',14);
hold on
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A','LineStyle','none');
subplot(2,1,2);
plot(numc,xval)
axis tight
title('Validation Samples','fontname','Times New Roman','fontsize',14); hold on
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.3,0.1,0.1],'String', 'B','LineStyle','none');
[ax,h1]=suplabel('Wavenumber (cm^-^1)','x',[.11 .11 .84 .84]); 
[ax,h2]=suplabel('First Derivative','y',[.08 .08 .84 .84]); 
set(h1,'FontSize',24,'FontName','Times New Roman') 
set(h2,'FontSize',24,'FontName','Times New Roman') 
%% 9) Cross-validation to determine the Number of Latent Variables
cvvc = my_cross_validation(xcal,ycal,10,10,4,0.8);

figure(4)
plot(cvvc.porc_am_class_cor')
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
legend ('Class 1', 'Class 2', 'Class 3', 'Class 4')
xlabel('Number of Latent Variables','FontSize',24,'FontName','Times New Roman');
ylabel('Percentage of correctly classified samples','FontSize',24,'FontName','Times New Roman');
disp(['Choose the Number of Latent Variables ','. ']);
A=input('');

%% 10) Analysing the presence of anomalous samples for the number of latent variables chosen
[model] = my_calc_qt_limits_cal (xcal,ycal,xval,A);
% Anomalous Samples
am_anomalasQ=find(model.Qres>=model.qlim);
am_anomalasT2=find(model.Thot>=model.tlim);
a=ismember(am_anomalasQ,am_anomalasT2);
am_anomalas=am_anomalasQ(1,a);

%% 11) Eliminating these samples and selecting again the optimal number of latent variables.
xcal(am_anomalas,:)=[];
ycal(am_anomalas,:)=[];

cvvc = my_cross_validation(xcal,ycal,10,10,4,0.8);

figure
plot(cvvc.porc_am_class_cor')
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
legend ('Classe 1', 'Classe 2', 'Classe 3', 'Classe 4')
xlabel('Número de variáveis latentes','FontSize',24,'FontName','Times New Roman');
ylabel('Porcentagem de amostras classificadas corretamente','FontSize',24,'FontName','Times New Roman');
disp(['Escolha o Número de variáveis latentes ','. ']);
A=input('');

%% 12) Analisando novamente a presença de amostras anômalas para o novo modelo PLS-DA nos conjuntos de calibração e validação.
[model_cal] = my_calc_qt_limits_cal (xcal,ycal,xval,A);
[model_val] = my_calc_qt_limits_val (xcal,ycal,xval,A);
% Desta vez não temos amostras anômalas no conjunto de calibração e
% validação

%% 13) Calcular o modelo final PLS-DA
[yprev_cal]=previsto_pls(xcal,ycal,xcal,0,A);
% Calculando o threshould
ts=[];
for ki=1:size(ycal,2)
    plsda_thres = plsdafindthr(yprev_cal(:,ki),ycal(:,ki));
    ts=[ts,plsda_thres.class_thr];
end
% Prevendo as amostras
for u = 1:size(ycal,2)
    yprev_calts(:,u)=yprev_cal(:,u)>=ts(:,u);
end
%% 14) Calcular os valores de previsão das amostras de validação
[yprev_val]=previsto_pls(xcal,ycal,xval,0,A);
for u = 1:size(yval,2)
    yprev_valts(:,u)=yprev_val(:,u)>=ts(:,u);
end

%% 15) Plotar o gráfico das amostras de calibração
Nome{1,1}='Oliva' ;  Nome{1,2} = 'k';
Nome{2,1}='Canola' ; Nome{2,2} = 'm';
Nome{3,1}='Milho' ;  Nome{3,2} = 'r';
Nome{4,1}='Soja' ;   Nome{4,2} = 'g';
figure(5)
for ki=1:size(yprev_cal,2)
    %figure(ki)
    subplot(2,2,ki)
    tp1=find(ycal(:,ki));tp2=setxor(1:length(ycal),tp1);
    plot(1:length(yprev_cal),yprev_cal(:,ki),'o','MarkerFaceColor','b'), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',14,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
    marc=strcat(Nome{ki,2},'o');
    plot(tp1,yprev_cal(tp1,ki),marc,'MarkerFaceColor',marc(:,1)), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',14,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
        hline(ts(ki),'r')
    title(Nome{ki,1},'FontSize',14,'FontName','Times New Woman')
    xlabel('Amostra','FontSize',18,'FontName','Times New Woman')
    ylabel(sprintf('Classe %g',ki),'FontSize',18,'FontName','Times New Woman')
    end
%% 16) Plotar o gráfico das amostras de validação
figure(6)
for ji=1:size(yprev_val,2)
    %figure(ji)
    subplot(2,2,ji)
    tv1=find(yval(:,ji));tv2=setxor(1:length(yval),tv1);
    plot(1:length(yprev_val),yprev_val(:,ji),'o','MarkerFaceColor','b'), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',14,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
    marc=strcat(Nome{ji,2},'o');
    plot (tv1,yprev_val(tv1,ji),marc,'MarkerFaceColor',marc(:,1)), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',12,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
    hline(ts(ji),'r')
    title(Nome{ji,1},'FontSize',14,'FontName','Times New Woman')
    xlabel('Amostra','FontSize',18,'FontName','Times New Woman')
    ylabel(sprintf('Classe %g',ji),'FontSize',18,'FontName','Times New Woman')
end
    
%% 17 Calculando a tabela de confusão para o conjunto de Calibração

[Tcal Tcal2] = my_ConfTable_cal (ycal,yprev_cal,ts)
[Tval Tval2] = my_ConfTable_val (yval,yprev_val,ts)

clearvars -except model_cal model_val class_param_cal class_param_val num numc xcal xval ycal yval YCAL YVAL yprev_cal yprev_cal_new yprev_calts yprev_val yprev_val_new yprev_valts Tcal2 Tcal Tval2 Tval
