function [cvv] = my_cross_validation(X,y,comp,cv_groups,n_class,perc_in);
% Rotina que realiza a valida��o cruzada selecionando aleat�riamente 80%
% das amostras de cada classe e com o restante das amostras realizada a
% valida��o cruzada.
% X � conjunto de calibracao
% y s�o os valores de classe da forma y = [1 0 0..
%                                          [0 1 0..]
% comp � o n�mero m�ximo de variaveis que ser�o testadas
% n_class � o n�mero de classes que contem o matriz Y.
% cv_groups � o n�mero de vezes que as amostras ser�o aleat�rizadas para
% fazer a valida��o cruzada.
% perc_in � 1-% de amostras que ser�o utilizadas na valida��o cruzada
% Rotina desenvolvida por Felipe Bachion para fins n�o comerciais
% contato: felipebachion@gmail.com

if nargin==5; perc_in=0.8; end

classes=[];

    h=waitbar(0,'Construindo os modelos PLS-DA');
for n=1:comp %1
    h=waitbar(n/comp);
       for i=1:cv_groups
        
        [x_in,y_in,x_out,y_out] = myPerm (X,y,n_class,perc_in);
        
        [yprev_val_cruzada]=previsto_pls(x_in,y_in,x_out,0,n);
        [yprev_cal_cruzada]=previsto_pls(x_in,y_in,x_in,0,n);
        cvv.pred_valcruzada(:,i,:)=yprev_val_cruzada;
        cvv.ref_valcruzada(:,i,:)=y_out;
        cvv.pred_calcruzada(:,i,:)=yprev_cal_cruzada;
        cvv.ref_calcruzada(:,i,:)=y_in;
         % calculo do threshould para cada vl nas amostras de valida��o
         % cruzada e j� encontrando o n�mero de amostras classificadas
         % corretamente
                for ki=1:size(cvv.pred_valcruzada,2)
                    for ko=1:size(y,2);
                plsda_thres=plsdafindthr(cvv.pred_calcruzada(:,ki,ko),cvv.ref_calcruzada(:,ki,ko));
                ts (n,ki,ko)=[plsda_thres.class_thr];
                     cv=yprev_val_cruzada;
                     cv=cv>=ts(n,ki,ko);
                     cv=double(cv);
                     aw1=cv==y_out;
                     aw=sum(aw1);
                     n_correc_class_val_cruzada(n,ki,ko)=sum(aw);
                    end
                end
                
     cvv.ts=ts;
     cvv.perc_correc_class_val_cruzada=100*(n_correc_class_val_cruzada/(n_class*size(y_out,1)));


        cvv.pred_valcruzada(:,i,:)=yprev_val_cruzada;
        cvv.ref_valcruzada(:,i,:)=y_out;
        cvv.pred_calcruzada(:,i,:)=yprev_cal_cruzada;
        cvv.ref_calcruzada(:,i,:)=y_in;
         % calculo do threshould para cada vl nas amostras de valida��o
         % cruzada e j� encontrando o n�mero de amostras classificadas
         % corretamente
                for ki=1:size(cvv.pred_valcruzada,2)
                    for ko=1:size(y,2);
                plsda_thres=plsdafindthr(cvv.pred_calcruzada(:,ki,ko),cvv.ref_calcruzada(:,ki,ko));
                ts (n,ki,ko)=[plsda_thres.class_thr];
                     cv=yprev_val_cruzada;
                     cv=cv>=ts(n,ki,ko);
                     cv=double(cv);
                     aw1=cv==y_out;
                     aw=sum(aw1);
                     n_correc_class_val_cruzada(n,ki,ko)=sum(aw);
                    end
                end
     cvv.ts=ts;
     cvv.perc_correc_class_val_cruzada=100*(n_correc_class_val_cruzada/(n_class*size(y_out,1)));
    end

end
close (h)

    for ll=1:comp
        for u=1:n_class
cvv.porc_am_class_cor(u,ll)=mean(cvv.perc_correc_class_val_cruzada(ll,:,u));
        end
        
    end

    end