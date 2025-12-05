clear *
clc
close all

dir_name = './DATA/';

% LOAD ALL DATA FROM AIA EXPERIMENTS AND GET INDIFFERENCE POINTS
k = 1;
s1 = 1;
s2 = 1;

k = 1;
lista = [];
lista = cell(124,4);
for i=1:31
   lista(k,1) = {['CONTROL_TS_',int2str(i),'_A_PRE.csv']};
   lista(k,2) = {'REAL'};
   lista(k,3) = {'PRE'}; 
   lista(k,4) = {sprintf('%02d',i)};k = k+1;
end
for i=1:31
   lista(k,1) = {['CONTROL_TS_',int2str(i),'_A_POST.csv']};
   lista(k,2) = {'REAL'};
   lista(k,3) = {'POST'};
   lista(k,4) = {sprintf('%02d',i)};k = k+1;
end
for i=1:31
   lista(k,1) = {['CONTROL_TS_',int2str(i),'_B_PRE.csv']};
   lista(k,2) = {'SHAM'};
   lista(k,3) = {'PRE'};
   lista(k,4) = {sprintf('%02d',i)};k = k+1;
end
for i=1:31
   lista(k,1) = {['CONTROL_TS_',int2str(i),'_B_POST.csv']};
   lista(k,2) = {'SHAM'};
   lista(k,3) = {'POST'};
   lista(k,4) = {sprintf('%02d',i)};k = k+1;
end

k = 1;
for i=1:size(lista,1)
    A_name = char(lista(i,1));
    if isempty(A_name)==0
        exp_A = import_dd([dir_name A_name]);    
        [exp_A.k_low, exp_A.k_hi, exp_A.r2_low, exp_A.r2_hi] = dd_plot(exp_A.delays, exp_A.value_low, exp_A.value_hi, ['./fig/AIA_sub',char(lista(i,4)),'_',char(lista(i,2)),'_',char(lista(i,3)),''], '','');

        for j=1:size(exp_A.delays,1)
            SUBJ_ID(k) = lista(i,4);
            TIMEPOINT(k) = lista(i,3);
            TES(k) = lista(i,2);
            DEL(k) = exp_A.delays(j);
            REW_LOW(k) = exp_A.value_low(j);
            REW_HI(k) = exp_A.value_hi(j);
            M_RT_LOW(k) = mean(exp_A.response_time_low(j));
            M_RT_HI(k) = mean(exp_A.response_time_hi(j));
            K_LOW(k) = exp_A.k_low;
            K_HI(k) = exp_A.k_hi;
            R2_LOW(k) = exp_A.r2_low;
            R2_HI(k) = exp_A.r2_hi;
            k = k+1;
        end
    end
    
end

DD_TAB_AIA = table(string(SUBJ_ID'),string(TIMEPOINT'),string(TES'),DEL',REW_LOW',REW_HI',M_RT_LOW',M_RT_HI',K_LOW',K_HI',R2_LOW',R2_HI','VariableNames',{'SUBJ_ID','TIMEPOINT','TES','DEL','REW_LOW','REW_HI','M_RT_LOW','M_RT_HI','K_LOW','K_HI','R2_LOW','R2_HI'});

save('DD_TAB_AIA.mat','DD_TAB_AIA');

%%

