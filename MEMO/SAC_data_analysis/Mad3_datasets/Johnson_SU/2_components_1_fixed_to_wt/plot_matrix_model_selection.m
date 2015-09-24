% plots matrix of calculated reduced models sorted according to BIC

clear all
close all

load optimization.mat
load model_sel.mat

TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize = 6;
set(0,TextSizes);

%% Find calculated models in S
index=[];
R_calc=[];
for i=1:length(S)
    if ~isempty(S{1,i}.model)
        if S{1,i}.parameters.MS.MAP.logPost ~= -inf;
            index=[index , i];
        end
    end
end

%%
S_calc=S(index);


R_calc=[[1,2], R(index(2:end)-1)]; % add R [1,2] for full model

%% Sort models according to BIC (ascending, best model first)
sortvec=[];

for i=1:length(S_calc)
    sortvec=[sortvec, S_calc{i}.parameters.BIC];
    
end



[BIC_sort,BIC_sort_ind]=sort (sortvec);



S_calc_sortBIC=S_calc(BIC_sort_ind);
R_calc_sortBIC=R_calc(BIC_sort_ind);




%% reduce number of models for plot
BIC_sort=BIC_sort(1:6);
BIC_sort_ind=BIC_sort_ind(1:6);
S_calc_sortBIC=S_calc_sortBIC(1:6);


%% PLOT table for models ranked according to BIC
%calculations

BIC_dif=BIC_sort- BIC_sort(1,1);

f1 = figure(4);


merker2=0;
merker6=0;
merker10=0;
ah1=subplot(2,1,2)
for i=1:length(S_calc_sortBIC)
    
    
    %Querstreifen
    plot([i-0.5, i-0.5],[0 length(M.experiment)+1],'k-'); hold on;
    exp_matrix_deleted=zeros(length(D),2);
    
    
    
    
    for m=1:length(M.experiment)
        % determine order of experimental conditions for plot
        order=[ 1 2 3 4 5 ];
        
        
        
        if length(S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index)==1
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==1  % kein Wt
                
                aa(1)=plot(i,m+0.05,'ko','Markersize',7); hold on;
                
            end
            
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==2  % keine SAC_pop
                aa(2)=plot(i,m+0.05,'ko','Markersize',7); hold on;
                
                set(aa(2), 'MarkerFaceColor', [0 ,0,0]);
                
                
            end
        else
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==1 && S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(2)==2 % beide POP
                w_plot=S_calc_sortBIC{1,i}.model.experiment(1,order(m)).w_fun{1,2}(S_calc_sortBIC{1,i}.parameters.MS.MAP.par);
                aa(3)=plot(i,m+0.05,'ks','Markersize',8); hold on;
                set(aa(3), 'MarkerFaceColor',[119/255,136/255,153/255]);
                
                
            end
            
        end
    end
    
    hold on;
    
    set(gcf, 'PaperPositionMode', 'auto');
end


hold on

axis equal
ylim([0.5 length(M.experiment)+0.5]);
xlim([0.5 length(S_calc_sortBIC)+0.5]);

pos0 = get(ah1,'Position');
pos0(1)=pos0(1).*0.6;
pos0(2)=pos0(2).*0.6;
pos0(3)=pos0(3).*0.6;
pos0(4)=pos0(4).*0.6;
set(ah1,'Position',pos0)

TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize = 6;
set(0,TextSizes);

for k=1:length(M.experiment)
    lab{k}=D{1,order(k)}.name;
end

set(gca,'yTick',[1:1:length(M.experiment)]);
set(gca,'yticklabel',lab);
set(gca,'tickdir','out');
xlabel('model i' );
pos1 = get(ah1,'Position');

pos1(3)=pos1(3).*0.8;
pos1(4)=pos1(4).*0.8;
set(ah1,'Position',pos1)

% determine plottet area of subplot 1
pos = plotboxpos(gca);

%% subplot 2 BIC

ah2=subplot(2,1,1)
plot([1:length(BIC_sort)], BIC_sort,'kx');
hold on
xlim([0.5 length(BIC_sort)+0.5]);
ylim([min(BIC_sort)-20 5.3858e+03]);

ylabel('BIC of model i');



%legend
h=legend(aa,'one population not WT','one population WT','two subpopulations','Location','Best');

%# find current position [x,y,width,height]
pos2 = get(ah2,'Position');
pos1 = get(ah1,'Position');

%# set width of second axes equal to first
pos2(3) = pos(3);
pos2(4) = pos1(4);
set(ah2,'Position',pos2)
set(ah2,'XTickLabel','')
pos2(2) = pos(2) + pos(4);
pos2(1) = pos(1);
set(ah2,'Position',pos2)




