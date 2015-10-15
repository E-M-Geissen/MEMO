clear all
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
%index=ind_S;
S_calc=S(index);

%S_calc=S;
R_calc=[[1,2], R(index(2:end)-1)]; % add R [1,2] for full model

%% Sort models according to BIC (ascending, best model first)
sortvec=[];

for i=1:length(S_calc)
    sortvec=[sortvec, S_calc{i}.parameters.BIC];
    
end



[BIC_sort,BIC_sort_ind]=sort (sortvec);

S_calc_sortBIC=S_calc(BIC_sort_ind);
R_calc_sortBIC=R_calc(BIC_sort_ind);

%%
for i=1:length(S_calc_sortBIC)
    for m=1:length(M.experiment)
        if length(S_calc_sortBIC{1,i}.model.experiment(1,m).component_index)==2
            C_index{i}(m,1)=S_calc_sortBIC{1,i}.model.experiment(1,m).component_index(1);
            C_index{i}(m,2)=S_calc_sortBIC{1,i}.model.experiment(1,m).component_index(2);
        else
            if S_calc_sortBIC{1,i}.model.experiment(1,m).component_index(1)==1
                C_index{i}(m,1)=1;
                C_index{i}(m,2)=0;
            else
                C_index{i}(m,1)=0;
                C_index{i}(m,2)=2;
            end
            
        end
    end
end


i=1;
m=1;

while  i<=length(S_calc_sortBIC)
    m=i+1;
    while m<=length(S_calc_sortBIC)
        if  isequal(C_index{1,i},C_index{1,m})
            S_calc_sortBIC(m)=[];
            C_index(m)=[];
            BIC_sort(m)=[];
            m=m-1;
        end
        m=m+1;
    end
    
    i=i+1;
end


%% reduce number of models for plot
BIC_sort=BIC_sort(1:6);
BIC_sort_ind=BIC_sort_ind(1:6);
S_calc_sortBIC=S_calc_sortBIC(1:6);


%% PLOT table for models ranked according to BIC
%calculations

BIC_dif=BIC_sort- BIC_sort(1,1);

f1 = figure(4);

ah1=subplot(2,1,2)
for i=1:length(S_calc_sortBIC)
    
    
    plot([i-0.5, i-0.5],[0 length(M.experiment)+1],'k-'); hold on;
    exp_matrix_deleted=zeros(length(D),2);
    
    
    
    
    for m=1:length(M.experiment)
        order=[ 4 3 2 1 5 ];
        
        
        
        if length(S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index)==1
            
            aa(1)=plot(i,m+0.05,'ko','Markersize',7); hold on;
            
            
            
        else
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==1 && S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(2)==2 % beide POP
                w_plot=S_calc_sortBIC{1,i}.model.experiment(1,order(m)).w_fun{1,2}(S_calc_sortBIC{1,i}.parameters.MS.MAP.par);
                aa(2)=plot(i,m+0.05,'ks','Markersize',8); hold on;
                set(aa(2), 'MarkerFaceColor',[119/255,136/255,153/255]);
                
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

%%
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
ylim([min(BIC_sort)-10 max(BIC_sort)+10]);
ylim([  5.1100e+03 max(BIC_sort)+20])
ylabel('BIC of model i');


h=legend(aa,'one subpopulation','two subpopulations','Location','Best');


pos2 = get(ah2,'Position');
pos1 = get(ah1,'Position');



%# set width of second axes equal to first

pos2(3) = pos(3);
pos2(4) = pos1(4);
set(ah2,'Position',pos2)
%
set(ah2,'XTickLabel','')

pos2(2) = pos(2) + pos(4);

pos2(1) = pos(1);


set(ah2,'Position',pos2)




