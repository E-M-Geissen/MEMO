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
BIC_sort=BIC_sort(1:8);
BIC_sort_ind=BIC_sort_ind(1:8);
S_calc_sortBIC=S_calc_sortBIC(1:8);
% R_calc_sortBIC=R_calc_sortBIC(1:8);

%% PLOT table for models ranked according to BIC
%calculations

BIC_dif=BIC_sort- BIC_sort(1,1);

f1 = figure(4);
%screen_size = get(0, 'ScreenSize');

%set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );


merker2=0;
merker6=0;
merker10=0;
ah1=subplot(2,1,2)
for i=1:length(S_calc_sortBIC)
    
    %if i~=1 
    %Querstreifen
    plot([i-0.5, i-0.5],[0 length(M.experiment)+1],'k-'); hold on;
    exp_matrix_deleted=zeros(length(D),2);
    %end
    
    
    
    for m=1:length(M.experiment)
        order=[ 4 3  2 1 8 7 6 5 9];
       % order=[ 9 5 6 7 8 1 2 3 4 ];
        
        
        if length(S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index)==1
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==2  % kein Wt
                %parindex=parindex+delta_parindex-1;
                aa(1)=plot(i,m+0.05,'ko','Markersize',7); hold on;
                %aa(1)=plot(i,m+0.05,'ko'); hold on;
                %set(aa(4), 'MarkerFaceColor', [0.5 ,0.5,0]);
                %set(aa(1), 'MarkerFaceColor', [0 ,0,0]);
            end
            
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==1  % keine SAC_pop
                aa(2)=plot(i,m+0.05,'ko','Markersize',7); hold on;
                %aa(2)=plot(i,m+0.05,'ko'); hold on;
                set(aa(2), 'MarkerFaceColor', [0 ,0,0]);
                %set(h, 'MarkerFaceColor', [0.5 ,0.5,0]);
                
            end
        else
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==1 && S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(2)==2 % beide POP
                w_plot=S_calc_sortBIC{1,i}.model.experiment(1,order(m)).w_fun{1,2}(S_calc_sortBIC{1,i}.parameters.MS.MAP.par);
                aa(3)=plot(i,m+0.05,'ks','Markersize',8); hold on;
                set(aa(3), 'MarkerFaceColor',[119/255,136/255,153/255]);
                %aa(3)=plot(i,m+0.05,'ks'); hold on;
                %                 rot=1;
                %                 gruen=1-w_plot;
                %                 blau =1-w_plot;
                rot=0;
                gruen=0;
                blau =1;
                %set(aa(5), 'MarkerFaceColor', [1 0 1 ]);
                %set(aa(3), 'MarkerFaceColor', [rot gruen blau]);
            end
            
        end
    end
    
    hold on;
    %
    %     if  BIC_dif(i) >=2 && merker2==0;
    %         aa(1)=plot([(i-1)+0.5, (i-1)+0.5],[0 length(M.experiment)+1],'k-','LineWidth', 4);
    %         merker2=1;
    %     end
    %
    %         if  BIC_dif(i) >=6 && merker6==0;
    %         aa(2)=plot([(i-1)+0.5, (i-1)+0.5],[0 length(M.experiment)+1],'k--','LineWidth', 4);
    %         merker6=1;
    %         end
    
    %         if  BIC_dif(i) >=10 && merker10==0;
    %         aa(3)=plot([(i-1)+0.5, (i-1)+0.5],[0 length(M.experiment)+1],'k-.','LineWidth', 4);
    %         merker10=1;
    %         end
    set(gcf, 'PaperPositionMode', 'auto');
end
%box off
%axes(gca,'Linewidth',5);
%set(gca, 'visible', 'off') ; 
%axis off
%set(gca,'ytick',[],'ycolor',get(ah1,'color'))
%set(gca,'xtick',[],'xcolor',get(ah1,'color'))


%%
%plot([i+0.5, i+0.5],[0 length(S_calc_sortBIC)+1],'k-');
hold on
%pos1 = get(ah1,'Position');
%pos1(3)=0.5
% set(ah1,'Position',pos1)
axis equal
ylim([0.5 length(M.experiment)+0.5]);
xlim([0.5 length(S_calc_sortBIC)+0.5]);

%%
TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize = 6;
set(0,TextSizes);

for k=1:length(M.experiment)
    lab{k}=D{1,order(k)}.name;
end
lab{1}='200% Mad2';
set(gca,'yTick',[1:1:length(M.experiment)]);
set(gca,'yticklabel',lab);
%set(gca,'XAxisLocation','top')
set(gca,'tickdir','out');
xlabel('model i' );
%ylabel('experiment');



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
ylim([min(BIC_sort)-30 max(BIC_sort)+10]);
ylim([8185 8225])
ylabel('BIC of model i');

% merker für senrechte LSignifikanzlinien
merker2=0;
merker6=0;
merker10=0;
for i=1:length(S_calc_sortBIC)
    %     if  BIC_dif(i) >=2 && merker2==0;
    %        aa(1) =plot([(i-1)+0.5, (i-1)+0.5],[0 17765],'k-','LineWidth', 4);
    %         merker2=1;
    %     end
    %
    %         if  BIC_dif(i) >=6 && merker6==0;
    %        aa(2) =plot([(i-1)+0.5, (i-1)+0.5],[0 17765],'k--','LineWidth', 4);
    %         merker6=1;
    %         end
    
    %         if  BIC_dif(i) >=10 && merker10==0;
    %        aa(3) =plot([(i-1)+0.5, (i-1)+0.5],[0 max(BIC_sort)+10],'k-.','LineWidth', 4);
    %         merker10=1;
    %         end
end

%%
%grid on
%h=legend(aa,'one population not WT','one population WT','two subpopulations','Location','Best');
%v = get(h,'title');
%set(v,'string','Legend Title');

%title( ['model EJ CJ vfv in EMM'])
%# find current position [x,y,width,height]
pos2 = get(ah2,'Position');
pos1 = get(ah1,'Position');

% %set height of ah1 to have same space in columns and rows
% pos1(4)=(pos1(3)/length(BIC_sort))*length(M.experiment)+0.1;
% set(ah1,'Position',pos1)

%# set width of second axes equal to first
%pos2(3) = pos1(3);
pos2(3) = pos(3);
pos2(4) = 0.5*pos1(4);
set(ah2,'Position',pos2)
%
set(ah2,'XTickLabel','')
%pos1(2) = pos2(2) - pos1(4);
pos2(2) = pos(2) + pos(4);
%pos1(1) = pos2(1);
pos2(1) = pos(1);

%set(ah1,'Position',pos1)
set(ah2,'Position',pos2)


%%
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 4 8])
print('-depsc','./BIC_ranking_bis8_WT');

