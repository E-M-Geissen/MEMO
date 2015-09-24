%clear all
%model_Mad_all_ElogN_CJ_EMM_vv;
%load ./all.mat % R S D

%%

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


R_calc=[[7,2], R(index(2:end)-1)]; % add R [1,2] for full model 

%% Sort models according to BIC (ascending, best model first)
sortvec=[];

for i=1:length(S_calc)
sortvec=[sortvec, S_calc{i}.parameters.BIC];
end



[BIC_sort,BIC_sort_ind]=sort (sortvec);



S_calc_sortBIC=S_calc(BIC_sort_ind);
R_calc_sortBIC=R_calc(BIC_sort_ind);

%% Delete models that are the same
m=0;
i=1;
while  i<=length(S_calc_sortBIC)
test=R_calc_sortBIC{1,i}(:,1);
m=i+1;
while m<=length(S_calc_sortBIC)
if length(R_calc_sortBIC{1,m}(:,1)) == length(test)
    if R_calc_sortBIC{1,m}(:,1) == test
        R_calc_sortBIC(m) =[];
        S_calc_sortBIC(m)=[];
        BIC_sort(m)=[];
        m=m-1;
    end
end
    m=m+1;
        
end
i=i+1;

end


%% reduce number of models for plot
BIC_sort=BIC_sort(1:8);
BIC_sort_ind=BIC_sort_ind(1:8);
 S_calc_sortBIC=S_calc_sortBIC(1:8);
R_calc_sortBIC=R_calc_sortBIC(1:8);

%% PLOT table for models ranked according to BIC



%calculations

BIC_dif=BIC_sort- BIC_sort(1,1);

f1 = figure(4);
screen_size = get(0, 'ScreenSize');

set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );


merker2=0;
merker6=0;
merker10=0;
 ah1=subplot(2,1,2)
for i=1:length(S_calc_sortBIC)
    
    %Querstreifen
    plot([i-0.5, i-0.5],[0 length(M.experiment)+1],'k-'); hold on;
    exp_matrix_deleted=zeros(length(D),2);
    
%     idx = sub2ind(size(exp_matrix_deleted), R_calc_sortBIC{1,i}(:,1), R_calc_sortBIC{1,i}(:,2));
%     exp_matrix_deleted(idx)=1;
%     exp_matrix_incl=exp_matrix-exp_matrix_deleted;
    

    
    for m=1:length(M.experiment)
        
    order=[ 4 3  2 1 8 7 6 5 9];
            
            
    if length(S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index)==1
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==2  % kein Wt
                %parindex=parindex+delta_parindex-1;
                aa(4)=plot(i,m+0.05,'ko','Markersize',13); hold on;
                set(aa(4), 'MarkerFaceColor', [0 ,0,0]);
            end

            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==1  % keine SAC_pop
                h=plot(i,m+0.05,'ko','Markersize',13); hold on;
                set(h, 'MarkerFaceColor', [0 ,0,0]);
            end
    else
            if S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(1)==1 && S_calc_sortBIC{1,i}.model.experiment(1,order(m)).component_index(2)==2 % beide POP
               w_plot=S_calc_sortBIC{1,i}.model.experiment(1,order(m)).w_fun{1,2}(S_calc_sortBIC{1,i}.parameters.MS.MAP.par);
                aa(5)=plot(i,m+0.05,'ks','Markersize',15); hold on;
                rot=1;
                gruen=1-w_plot;
                blau =1-w_plot;
                set(aa(5), 'MarkerFaceColor', [1 1 1 ]);
            end
  
    end
    end
    
    hold on;
    
    if  BIC_dif(i) >=2 && merker2==0;
        aa(1)=plot([(i-1)+0.5, (i-1)+0.5],[0 length(M.experiment)+1],'k-','LineWidth', 4);
        merker2=1;
    end
    
        if  BIC_dif(i) >=6 && merker6==0;
        aa(2)=plot([(i-1)+0.5, (i-1)+0.5],[0 length(M.experiment)+1],'k--','LineWidth', 4);
        merker6=1;
        end
        
        if  BIC_dif(i) >=10 && merker10==0;
        aa(3)=plot([(i-1)+0.5, (i-1)+0.5],[0 length(M.experiment)+1],'k-.','LineWidth', 4);
        merker10=1;
        end
    set(gcf, 'PaperPositionMode', 'auto');
end
%%
plot([i+0.5, i+0.5],[0 length(S_calc_sortBIC)+1],'k-');
hold on
%pos1 = get(ah1,'Position');
 %pos1(3)=0.5
% set(ah1,'Position',pos1)
axis equal
ylim([0.5 length(M.experiment)+0.5]);
xlim([0.5 length(S_calc_sortBIC)+0.5]);

%%
for k=1:length(M.experiment)
    lab{k}=D{1,order(k)}.name;
end
set(gca,'yTick',[1:1:length(M.experiment)]);
set(gca,'yticklabel',lab);
%set(gca,'XAxisLocation','top')
 set(gca,'tickdir','out');
xlabel('model i' );
ylabel('experiment');



% colorbar für SAC-Anteil bei 2 Populationen
% 
% %colormap für colorbar
% mycmap=zeros([11 3]);
% mycmap(:,1)=1;
% mycmap(:,2)=1-[0:0.1:1];
% mycmap(:,3)=1-[0:0.1:1];
% mycmap=zeros([21 3]);
% mycmap(:,1)=1;
% mycmap(:,2)=1-[0:0.05:1];
% mycmap(:,3)=1-[0:0.05:1];
% 
% colormap(mycmap);
% %Ylim([0 20])
% hcb=colorbar;
% %hcb=lcolorbar(labels)
% set(hcb,'YTick',[1 22],'YTickLabel',{'0','100'}) 
% ylabel(hcb,'% lower population') 
% labels={'0% lower population','100% lower population'};



% determine plottet area of subplot 1 
pos = plotboxpos(gca); 
%% subplot 2 BIC

 ah2=subplot(2,1,1)
 plot([1:length(BIC_sort)], BIC_sort,'kx');
 hold on
xlim([0.5 length(BIC_sort)+0.5]);
ylim([min(BIC_sort)-30 max(BIC_sort)+10]);
ylabel('BIC of model i');

% merker für senrechte LSignifikanzlinien
merker2=0;
merker6=0;
merker10=0;
for i=1:length(S_calc_sortBIC)
    if  BIC_dif(i) >=2 && merker2==0;
       aa(1) =plot([(i-1)+0.5, (i-1)+0.5],[0 max(BIC_sort)+10],'k-','LineWidth', 4);
        merker2=1;
    end
    
        if  BIC_dif(i) >=6 && merker6==0;
       aa(2) =plot([(i-1)+0.5, (i-1)+0.5],[0 max(BIC_sort)+10],'k--','LineWidth', 4);
        merker6=1;
        end
        
        if  BIC_dif(i) >=10 && merker10==0;
       aa(3) =plot([(i-1)+0.5, (i-1)+0.5],[0 max(BIC_sort)+10],'k-.','LineWidth', 4);
        merker10=1;
        end
end
%ylim([17692 17765])
%ylim([16130 16220])
%%
%grid on
h=legend(aa,'Left of line delta BIC < 2 (not worth mentioning)', 'Left of line delta BIC < 6 (weak)','Left of line delta BIC < 10 (strong)','one population','two subpopulations','Location','Best');   v = get(h,'title');
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


% 

 print('-depsc','./BIC_ranking_Mad2_vv');
 
%  %%
%  figure
%   plot([1:length(BIC_sort)], BIC_sort,'kx');
%   hold on;
% xlim([0.5 11.8]);
% ylim([BIC_sort(1)-1 BIC_sort(11)+1])
% ylabel('BIC');
% xlabel('model');
% set(gca,'xTick',[1:1:11]);
% merker2=0;
% merker6=0;
% merker10=0;
% for i=1:length(S_calc_sortBIC)
%     if  BIC_dif(i) >=2 && merker2==0;
%         aa(1)=plot([(i-1)+0.5, (i-1)+0.5],[0 max(BIC_sort(1:11))+2],'k-','LineWidth', 4);
%         merker2=1;
%     end
%     
%         if  BIC_dif(i) >=6 && merker6==0;
%          aa(2)=plot([(i-1)+0.5, (i-1)+0.5],[0 max(BIC_sort(1:11))+2],'k--','LineWidth', 4);
%         merker6=1;
%         end
%         
%         if  BIC_dif(i) >=10 && merker10==0;
%          aa(3)=plot([(i-1)+0.5, (i-1)+0.5],[0 max(BIC_sort(1:11))+2],'k-.','LineWidth', 4);
%         merker10=1;
%         end
% end
% %h=legend(aa,'Left of line delta BIC < 2 (not worth mentioning)', 'Left of line delta BIC < 6 (weak)','Left of line delta BIC < 10 (strong)','Location','Best');
% ylim([min(BIC_sort(1:11))-2 max(BIC_sort(1:11))+2]);
%  print('-depsc','./zoom');
% %set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11'}) 
% %% Sort models according to AIC (ascending, best model first)
% sortvec=[];
% 
% for i=1:length(S_calc)   
% sortvec=[sortvec, S_calc{i}.parameters.AIC];
% end
% 
% [AIC_sort,AIC_sort_ind]=sort (sortvec);
% 
% S_calc_sortAIC=S_calc(AIC_sort_ind);
% R_calc_sortAIC=R_calc(AIC_sort_ind);
% 
% %% PLOT table for mdels ranked according to AIC
% exp_matrix=[1 0; 1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1];
% exp_matrix_deleted=[0 0; 0 0;0 0;0 0;0 0;0 0;0 0;0 0; 0 0;0  0;0 0];
% 
% 
% 
% 
% screen_size = get(0, 'ScreenSize');
% f2 = figure(2);
% set(f2, 'Position', [0 0 screen_size(3) screen_size(4)/2 ] );
% % maximal kann man über titelbar gehen
% 
% 
% 
% 
% for i=1:length(S_calc_sortAIC)
%     plot([i-0.5, i-0.5],[0 length(S_calc_sortAIC)+1],'k-'); hold on;
%     exp_matrix_deleted=[0 0; 0 0;0 0;0 0;0 0;0 0;0 0;0 0; 0 0;0  0;0 0];
%     
%     idx = sub2ind(size(exp_matrix_deleted), R_calc_sortAIC{1,i}(:,1), R_calc_sortAIC{1,i}(:,2));
%     exp_matrix_deleted(idx)=1;
%     exp_matrix_incl=exp_matrix-exp_matrix_deleted;
%     
%     parindex=0;
%             switch (S{1,1}.model.mixture.type)
%         
%                 case 'log-normal'
%                     parindex=2;
%                 case 'johnson'
%                     parindex=4;
%             end
%     
%     for m=1:11
%         
%             switch (S_calc_sortAIC{1,i}.model.mixture.type)
%         
%                 case 'log-normal'
%                     delta_parindex=3;
%                 case 'johnson'
%                     delta_parindex=5;
%             end
%             
%         if m==1
%             h=plot(i,m+0.05,'k^','Markersize',15); hold on;
%         else
%     
%             if (exp_matrix_incl(m,1)==0) && (exp_matrix_incl(m,2) ==1) % kein Wt
%                 parindex=parindex+delta_parindex-1;
%                 h=plot(i,m+0.05,'k^','Markersize',15); hold on;
%                 set(h, 'MarkerFaceColor', [1 ,0,0]);
%             end
% 
%             if (exp_matrix_incl(m,1)==1) && (exp_matrix_incl(m,2) ==0) % keine SAC_pop
%                 h=plot(i,m+0.05,'k^','Markersize',15); hold on;
%             end
% 
%             if (exp_matrix_incl(m,1)==1) && (exp_matrix_incl(m,2) ==1) % beide POP
%                 parindex=parindex+delta_parindex; 
%                 h=plot(i,m+0.05,'k^','Markersize',15); hold on;
%                 rot=1;
%                 gruen=1-S_calc_sortAIC{1,i}.parameters.MS.MAP.par(parindex,1);
%                 blau =1-S_calc_sortAIC{1,i}.parameters.MS.MAP.par(parindex,1);
%                 set(h, 'MarkerFaceColor', [rot gruen blau]);
%             end
%         end    
% 
%     end
%     hold on;
%     set(gcf, 'PaperPositionMode', 'auto');
% end
% plot([i+0.5, i+0.5],[0 length(S_calc_sortAIC)+1],'k-');
% hold on
% ylim([0.8 12]);
% xlim([0.5 length(S_calc_sortAIC)+0.5]);
% for k=1:11
%     lab{k}=D{1,k}.name;
% end
% set(gca,'yTick',[1:1:11]);
% set(gca,'yticklabel',lab);
% set(gca,'XAxisLocation','top')
%  set(gca,'tickdir','out');
% xlabel('calculated models ranked according to AIC (ascending)' );
% ylabel('experiment');
% mycmap=zeros([11 3]);
% mycmap(:,1)=1;
% mycmap(:,2)=1-[0:0.1:1];
% mycmap(:,3)=1-[0:0.1:1];
% 
% 
% 
% colormap(mycmap);
% hcb=colorbar;
% %labels={'0% SAC-population','100% SAC-population'};
% 
% %hcb=lcolorbar(labels)
% set(hcb,'YTick',[1 12],'YTickLabel',{'0','100'}) 
% ylabel(hcb,'% SAC-population') 
%  print('-depsc','./AIC_ranking');
% 
% %% Plot matrix of p values for all compareable calculated models
% 
% 
% 
% 
% 
% for i=1: length(S_calc)
%     %Modellkonfuguration modell i
%     exp_matrix_fullm=[1 0; 1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1];
%     exp_matrix_deleted_i=[0 0; 0 0;0 0;0 0;0 0;0 0;0 0;0 0; 0 0;0  0;0 0];
%     
%     idx = sub2ind(size(exp_matrix_deleted), R_calc{1,i}(:,1), R_calc{1,i}(:,2));
%     exp_matrix_deleted_i(idx)=1;
%     exp_matrix_modeli=exp_matrix_fullm-exp_matrix_deleted_i;
%     
%     for j=1:length(S_calc)
%             exp_matrix_fullm=[1 0; 1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1;1 1];
%             exp_matrix_deleted_j=[0 0; 0 0;0 0;0 0;0 0;0 0;0 0;0 0; 0 0;0  0;0 0];
%             
%             idx = sub2ind(size(exp_matrix_deleted), R_calc{1,j}(:,1), R_calc{1,j}(:,2));
%             exp_matrix_deleted_j(idx)=1;
%             exp_matrix_modelj=exp_matrix_fullm-exp_matrix_deleted_j;
%             
%             % ist j mit i vergleichbar? zeilenweise durchgehen
%             for m =1:size(exp_matrix_fullm,1)
%             flag(m)=(( exp_matrix_modeli(m,1) & exp_matrix_modelj(m,1) & ~exp_matrix_modelj(m,2)) | (exp_matrix_modeli(m,2) & exp_matrix_modelj(m,1) & ~exp_matrix_modelj(m,2)) | (exp_matrix_modelj(m,1) & exp_matrix_modeli(m,1) & exp_matrix_modeli(m,2)) | (exp_matrix_modelj(m,2) & exp_matrix_modeli(m,1) & exp_matrix_modeli(m,2)) | (exp_matrix_modeli(m,2) & exp_matrix_modelj(m,2) & ~exp_matrix_modelj(m,1)));
%             flagsum=sum(flag);
%             end
%     
%            if flagsum==11
%                LR.Lambda(i,j) = exp(S_calc{1,j}.parameters.MS.MAP.logPost - S_calc{1,i}.parameters.MS.MAP.logPost);
%                LR.p(i,j) = 1-chi2cdf(-2*(S_calc{1,j}.parameters.MS.MAP.logPost - S_calc{1,i}.parameters.MS.MAP.logPost),S_calc{1,i}.parameters.number - S_calc{1,j}.parameters.number);
%            else
%                LR.Lambda(i,j) =NaN;
%                LR.p(i,j) = NaN;
%            end 
%           
%     end   
% 
%      
% end
% 
% f3 = figure(3);
% set(f3, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% 
% for i=1:length(S_calc)
%     plot([i-0.5, i-0.5],[0 length(S_calc_sortBIC)+1],'k-');
% hold on;
%  for j=1:length(S_calc)
%      plot([0 length(S_calc_sortBIC)+1],[j-0.5, j-0.5],'k-');
%      hold on;
%      if ~ isnan(LR.p(i,j))
%         % if LR.p(i,j)==0
%                       
%                
%               % plot(i,j,'wo');
%                      
%         % else
%              
%         
%         
%         fill([i-0.5,i-0.5,i+0.5, i+0.5],[j-0.5, j+0.5,j+0.5,j-0.5],[1,1-LR.p(i,j),1-LR.p(i,j)])%,'EdgeColor','none');
%         % end
%      else
%               
%                
%                fill([i-0.5,i-0.5,i+0.5, i+0.5],[j-0.5, j+0.5,j+0.5,j-0.5],[0.8,0.8,0.8])%,'EdgeColor','none');
%               
%               
%               
%      end
%               
% 
%      hold on;
%      
%  end
%  ylim([0.5  length(S_calc)+0.5]);
%  xlim ([0.5  length(S_calc)+0.5]);
% end
%  plot([i+0.5, i+0.5],[0 length(S_calc_sortBIC)+1],'k-');
%   plot([0 length(S_calc_sortBIC)+1],[j+0.5, j+0.5],'k-');
%   
%   set(gca,'tickdir','out');
% xlabel('model index i' );
% ylabel('model index j');
% set(f3, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% colormap(mycmap);
% hcb=colorbar;
% set(hcb,'YTick',[0 1],'YTickLabel',{'0','1'}) 
% ylabel(hcb,'p-value  H_0 = model j') 
%  print('-depsc','./p_values');