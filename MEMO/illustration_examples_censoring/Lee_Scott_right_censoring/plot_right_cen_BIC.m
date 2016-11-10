clear all
close all
load all
clear x
clear y
mus=1;
%% set figure size
width=5
height=5

TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize =6;
set(0,TextSizes);

%% VISUALIZATION OF RESULTS


% BIC:
%figure('name','BIC');
for i = 1:2
    
    switch i
        case 1, perc = N_correct_BIC_truncated;  ts = 'right truncation';name='v1'
        case 2, perc = N_correct_BIC_censored; ts = 'right censoring'; name='v2';
    end
    
    hFig(i)=figure;
    set(hFig(i), 'Units','centimeters','Position', [10 10 width height])
    
    t=1;
    for p=1:length(sigma)
        for j=1:length(cen_time)
            x(t)=cen_time(j);
            y(t)=sqrt((exp(sigma(p)^2)-1)*exp(2*mus+sigma(p)^2));
            %y(t)=sigma(p)^2;
            z(t)=perc(p,j);
            t=t+1;
        end
    end
    
    
    
    s=zeros(size(x));
    s(:,:)=20;
    
    
    colormap gray
    %colormap(flipud(colormap))
    scatter(x, y,s,z,'fill');
    hold on
    %cc=zeros(length(x),3);
%     c(1,:)=[0 0.749019622802734 0.749019622802734]; 
%     for cm=2:5
%       
%       c(cm-(cm-1)+((cm-1)*20+cm-1),:)=[0 0.749019622802734 0.749019622802734];  
%     end
   
    %h=scatter(x, y,s,cc,'LineWidth',1);
    h=scatter(x, y,s,[0,0,0],'LineWidth',1);
    hold all
    plot(cen_time(1:1:5),sqrt((exp(sigma(1:2:9).^2)-1).*exp(2*mus+sigma(1:2:9).^2)),'-','Color',[0 0.749019622802734 0.749019622802734],'LineWidth',1);
   set(gca, 'LooseInset', [0,0,0,0]);
    hold all
     ylim([-0.5 21.5])
     xlim([2.2 14.5])
    axis square
        %pos=plotboxpos(gca)
    %pos(3)=pos(4)
     %set(gca,'Position',pos)
    %colorbar
    ylabel('standard deviation');
    xlabel('end of measurement');
    caxis([0 100]);
    %axis equal
   
    hc = colorbar;
    ylabel(hc,'% identified correctly')
    % xlim([2.5 8.5])
    

    
    
    
    title(ts);
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
    
    %set(hFig(i),'renderer','opengl'); % the magic line 
%saveas(gcf,['./figs/right_cen_',name],'epsc2');
    print('-depsc2','-r1000',['./figs/BIC_right_cen_',name]);
  
    %fix_pcolor_eps(['./figs/right_cen_',name,'.eps'])
end
