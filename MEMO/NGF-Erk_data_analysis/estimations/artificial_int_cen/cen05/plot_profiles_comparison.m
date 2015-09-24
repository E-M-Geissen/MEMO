load('par.mat') 
for k=1:length(parameter_ar)
for i=24:30
parameter_ar{1,k}.name{i}= ['log(',parameter_ar{1,k}.name{i},')']
end
end
options.plot_type='pdf';
options.plot_data.type= 'hist';
options.plot_data.bins=50;
options.hold_on='true';
options.plot_subpop = 'true';
    options.width=18;
    options.height= 10;
plotP_mult(parameter_ar,[], [24:30], [{['cont. data']},{['\Deltax = 10^{0.2}']},{['\Deltax = 10^{0.5}']}],options)
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 options.width options.height])
print('-depsc2','-r1000',['./profiles_comp']);