%%%%%%%%%%%%%%
% INITIALIZE %
%%%%%%%%%%%%%%

clearvars;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFORMATION FOR SAVING PLOT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder_Plots = '../../PlotOutput/GrainSize/'; %folder for plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFORMATION FOR DETERMINING WHERE TO PLACE FACTORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = linspace(0.2,2,100);
% y_logincrease = 0.5*log(x)+0.2554;
% y_logdecrease = -y_logincrease;

x = unique([logspace(log10(0.2),log10(0.4),7),logspace(log10(0.4),log10(0.8),7),logspace(log10(0.8),log10(2),10)]);
y_logincrease = unique([linspace(-0.5493,-0.2027,7),linspace(-0.2027,0.1438,7),linspace(0.1438,0.6020,10)]);
y_logdecrease = fliplr(unique([linspace(0.5493,0.2027,7),linspace(0.2027,-0.1438,7),linspace(-0.1438,-0.6020,10)]));

% y_logincrease = linspace(-0.5,0.5,50);
% y_logdecrease = linspace(0.5,-0.5,50);

ind_fine = find(x<=0.4); N_fine = length(ind_fine);
ind_medium = find(x>= 0.4 & x<=0.8); N_medium = length(ind_medium);
ind_coarse = find(x>=0.8); N_coarse = length(ind_coarse);

% y_linearincrease = linspace(-0.2,0.4,100);
% y_lineardecrease = -y_linearincrease;
% y_decaydecrease = 0.0002*(1./(x-0.07)).^4;
% y_decayincrease = -0.01*(1./(x+0.1)).^4;
% y_parabolic = 1-linspace(-1,1,20).^2;

%% factors contributing to promotion of saltation
% %ejection speed
% ind_ejectionspeed = unique([ind_fine, ind_medium, ind_coarse]);
% y_ejectionspeed = y_logdecrease/2;
% y_ejectionspeed(ind_fine) = y_ejectionspeed(ind_medium(1));

ind_ejectionspeed = unique([ind_fine, ind_medium, ind_coarse]);
y_fine = y_logdecrease(ind_fine)*0;
y_medium = y_logdecrease(ind_medium)/4;
dy_finemedium = y_medium(1) - y_fine(end);
y_medium = y_medium-dy_finemedium;
y_coarse = y_logdecrease(ind_coarse)/2;
dy_mediumcoarse = y_coarse(1) - y_medium(end);
y_coarse = y_coarse-dy_mediumcoarse;
y_ejectionspeed = zeros(size(x));
y_ejectionspeed(ind_fine) = y_fine;
y_ejectionspeed(ind_medium) = y_medium;
y_ejectionspeed(ind_coarse) = y_coarse;
y_ejectionspeed = y_ejectionspeed+0.05;

%horizontal acceleration by wind
ind_windaccel = unique([ind_fine, ind_medium, ind_coarse]);
y_fine = y_logdecrease(ind_fine)/10;
y_medium = y_logdecrease(ind_medium)/3;
dy_finemedium = y_medium(1) - y_fine(end);
y_medium = y_medium-dy_finemedium;
y_coarse = y_logdecrease(ind_coarse)/1.5;
dy_mediumcoarse = y_coarse(1) - y_medium(end);
y_coarse = y_coarse-dy_mediumcoarse;
y_windaccel = zeros(size(x));
y_windaccel(ind_fine) = y_fine;
y_windaccel(ind_medium) = y_medium;
y_windaccel(ind_coarse) = y_coarse;
y_windaccel = y_windaccel+0.1;
 
%potential for rebound
ind_restitution = unique([ind_fine, ind_medium, ind_coarse]);
y_restitution = y_logdecrease/10;

%% factors contributing to travel length
ind_bedcapture = unique([ind_fine, ind_medium, ind_coarse]);
y_bedcapture = y_logdecrease/6;

%% additional factors contributing to saltation height

% %vertical drag inhibition of trajectories
% ind_verticaldrag = unique([ind_fine, ind_medium]);
% y_verticaldrag = zeros(size(x));
% y_verticaldrag(ind_fine) = y_logincrease(ind_fine) - y_logincrease(ind_medium(end));
% y_verticaldrag(ind_medium) = y_logincrease(ind_medium) - y_logincrease(ind_medium(end));

%vertical drag inhibition of trajectories
ind_verticaldrag = unique([ind_fine, ind_medium, ind_coarse]);
y_fine = y_logincrease(ind_fine)*1.5;
y_medium = y_logincrease(ind_medium)*1.18;
dy_finemedium = y_medium(1) - y_fine(end);
y_medium = y_medium-dy_finemedium;
y_coarse = y_logincrease(ind_coarse)/4;
dy_mediumcoarse = y_coarse(1) - y_medium(end);
y_coarse = y_coarse-dy_mediumcoarse;
y_verticaldrag = zeros(size(x));
y_verticaldrag(ind_fine) = y_fine;
y_verticaldrag(ind_medium) = y_medium;
y_verticaldrag(ind_coarse) = y_coarse;
y_verticaldrag = y_verticaldrag-y_verticaldrag(ind_coarse(1));
y_verticaldrag = y_verticaldrag + 0.2;

%suspension enhancement of fine trajectores
ind_turbulence = ind_fine;
y_turbulence = zeros(size(x));
y_turbulence(ind_fine) = (y_logdecrease(ind_fine) - y_logdecrease(ind_fine(end)))*0.9;

%response time
%ind_responsetime = unique([ind_fine, ind_medium, ind_coarse]);
ind_responsetime = unique([ind_coarse]);
y_responsetime = zeros(size(x));
y_responsetime(ind_coarse) = y_logdecrease(ind_coarse)/2;
y_responsetime(ind_fine) = y_responsetime(ind_coarse(1));
y_responsetime(ind_medium) = y_responsetime(ind_coarse(1));
y_responsetime = y_responsetime-y_responsetime(1);

%% factors contributing to size fraction

% %ejection rate
% ind_ejectionrate = unique([ind_fine, ind_medium, ind_coarse]);
% y_ejectionrate = y_logincrease/4;

%ejection rate
ind_ejectionrate = unique([ind_coarse]);
y_fine = y_logincrease(ind_fine)*0;
y_medium = y_logincrease(ind_medium)*0;
dy_finemedium = y_medium(1) - y_fine(end);
y_medium = y_medium-dy_finemedium;
y_coarse = y_logincrease(ind_coarse)/4;
dy_mediumcoarse = y_coarse(1) - y_medium(end);
y_coarse = y_coarse-dy_mediumcoarse;
y_ejectionrate = zeros(size(x));
y_ejectionrate(ind_fine) = y_fine;
y_ejectionrate(ind_medium) = y_medium;
y_ejectionrate(ind_coarse) = y_coarse;
y_ejectionrate = y_ejectionrate-y_ejectionrate(ind_coarse(1));

%promotion to saltation
ind_saltationpromotion = unique([ind_fine, ind_medium, ind_coarse]);
y_saltationpromotion = y_ejectionspeed+y_windaccel+y_restitution;

%saltation height
ind_height = unique([ind_fine, ind_medium, ind_coarse]);
y_height = y_verticaldrag+y_windaccel+y_restitution+y_turbulence;

%trajectory length
ind_length = unique([ind_fine, ind_medium, ind_coarse]);
y_length = y_height+y_restitution;

%% compute size fraction
ind_sizefraction = unique([ind_fine, ind_medium, ind_coarse]);
y_sizefraction = y_ejectionrate + y_saltationpromotion + y_length;

%% information for plotting
PlotFont = 12; %font for labels
color_ejectionspeed = [0    0.4470    0.7410];
marker_ejectionspeed = '+';
color_ejectionrate = [0    0.4470    0.7410];
marker_ejectionrate = 'o';
color_windaccel = [0.7    0.6940    0.1250];
marker_windaccel = '*';
color_responsetime = [0.7    0.6940    0.1250];
marker_responsetime = 'x';
color_restitution = [0.5250    0.4496    0.4036];
marker_restitution = 's';
color_bedcapture = [0.9290    0.6940    0.1250];
marker_bedcapture = 'diamond';
color_height = [0.4940    0.1840    0.5560];
marker_height = '^';
color_verticaldrag = [0.6350    0.0780    0.1840];
marker_verticaldrag = 'v';
color_turbulence = [0.3010    0.7450    0.9330];
marker_turbulence = 'pentagram';
color_saltationpromotion = [0.4660    0.6740    0.1880];
marker_saltationpromotion = '>';
color_length = [0.8500    0.3250    0.0980];
marker_length = '<';
color_sizefraction = [0 0 0];
marker_sizefraction = 'hexagram';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFORMATION ABOUT SUBPLOT POSITIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot_positions = {...
    [0.15 0.7 0.82 0.29],...
    [0.15 0.37 0.31 0.27],...
    [0.64 0.37 0.31 0.27],...
    [0.15 0.06 0.82 0.24]};

%%%%%%%%%%%%%%%%%%%%%
% GENERATE SUBPLOTS %
%%%%%%%%%%%%%%%%%%%%%

%initialize plot and subplots
figure(1);
h_subplot = gobjects(4,1);

for i = 1:4
    h_subplot(i) = subplot('position',subplot_positions{i}); hold on;
    plot([0.2 2],[0 0],'k','LineWidth',1);
    plot([0.4 0.4],[-1 1],'k--','LineWidth',1);
    plot([0.8 0.8],[-1 1],'k--','LineWidth',1);
    %if i==1 || i==4
    if i<=4
        text(0.28,0.88,'Fine','HorizontalAlignment','Center','FontSize',PlotFont);
        text(0.56,0.88,'Medium','HorizontalAlignment','Center','FontSize',PlotFont);
        text(1.05,0.88,'Coarse','HorizontalAlignment','Center','FontSize',PlotFont);
    else
        text(0.28,-0.88,'Fine','HorizontalAlignment','Center','FontSize',PlotFont);
        text(0.56,-0.88,'Medium','HorizontalAlignment','Center','FontSize',PlotFont);
        text(1.05,-0.88,'Coarse','HorizontalAlignment','Center','FontSize',PlotFont);
    end
    xlim([0.2 1.5]);
    ylim([-1 1]);
    set(gca,'XScale','log','YTick',[],'Box','On','FontSize',PlotFont);
    set(gca,'xtick',[0.2:0.1:1, 1.1:0.1:1.5]);
    set(gca,'xticklabels',{'0.2','','0.4','','','','0.8','','','','','','','1.5'});

    %subplot for size fraction
    if i == 1
        scalingfactor = 2;
        h1_1 = plot(x(ind_ejectionrate),y_ejectionrate(ind_ejectionrate)*scalingfactor,'Marker',marker_ejectionrate,'MarkerSize',5,'LineWidth',1.0,'Color',color_ejectionrate); %ejection rate
        h1_2 = plot(x(ind_saltationpromotion),y_saltationpromotion(ind_saltationpromotion)*scalingfactor,'Marker',marker_saltationpromotion,'MarkerSize',5,'LineWidth',1.0,'Color',color_saltationpromotion); %promotion to saltation
        h1_3 = plot(x(ind_length),y_length(ind_length)*scalingfactor,'Marker',marker_length,'MarkerSize',5,'LineWidth',1.0,'Color',color_length); %trajectory length
        h1_4 = plot(x(ind_sizefraction),y_sizefraction(ind_sizefraction)*scalingfactor,'LineWidth',2.0,'Color',color_sizefraction); %combined - size fraction
        %h1_4 = plot(x(ind_sizefraction),y_sizefraction(ind_sizefraction)*scalingfactor,'Marker',marker_sizefraction,'MarkerSize',4,'LineWidth',2.0,'Color',color_sizefraction); %combined - size fraction
                
        h1_legend = legend([h1_1 h1_2 h1_3 h1_4],...
            'Ejection rate','Promotion to saltation','Travel length','Combined effect');
        
        set(h1_legend,'Location','EastOutside','FontSize',10);       
        text(0.175,0.9,'(a)','FontSize',PlotFont);
%         ylabel({'Enhance';'';'          $$\Uparrow$$';'';...
%             '\textbf{Effect on}';'\textbf{size-sel.}';'\textbf{mobility}';...
%             '';'          $$\Downarrow$$';'';'Reduce'},...
%             'Interpreter','Latex');
        ylabel({'Enhance';'          $$\Uparrow$$';...
            '\textbf{Effect on}';'\textbf{size-sel.}';'\textbf{mobility}';...
            '          $$\Downarrow$$';'Reduce'},...
            'Interpreter','Latex','FontSize',PlotFont);
        set(get(gca,'YLabel'),'Rotation',0,'VerticalAlignment','Middle');
        set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position')+[-0.05 0 0])
    
    elseif i == 2
        scalingfactor = 2;
        h2_1 = plot(x(ind_ejectionspeed),y_ejectionspeed(ind_ejectionspeed)*scalingfactor,'Marker',marker_ejectionspeed,'MarkerSize',5,'LineWidth',1.0,'Color',color_ejectionspeed); %ejection speed
        h2_2 = plot(x(ind_windaccel),y_windaccel(ind_windaccel)*scalingfactor,'Marker',marker_windaccel,'MarkerSize',5,'LineWidth',1.0,'Color',color_windaccel); %wind-driven accel.       
        h2_3 = plot(x(ind_restitution),y_restitution(ind_restitution)*scalingfactor,['-.',marker_restitution],'MarkerSize',5,'LineWidth',1.0,'Color',color_restitution); %restitution
        h2_4 = plot(x(ind_saltationpromotion),y_saltationpromotion(ind_saltationpromotion)*scalingfactor,'LineWidth',2.0,'Color',color_saltationpromotion); %saltation promotion
        %h2_4 = plot(x(ind_saltationpromotion),y_saltationpromotion(ind_saltationpromotion)*scalingfactor,'Marker',marker_saltationpromotion,'MarkerSize',4,'LineWidth',2.0,'Color',color_saltationpromotion); %saltation promotion
        
        h2_legend = legend([h2_1 h2_2 h2_3 h2_4],...
            'Ejection speed','Wind-driven accel.','Restitution coef.','Combined effect');
        
        set(h2_legend,'Location','SouthWest','FontSize',10);       
        text(0.16,0.8,'(b)','FontSize',PlotFont);
        ylabel({'Enhance';'          $$\Uparrow$$';...
            
            '\textbf{Effect on}';'\textbf{promotion}';'\textbf{to saltation}';...
            '          $$\Downarrow$$';'Reduce'},...
            'Interpreter','Latex','FontSize',PlotFont);
        set(get(gca,'YLabel'),'Rotation',0,'VerticalAlignment','Middle');
        set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position')+[-0.07 0 0])
    
    elseif i == 3
        scalingfactor = 2;
        h3_1 = plot(x(ind_height),y_height(ind_height)*scalingfactor,'Marker',marker_height,'MarkerSize',5,'LineWidth',1.0,'Color',color_height);
        h3_2 = plot(x(ind_bedcapture),y_bedcapture(ind_bedcapture)*scalingfactor,['-.',marker_bedcapture],'MarkerSize',5,'LineWidth',1.0,'Color',color_bedcapture);
        h3_3 = plot(x(ind_length),y_length(ind_length)*scalingfactor,'LineWidth',2.0,'Color',color_length);
        %h3_3 = plot(x(ind_length),y_length(ind_length)*scalingfactor,'Marker',marker_length,'MarkerSize',4,'LineWidth',2.0,'Color',color_length);
        
        h3_legend = legend([h3_1 h3_2 h3_3],...
            'Saltation height','Bed capture prob.','Combined effect');
    
        set(h3_legend,'Location','SouthEast','FontSize',10);
        text(0.16,0.8,'(c)','FontSize',PlotFont);
        ylabel({'Enhance';'          $$\Uparrow$$';...
            '\textbf{Effect on}';'\textbf{travel}';'\textbf{length}';...
            '          $$\Downarrow$$';'Reduce'},...
            'Interpreter','Latex','FontSize',PlotFont);
        set(get(gca,'YLabel'),'Rotation',0,'VerticalAlignment','Middle');
        set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position')+[-0.07 0 0])
   
    elseif i == 4
        scalingfactor = 1.5;
        h4_1 = plot(x(ind_verticaldrag),y_verticaldrag(ind_verticaldrag)*scalingfactor,'Marker',marker_verticaldrag,'MarkerSize',5,'LineWidth',1.0,'Color',color_verticaldrag);
        h4_2 = plot(x(ind_restitution),y_restitution(ind_restitution)*scalingfactor,['-.',marker_restitution],'MarkerSize',5,'LineWidth',1.0,'Color',color_restitution);
        h4_3 = plot(x(ind_responsetime),y_responsetime(ind_responsetime)*scalingfactor,'Marker',marker_responsetime,'MarkerSize',5,'LineWidth',1.0,'Color',color_responsetime);
        h4_4 = plot(x(ind_turbulence),y_turbulence(ind_turbulence)*scalingfactor,'Marker',marker_turbulence,'MarkerSize',5,'LineWidth',1.0,'Color',color_turbulence);
        h4_5 = plot(x(ind_height),y_height(ind_height),'LineWidth',2.0,'Color',color_height);
        %h4_5 = plot(x(ind_height),y_height(ind_height),'Marker',marker_height,'MarkerSize',4,'LineWidth',2.0,'Color',color_height);

        h4_legend = legend([h4_1 h4_3 h4_2 h4_4 h4_5],...
            'Vertical drag','Response time','Restitution coeff.','Turbulence effects','Combined effect');
        
        set(h4_legend,'Location','EastOutside','FontSize',10);
        text(0.175,0.9,'(d)','FontSize',PlotFont);
%         ylabel({'Enhance';'';'          $$\Uparrow$$';'';...
%             '\textbf{Effect on}';'\textbf{saltation}';'\textbf{height}';...
%             '';'          $$\Downarrow$$';'';'Reduce'},...
%             'Interpreter','Latex','FontSize',PlotFont);
        ylabel({'Enhance';'          $$\Uparrow$$';...
            '\textbf{Effect on}';'\textbf{saltation}';'\textbf{height}';...
            '          $$\Downarrow$$';'Reduce'},...
            'Interpreter','Latex','FontSize',PlotFont);
        set(get(gca,'YLabel'),'Rotation',0,'VerticalAlignment','Middle');
        set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position')+[-0.05 0 0])
    end   

    xlabel('Normalized grain size, $$d_i / d_{50,bed}$$','Interpreter','Latex','FontSize',PlotFont);
end

set(gcf,'PaperUnits','inches','PaperSize',[7 9],'PaperPosition',[0 0 7 9],'PaperPositionMode','Manual');
print([folder_Plots,'ConceptualFigure.png'],'-dpng');
