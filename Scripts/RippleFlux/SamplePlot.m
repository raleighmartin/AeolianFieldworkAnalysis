%initialize script
clearvars; %clears all variables
close all; %closes all figure windows

%creating our synthetic timeseries
t = 1:1000;
z1 = sin(t/10)+rand(1,1000); %sine wave with random noise
z2 = sin(t/10-20)+rand(1,1000); %another sine wave with random noise

%find zero points of z1
ind_z1_0 = find(round(z1,1)==0); %indices of points that meet criteria of z1=0 (rounded to nearest 10th)
t_z1_0 = t(ind_z1_0); %times corresponding to these indices
z_z1_0 = z1(ind_z1_0); %z1's corresponding to these indices

%find the point where z1 = its maximum for each 100 points
wavelength = round(2*pi*10);
z1max_wavelength = []; %initialize list of z1max's
t_z1max_wavelength = []; %initialize list of corresponding times

for i = 1:(round(1000/wavelength)-1)
   ind_range = (i-1)*wavelength+(1:wavelength);
   t_range = t(ind_range);
   z1_range = z1(ind_range);
   z1max_range = max(z1_range);
   ind_z1max_range = find(z1_range==z1max_range);
   t_z1max_range = t_range(ind_z1max_range);
   
   %add these values to lists
   z1max_wavelength = [z1max_wavelength z1max_range];
   t_z1max_wavelength = [t_z1max_wavelength t_z1max_range];
end

%create plot
figure; %creates a new figure
hold on; %prevents the figure from erasing previous lines
plot(t,z1);
plot(t,z2);
plot(t_z1_0, z_z1_0, 'o','MarkerSize',10);
plot(t_z1max_wavelength, z1max_wavelength, 'x','MarkerSize',20);

ylim([-1 3]);
xlim([0 1000]);

%plot lines corresponding to zero crossings
ylimits = ylim;
for i = 1:length(t_z1_0)
    plot([t_z1_0(i) t_z1_0(i)],ylimits,'k')
end

legend('z1','z2','z1 zero points','z1 maxima','Location','NorthWest');
xlabel('time (s)');
ylabel('elevation (mm)');
title('Synthetic data');
set(gca,'FontSize',14);

%print plot
print('sampleplot.png','-dpng');