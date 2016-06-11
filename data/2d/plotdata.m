function [temp_plot, energy_plot] = plotdata(temp, energy)
     tosort = [temp', energy'];
     sorted = sortrows(tosort);
     temp_sorted = sorted(:,1);
     energy_sorted = sorted(:,2);
     temp_plot = temp_sorted;
     energy_plot = energy_sorted;
     plot(temp_plot, energy_plot);
