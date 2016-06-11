%% Save all the data in one matrix
clear
data_dirs = dir('.');
count = 0;
folder = {};
for i = 3:length(data_dirs)
    if data_dirs(i).isdir == 1 && ~isempty(strfind(data_dirs(i).name,'lattice'))
        count = count + 1;
        folder{count} = data_dirs(i).name;
        folder_name = char(folder(count));
        F = dir(strcat(folder_name, '/*.txt'));
        num_files = length(F);
        a = length(csvread(strcat(folder_name, '/', F(1).name)));
        for j = 1:num_files
            data(((count-1)*a+1):count*a,j) = csvread(strcat(folder_name, '/', F(j).name));
        end
    end
end

%% calculate the steepest position and get the critical temperature
for i = 1:count
    lattice_length(i,:) = data(2 + 6*(i-1));
    data_lattice = data(3 + (i-1) * 6: 6 + (i-1) * 6,:);
    num_temp = size(data,2);
    for j = 1:num_temp-1
        temp_diff = data(1,j+1) - data(1,j);
        matrix_slope(1+(i-1)*4:4 +(i-1)*4,j) = (data_lattice(:, j+1) - data_lattice(:, j))/temp_diff;
    end
end
[slope, index] = max(abs(matrix_slope),[],2);
temp = data(1,:);
for i = 1:length(index)
    critical_temp(i,:) = temp(index(i));
end

%% Plot four subplots, each one is the comparision of different size of lattice
set(0,'DefaultFigureVisible','off');
foldername = 'plots';
if ~exist(foldername,'dir')
    mkdir(foldername);
end
figure();
for i = 1:count
    x = data(1 + a*(i-1),:);
    y = data(3 + a*(i-1),:);
    plotdata(x,y);
    hold on;
end
legend(folder);
xlabel('kT/J');
ylabel('Average Energy');
title('Temperature vs. Average Energy');
hold off;
saveas(gcf, strcat(foldername, '/Average Energy'),'png');
clf;

figure();
for i = 1:count
    x = data(1 + a*(i-1),:);
    y = data(4 + a*(i-1),:);
    [temp_filter, energy_filter] = plotdata(x,y);
    [max_energy, index] = max(energy_filter);
    critical_temp(2 + 4*(i-1),:) = temp_filter(index);
    hold on;
end
legend(folder);
xlabel('kT/J');
ylabel('Heat Compacity');
hold off;
title('Temperature vs. Heat Compacity');
saveas(gcf, strcat(foldername, '/Heat Compacity'),'png');
clf;

figure();
for i = 1:count
    x = data(1 + a*(i-1),:);
    y = data(5 + a*(i-1),:);
    plotdata(x,y);
    hold on;
end
legend(folder);
xlabel('kT/J');
ylabel('Average Magnetization');
title('Temperature vs. Average Magnetization');
hold off;
saveas(gcf, strcat(foldername, '/Average Magnetization'),'png');
clf;

figure();
for i = 1:count
    x = data(1 + a*(i-1),:);
    y = data(6 + a*(i-1),:);
    [temp_filter, energy_filter] = plotdata(x,y);
    [max_energy, index] = max(energy_filter);
    critical_temp(4 + 4*(i-1),:) = temp_filter(index);
    hold on;
end
legend(folder);
xlabel('kT/J');
ylabel('Magnetization Susceptibility');
title('Temperature vs. Magnetization Susceptibility');
hold off;
saveas(gcf, strcat(foldername, '/Magnetization Susceptibility'),'png');
clf;

finalresult = [repmat(lattice_length, 4, 1), critical_temp];
finalresult_sort = sortrows(finalresult);
%% compare the average magnetization with theoretical magnetization
beta = 1./temp;
J = 1;
theoretical_magnetization = (1 - sinh(2*beta*J).^(-4)).^(1/8);
figure();
for i = 1:count
    x = data(1 + a*(i-1),:,:);
    y = data(5 + a*(i-1),:);
    plotdata(x,y);
    hold on;
end
plotdata(temp(1:25), theoretical_magnetization(1:25));
legend(folder, 'theoretical');
xlabel('kT/J');
ylabel('Magnetization');
title('Theoretical Magnetization vs. Average Magnetization');
hold off;
saveas(gcf, strcat(foldername, '/Magnetization Comparison'),'png');
clf;