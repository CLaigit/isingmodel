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

%% Plot four subplots, each one is the comparision of different size of lattice
set(0,'DefaultFigureVisible','off');
foldername = 'plots';
if ~exist(foldername,'dir')
    mkdir(foldername);
end
figure();
for i = 1:count
    scatter(data(1,:), data(3 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel('Average Energy');
hold off;
saveas(gcf, strcat(foldername, '/Average Energy'),'png');
clf;

figure();
for i = 1:count
    scatter(data(1,:), data(4 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel('Heat Compacity');
hold off;
saveas(gcf, strcat(foldername, '/Heat Compacity'),'png');
clf;

figure();
for i = 1:count
    scatter(data(1,:), data(5 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel('Average Magnetization');
hold off;
saveas(gcf, strcat(foldername, '/Average Magnetization'),'png');
clf;

figure();
for i = 1:count
    scatter(data(1,:), data(6 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel('Magnetization Susceptibility');
hold off;
saveas(gcf, strcat(foldername, '/Magnetization Susceptibility'),'png');
clf;

%% calculate the steepest position and get the critical temperature
data_sorted = sort(data,2);
for i = 1:count
    lattice_length(i,:) = data_sorted(2 + 6*(i-1));
    data_lattice = data_sorted(3 + (i-1) * 6: 6 + (i-1) * 6,:);
    num_temp = size(data_sorted,2);
    for j = 1:num_temp-1
        temp_diff = data_sorted(1,j+1) - data_sorted(1,j);
        matrix_slope(1+(i-1)*4:4 +(i-1)*4,j) = (data_lattice(:, j+1) - data_lattice(:, j))/temp_diff;
    end
end
[slope, index] = max(matrix_slope,[],2);
temp = data_sorted(1,:);
for i = 1:length(index)
    critical_temp(i,:) = temp(index(i));
end
finalresult = [repmat(lattice_length, count,1),slope , critical_temp];