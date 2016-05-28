clear
F = dir('*.txt');

for i = 1:length(F)
    fid = fopen(F(i).name);
    data(:,i) = csvread(F(i).name);
end

figure
subplot(2,2,1);
scatter(data(1,:), data(3,:));


subplot(2,2,2);
scatter(data(1,:), data(4,:));

subplot(2,2,3);
scatter(data(1,:), data(5,:));

subplot(2,2,4);
scatter(data(1,:), data(6,:));