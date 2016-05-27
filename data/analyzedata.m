F = dir('*.txt');

for i = 1:length(F)
    fid = fopen(F(i).name);
    data(:,i) = csvread(F(i).name);
end