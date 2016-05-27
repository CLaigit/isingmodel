clear

data = csvread('isingT4.5.dat');
lattice_length = 20;

[column, row] = size(data);

numberOfmatrix = column/row;

one = zeros(row,row,numberOfmatrix);
for i = 1:numberOfmatrix
    one(:,:,i) = data(1+(i-1)*lattice_length:i*lattice_length, 1:row);
end

E = 0;
for i = 1:numberOfmatrix
    E = E + energy(one(:,:,i));
end
E = E / numberOfmatrix;

% 
% k = 8007
% hold on
% for i = 1:lattice_length
%     for j = 1:lattice_length
%         if one(i,j,k) == 1
%             scatter(i,j,'b',...
%                     'MarkerFaceColor',[0.7 0 .7],...
%                     'LineWidth',10)
%         else
%             scatter(i,j,'y',...
%                     'MarkerFaceColor',[0 .7 .7],...
%                     'LineWidth',10)
%         end
%     end
% 
% end
