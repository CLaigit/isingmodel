function T = saveFrames(temperature)
   % read data and store in a matrix
   filename = strcat('ising',temperature,'.txt');
   lattice_matrix = csvread(filename);
   [nrow, ncol] = size(lattice_matrix);
   % how many figures to plot
   T = nrow/ncol;
   % create directory to store the images
   foldername = strcat('Images/',temperature);
   if ~exist(foldername,'dir')
       mkdir(foldername);
   end
   set(0, 'DefaultFigureVisible', 'off')
   parfor i = 1:T
     toplot = lattice_matrix((i-1)*ncol+1:i*ncol,:) % the lattice to print
     [X, Y] = meshgrid(1:1:ncol);
     figure;
     colormap([1.0 0 0
               0 0 1.0])
     pcolor(X,Y,toplot);
     saveas(gcf,strcat(foldername,'/time_',int2str(i)),'png');
     hold off;
     clf;
   end
