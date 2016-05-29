function produceVideo(temperature, T)
   framepath = strcat('Images/',temperature);
   path_save = strcat('Videos/', temperature);
   v = VideoWriter(path_save);
   v.FrameRate = 15;
   open(v);
   for i = 1:T
       frame = imread(strcat(framepath, '/time_', int2str(i), '.png'));
       writeVideo(v, frame);
   end