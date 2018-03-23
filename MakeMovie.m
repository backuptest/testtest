writerObj = VideoWriter('StrucSparsityCompare1.avi'); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj); 
figure;

for i=1:100    
    img1 = Fg_cluster(:,i);
img1 = reshape(img1, imSize);

img2 = Fg_repro(:,i);
img2 = reshape(img2,imSize);
%img = uint8(img);
set(gca, 'position', [0 0 1 1]);
subplot(1,2,1);
imshow(img1,[]);
title('graph-sparsity')
subplot(1,2,2);
imshow(img2,[]);
title('ReProCS')

    
   
    
    
 
    

        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
  
 
end
hold off
close(writerObj); % Saves the movie.
