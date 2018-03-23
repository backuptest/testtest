
for i = 1:100


img1 = Fg_cluster(:,i);
img1 = reshape(img1, imSize);

img2 = Fg_repro(:,i);
img2 = reshape(img2,imSize);
%img = uint8(img);
set(gca, 'position', [0 0 1 1]);
subplot(1,2,1);
imshow(img1,[]);
subplot(1,2,2);
imshow(img2,[]);
N(i)=getframe;
end


