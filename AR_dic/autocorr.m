clear all; close all; clc
%set dic images
%===========================

%% images


imagepath = strcat('../AR')
imagefiles = dir(fullfile(imagepath,'*.jpg'));% .bmp or .tif 
imagefiles = {imagefiles.name}';
images = strcat(imagepath,filesep,imagefiles);

%%
% Computes autocorrelation of a 2D matrix.
% 
% The autocorrelation si the correlation of a signal with a delayed copy of 
% itself.
% This gives us a measure of "quality" of our speckle:
% A good speckle should give us an autocorrelation length of ~5 pixels

filename = images{1};

% initial fiugre
img = imread(filename);
img = rgb2gray(img);
[nbx,nby]=size(img);

x1 = 410;
y1 = 1;
x2 = 645;
y2 = 1700;

crop.x1=x1;
crop.x2=x2;
crop.y1=y1;
crop.y2=y2;

fig = figure
ttl=split(pwd,'/');
ttl=strrep(ttl{end},'_',' ');

subplot(2,2,[1,2])
title({ttl,'reference figure with ROI'})
hold on
imagesc(img)
plot([y1,y2,y2,y1,y1],[x1,x1,x2,x2,x1])

daspect([1,1,1])
colormap(gca,'gray')
xlim([0, nby])
ylim([0, nbx])
hold off
img = img(x1:x2,y1:y2);
caxis([0,2^8])

% autocorrelation
subplot(2,2,3)
hold on
title('autocorrelation')
[ny,nx]=size(img);
x=linspace(0,nx,nx)-nx/2;
y=linspace(0,ny,ny)-ny/2;
imagesc('XData',x,'YData',y,'CData',autocorrr(img))
colormap(gca,'Parula')
xlim([-20, 20])
ylim([-20, 20])
hold off

% histogram
subplot(2,2,4)
hold on
title('histogram of graylevels')
imgflat = reshape(img,1,[]);
hist(single(imgflat),1000)%2^12)
xlim([0,2^8]);
hold off

min(imgflat)
max(imgflat)

saveas(fig,'DIC_images_ref.png');

save('DIC_images.mat','crop','images')


function [R]=autocorrr(im)
Fr = fft2(im);
S = Fr.*conj(Fr);
R = ifft2(S);
R = fftshift(R);
R = real(R);

end