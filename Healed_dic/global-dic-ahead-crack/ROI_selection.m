% dic_clip_position
%
%
% ROI coordinate
clear all; close all; clc


set(0,'defaultTextInterpreter','none'); %trying to set the default
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',10);

imagepath = '../../Healed';
imagefiles = dir(fullfile(imagepath,'*.jpg'));% .bmp or .tif //uncompressed
imagefiles = {imagefiles.name}';
images = strcat(imagepath,filesep,imagefiles);
i=23
f = imread(images{i});%-imread(images{i-1});;
f = rgb2gray(f);

save_file = ['ROI_selection.mat'];

try
load(save_file);
end
if ~exist('P')
    figure
    hold on
    title({'p1,p2,p3 = notch tip, notch top, notch bot',' p4,p5,p6 = right side top, right side center, right side bot'})
    imshow(f)
    hold off
    P = ginput(6)
end

%%

[top,bot,global_roi,global_roi_ahead,info] = compute_dic_clip_position(P);

fig = figure

imshow(f)
hold on
sl =[1,3,2,4,1];
plot(top.p(sl,1),top.p(sl,2),'-')
plot(bot.p(sl,1),bot.p(sl,2),'-')
plot(global_roi.p(sl,1),global_roi.p(sl,2),'-')
plot(global_roi_ahead.p(sl,1),global_roi_ahead.p(sl,2),'-')
sl =[1,2];
plot(top.p(sl,1),top.p(sl,2),'ko')
plot(bot.p(sl,1),bot.p(sl,2),'ko')
plot(global_roi.p(sl,1),global_roi.p(sl,2),'ko')

plot(global_roi_ahead.p(sl,1),global_roi_ahead.p(sl,2),'ko')

plot(P(:,1),P(:,2),'go')

hold off

saveas(fig,'ROI_selection.png')
save(save_file,'P','global_roi','global_roi_ahead','top','bot','info')
%save([wdir,'global/roi.mat'],'global_roi')
%save([wdir,'COD/top/roi.mat'],'top')
%save([wdir,'COD/bot/roi.mat'],'bot')

function [top,bot,global_roi,global_roi_ahead,info] = compute_dic_clip_position(P)

notch_tip = P(1,:);
notch_top = P(2,:);
notch_bot = P(3,:);
r1 = P(4,:);
r2 = P(5,:);
r3 = P(6,:);

up1 = r2-r1;
up = up1/norm(up1);

up2 = r3-r2;
up2 = up2/norm(up2);

up3 = r3-r1;
up3 = up3/norm(up3);

up = up1+up2+up3;
up = up/norm(up);

side = [-up(2),up(1)]
rightavg = (r1 + r2 + r3)/3;


scale = dot(rightavg-notch_tip,side)/30.0 % pixel/mm 
up = up*scale;
side =side*scale;
%
ROI_x = 5; % mm
ROI_y = 2.5;% mm
dist_up = 1; % (mm) istance from crack
% 2 small ROI for COD
ROI_side = side*ROI_x;
ROI_up = up*ROI_y;

top.p=zeros(4,2);
bot.p=zeros(4,2);

top.p(1,:) = notch_tip  + ROI_up+dist_up*up;
top.p(2,:) = top.p(1,:) - ROI_up+ROI_side;
top.p(3,:) = top.p(1,:) - ROI_up;
top.p(4,:) = top.p(1,:) + ROI_side;

bot.p(1,:) = notch_tip - dist_up*up;
bot.p(2,:) = bot.p(1,:) - ROI_up + ROI_side;
bot.p(3,:) = bot.p(1,:) - ROI_up;
bot.p(4,:) = bot.p(1,:) + ROI_side;

% 1 big ROI for global analysis
%
ROI_x = 30; % mm
ROI_y = 4.4;% mm
dist_up = ROI_y/2; % (mm) istance from crack

ROI_side = side*ROI_x;
ROI_up = up*ROI_y;

global_roi.p=zeros(4,2);
global_roi.p(1,:) = notch_tip + ROI_up-dist_up*up;
global_roi.p(2,:) = global_roi.p(1,:) - ROI_up + ROI_side;
global_roi.p(3,:) = global_roi.p(1,:) - ROI_up;
global_roi.p(4,:) = global_roi.p(1,:) + ROI_side;

% 1 small roi ahead of crack

ROI_x = 10; % mm
ROI_y = 3.0;% mm
dist_up = ROI_y/2; % (mm) istance from crack

ROI_side = side*ROI_x;
ROI_up = up*ROI_y;
ROI_shift = side*2.2;
global_roi_ahead.p=zeros(4,2);
global_roi_ahead.p(1,:) = ROI_shift+notch_tip + ROI_up-dist_up*up;
global_roi_ahead.p(2,:) = global_roi_ahead.p(1,:) - ROI_up + ROI_side;
global_roi_ahead.p(3,:) = global_roi_ahead.p(1,:) - ROI_up;
global_roi_ahead.p(4,:) = global_roi_ahead.p(1,:) + ROI_side;



info.up=up;
info.side=side;
info.scale=scale;
info.pixelsize=1.0/scale/1000;
info.notch_tip=notch_tip;
info.notch_top=notch_top;
info.notch_bot=notch_bot;
end
