clear all
close all
clc
%
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',18);
set(0,'defaultPatchEdgeColor',[0.5,0.5,0.5]);
set(0,'defaultPatchFaceColor','none');
set(0,'defaultFigureColormap',parula);
%
flag_mes_cur = 1;
file_save = ['DI_results_mesh',num2str(flag_mes_cur),'_D.mat'];
eval(['load(''',file_save,''')'])
radix = strrep(file_save,'_D.mat','');
file_strain = [radix,'_strains.mat'];
%
Udic = Mesh.Udic;
Rdic = Mesh.R;
Xdic = Mesh.coordinates(:,1);
Ydic = Mesh.coordinates(:,2);
Odic = ones(size(Ydic));
%

imagepath = '../../AR';
imagefiles = dir(fullfile(imagepath,'*.jpg'));% .bmp or .tif //uncompressed
imagefiles = {imagefiles.name}';
images = strcat(imagepath,filesep,imagefiles);

frgb = imread(images{1});
[n,m,~] = size(frgb);
% frgb(:,:,1) = fz;
% frgb(:,:,2) = fz;
% frgb(:,:,3) = fz;
% fz = single(fz);
% dynamicrange = max(fz(:)) - min(fz(:));
% fprintf('Dynamic range of f: %g [GV] \n',dynamicrange);
imagesc(frgb)
hold on 
axis image

dec  = 1;
imax = size(Rdic,2);
% m    = 1024;
% n    = 1024;
%
co = 0;
for idx = 1:dec:imax
    co = co+1;
    h = figure;
    box on
    set(h,'Position',[350 150 700 500])
    correli_fig(Mesh,Udic{idx},Udic{idx}(:,2),m,n)
    title(sprintf('Picture # %03d -- u_x, pixels',idx))
    %pause(1)
    frame = getframe(h);
    img{co} = frame2im(frame);
    close(h)
end
%
file_print = strrep(file_save,'.mat','_ux.gif');
for idx = 1:co
    [A,map] = rgb2ind(img{idx},256);
    if idx == 1
        imwrite(A,map,file_print,'gif','LoopCount',Inf,'DelayTime',.25);
    else
        imwrite(A,map,file_print,'gif','WriteMode','append','DelayTime',.25);
    end
end
%
co = 0;
for idx = 1:dec:imax
    co = co+1;
    h = figure;
    box on
    set(h,'Position',[350 150 700 500])
    correli_fig(Mesh,Udic{idx},Udic{idx}(:,1),m,n)
    title(sprintf('Picture # %03d -- u_y, pixels',idx))
    %pause(1)
    frame = getframe(h);
    img{co} = frame2im(frame);
    close(h)
end
%
file_print = strrep(file_save,'.mat','_uy.gif');
for idx = 1:co
    [A,map] = rgb2ind(img{idx},256);
    if idx == 1
        imwrite(A,map,file_print,'gif','LoopCount',Inf,'DelayTime',.25);
    else
        imwrite(A,map,file_print,'gif','WriteMode','append','DelayTime',.25);
    end
end
%
% co = 0;
% for idx = 1:dec:imax
%     co = co+1;
%     h = figure;
%     coef = [Odic(:) Xdic(:) Ydic(:)]\Udic{idx}(:,2);
%     Unaff = Udic{idx}(:,1) - [Odic(:) Xdic(:) Ydic(:)]*coef;
%     set(h,'Position',[350 150 700 500])
%     correli_fig(Mesh,Udic{idx},Unaff,m,n)
%     title(sprintf('Picture # %03d -- u_x (non affine), pixels',idx))
%     pause(1)
%     frame = getframe(h);
%     img{co} = frame2im(frame);
%     close(h)
% end
% %
% file_print = strrep(file_save,'.mat','_ux_naff.gif');
% for idx = 1:co
%     [A,map] = rgb2ind(img{idx},256);
%     if idx == 1
%         imwrite(A,map,file_print,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,file_print,'gif','WriteMode','append','DelayTime',1);
%     end
% end
% %
co = 0;
for idx = 1:dec:imax
    co = co+1;
    h = figure;
    set(h,'Position',[350 150 700 500])
    correli_fig(Mesh,Udic{idx},Rdic{idx},m,n)
    title(sprintf('Picture # %03d -- residual, gray levels',idx))
    %pause(1)
    frame = getframe(h);
    img{idx} = frame2im(frame);
    close(h)
end
%
file_print = strrep(file_save,'.mat','_res.gif');
for idx = 1:co
    [A,map] = rgb2ind(img{idx},256);
    if idx == 1
        imwrite(A,map,file_print,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,file_print,'gif','WriteMode','append','DelayTime',1);
    end
end
%
% %%
% if ~exist(file_strain,'file')
%     Edic = mesh_strain(Mesh,Mesh.Udic,'StrainDefinition','GreenLagrange');
%     eval(['save(''',radix,'_strains.mat'',''Edic'')'])
% else
%     eval(['load(''',file_strain,''')'])
% end
% %
% co = 0;
% for idx = 1:dec:imax
%     co = co+1;
%     h = figure;
%     box on
%     set(h,'Position',[350 150 1000 600])
% %     maxe = 5*std(Edic{idx}(:,2));
% %     ind = find(abs(Edic{idx}(:,2))>maxe);
% %     Edic{idx}(ind,2) = NaN;
%     correli_fig(Mesh,Udic{idx},Edic{idx}(:,2),m,n)
%     title(sprintf('Picture # %03d -- E_{xx}',idx))
%     pause(1)
%     frame = getframe(h);
%     img{co} = frame2im(frame);
%     close(h)
% end
% %
% file_print = strrep(file_save,'.mat','_exx.gif');
% for idx = 1:co
%     [A,map] = rgb2ind(img{idx},256);
%     if idx == 1
%         imwrite(A,map,file_print,'gif','LoopCount',Inf,'DelayTime',.25);
%     else
%         imwrite(A,map,file_print,'gif','WriteMode','append','DelayTime',.25);
%     end
% end