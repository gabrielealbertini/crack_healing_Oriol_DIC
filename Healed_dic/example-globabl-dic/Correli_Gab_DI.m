clear all; close all; clc
% setting some default plot options
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',18);
set(0,'defaultPatchEdgeColor',[0.5,0.5,0.5]);
set(0,'defaultPatchFaceColor','none');
set(0,'defaultFigureColormap',hot);

% specify the path to the Correli folder
Correlipath = '/home/ga288/sources/correli/'; %fullfile(pwd,'..','..','..');

% this adds the correct folder to the matlab path
addpath(fullfile(Correlipath,'ml_lib_core'));

plot_movie = 1;

% ==========================
% DIC Options
% ==========================
% global options (used for each coarse grain step)
options.verbose = 0;      % how much info to write during correlation
options.plotres = 1;      % show intermediate residual plots

options.mechreg     = 10*[1 1 1];
options.initreg     = 0;
options.trustregion = 0;
options.bestit      = false;
options.divmax      = 0;

% ====================================================
% Experimental Data
% ====================================================

% images
imagepath = '../../Healed';
imagefiles = dir(fullfile(imagepath,'*.jpg'));
imagefiles = {imagefiles.name}';
images = strcat(imagepath,filesep,imagefiles);
imini = 1;
imdec = 1;
imfin = numel(images)-imdec;

% create node coordinates from a to b, with the edge elements a factor c
% bigger than the internal elements

load('ROI_selection.mat');


roi = global_roi;
%  p1 --------- p3
%   |            |
%   |            |
%  p4 --------- p2
x1 = min(roi.p(2,2),roi.p(1,2));
y1 = min(roi.p(2,1),roi.p(1,1));

x2 = max(roi.p(2,2),roi.p(1,2));
y2 = max(roi.p(2,1),roi.p(1,1));

x3 = roi.p(3,2);
y3 = roi.p(3,1);
x4 = roi.p(4,2);
y4 = roi.p(4,1);

lx = norm(roi.p(3,:)-roi.p(1,:));
ly = norm(roi.p(4,:)-roi.p(1,:));
dir1 = roi.p(3,:)-roi.p(1,:);
ang = -atan(dir1(1)/dir1(2));


% ang = atan2(21,1727);
%
% Initialization option
flag_ini     = 1; %0 = no existing result,
flag_mes_ini = 3; % mesh.tag of initial mesh
flag_mes_cur = 4; % mesh.tag of current mesh

file_init = ['DI_results_mesh',num2str(flag_mes_ini),'_D.mat'];
if flag_ini ~= 0
    eval(['load(''',file_init,''')'])
    Mesh_old = Mesh;
    Udic_old = Mesh.Udic;
end

file_save = ['DI_results_mesh',num2str(flag_mes_cur),'_D.mat'];
if ~exist(file_save,'file')
    make_mesh = 1;
else
    make_mesh = 0;
    eval(['load(''',file_save,''')'])
end

% Material parameters for titanium
nu=0.3;              % Poisson's ratio
Young=1;              % Young's modulus (MPa)
hyp = 'Pstress-iso';
C = mechreg_stiff_tensor(hyp,nu);
% hypL = 'Laplacian-2D';
% C = mechreg_stiff_tensor(hypL);

% load the reference image
f = rgb2gray(imread(images{1}));
[n,m,~] = size(f);
frgb(:,:,1) = f;
frgb(:,:,2) = f;
frgb(:,:,3) = f;
f = single(f);
dynamicrange = max(f(:)) - min(f(:));
fprintf('Dynamic range of f: %g [GV] \n',dynamicrange);
imagesc(frgb)
axis image
%
% ==========================
% Mesh
% ==========================

if (flag_ini == 0) || (make_mesh == 1)
    
    % number of nodes (x,y)
    if flag_mes_cur == 5
        Nnod = [160,160];
        Mesh.tag  = 5;
        c = 1;%sqrt(2);
        % number of nodes (x,y)
    elseif flag_mes_cur == 4
        Nnod = [40,80];
        Mesh.tag  = 4;
        c = 1;%sqrt(2);
        % number of nodes (x,y)
    elseif flag_mes_cur == 3
        Nnod = [20,40];
        Mesh.tag  = 3;
        c = 1;%sqrt(2);
    elseif flag_mes_cur == 2
        Nnod = [10,20];
        Mesh.tag  = 2;
        c = 1;
    elseif flag_mes_cur == 1
        Nnod = [5,10];
        Mesh.tag  = 1;
        c = 1;
    end
    %
    if flag_mes_cur ~= 10
        Mesh.Nnod= Nnod;
        L = c./(2+Nnod-3);
        nodex = x1 + (x2-x1)*[0, linspace(L(1),1-L(1),Nnod(1)-2), 1];
        nodey = y1 + (y2-y1)*[0, linspace(L(2),1-L(2),Nnod(2)-2), 1];
        
        % Let's create a structured mesh for regular DIC
        Mesh = mesh_gen_structured(nodey,nodex,'T3');
%         Mesh.coordinates(:,2) = Mesh.coordinates(:,2)+Mesh.coordinates(:,1)*tan(ang);
    else
        load('BM_E_mesh_1.mat')
        x1 = min(Mesh.coordinates(:,2));
        x2 = max(Mesh.coordinates(:,2));
        y1 = min(Mesh.coordinates(:,1));
        y2 = max(Mesh.coordinates(:,1));
    end
    %
    % plot the mesh on the image
    figure('Position',[50 50 600 400])
    subplot(121)
    imagesc(frgb)
    h = mesh_plot(Mesh,'Show','NodeArea','Wireframe',true);
    daspect([1 1 1]);
    title('Element Area per node')
    colorbar
    subplot(122)
    h = mesh_plot(Mesh,'Indices',true);
    daspect([1 1 1]);
    title('Numbering')
    colorbar
    disp('NEW MESH CREATED')
end

xcen = (x1+x2)/2;
ycen = (y1+y2)/2;
xsiz = 2^(floor(log2(x2-x1)));
ysiz = 2^(floor(log2(y2-y1)));
opt.subsetsize = [ysiz,xsiz];
opt.subpixel   = 0;

% options.Kelast      = mechreg_stiff_bulk(Mesh,C);
% 
EdgeL = find(abs(Mesh.coordinates(:,1)-y1)<2);
EdgeR = find(abs(Mesh.coordinates(:,1)-y2)<2);
EdgeT = find(abs(Mesh.coordinates(:,2)-x1)<2);
EdgeB = find(abs(Mesh.coordinates(:,2)-x2)<2);
Mesh.groupsnodes.mechreg_corner = [EdgeL]';
Mesh.groupsnodes.mechreg_edge   = [EdgeR]';
Mesh.groupsnodes.mechreg_bulk   = setdiff(Mesh.nodes,[EdgeL;EdgeR]');
% Mesh.groupsnodes.mechreg_bulk   = Mesh.nodes;

% Initialize
% =====================
if flag_mes_ini ~= flag_mes_cur && (flag_ini ~= 0)
    Udic = mesh_interp(Mesh_old,Udic_old,Mesh,'linear');
elseif flag_ini ~= 0
    Udic = Udic_old;
end

Mesh_ini = Mesh;
fz = f;


%%
% Direct calculation

% options.mechreg = [100 100 100]; % 0 but 10 for last TWO meshes

% file_save = ['zoom_results_mesh',num2str(flag_mes_cur),'_D.mat'];
% file_rest = ['knitting_results_mesh',num2str(flag_mes_ini),'_D.mat'];
% if exist(file_rest,'file')
%     eval(['load(''',file_rest,''')'])
%     Udic = Mesh.Udic;
% end

co = 0;
for im = imini:imdec:imfin %numel(images)-1
    co = co+1;
    % load the images
    go = single(rgb2gray(imread(images{im})));
    go = go(:,:,1);
    g = single(rgb2gray(imread(images{im+imdec})));
    g = g(:,:,1);
    if flag_ini == 0
        if co==1
            Uinit = zeros(length(Mesh.nodes),3);
%             Uinit = zeros(prod(Nnod),3);
        else
            Uinit = Udic_d{co-1};
        end
    else
        Uinit = Udic{co};
    end
    % ==========================
    % Perform DIC
    % ==========================
    
    options.itermax = 1000;
    options.conv_limit = 1.e-2;
    options.outputs = 'FRN';    % store the residual image R
    cor = correlate_2D(fz,g,Mesh_ini,Uinit,options);
    
    % readable convergence state
    cstate{1,1} = '';
    cstate{2,1} = ', max. iter.';
    cstate{3,1} = ', div.';
    
    % print a status line (red if not converged properly)
    fid = (cor.s(1) ~= 0) + 1;
    fprintf(fid,'%3d   it%3d/%3d, r%9.2e, dU%9.2e%s\n',im,cor.s(3),cor.s(2),cor.s(4),cor.s(6),cstate{cor.s(1)+1});
    
    Udic_d{co} = cor.U;
    Rdic_d{co} = cor.N;
    RMSres_d(co) = cor.r;
    resi_d{co}   = cor.R;
    
end
%
%%

figure
plot(RMSres_d,'r+-')
xlabel('Picture #')
ylabel('RMS residual, gray levels')
file_print = strrep(file_save,'.mat','_RMSres_D.png');
eval(['print -dpng ''',file_print,''''])

Mesh        = Mesh_ini;
Mesh.Udic   = Udic_d;
Mesh.R      = Rdic_d;
Mesh.RMSres = RMSres_d;
% Mesh.resi   = resi_d;

save(file_save,'Mesh')
%
%%
savename = strrep(file_save,'.mat','_res.mat');
save(savename,'-v7.3','resi_d');