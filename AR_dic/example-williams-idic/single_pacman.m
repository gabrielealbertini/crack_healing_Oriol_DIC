clear ; close all ; clc

% setting some default plot options
dpi = 200;
fontsize = 18;
markersize = 12;
Ncolor = 21;
set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',fontsize);
set(0,'DefaultLineMarkerSize',markersize);
set(0,'DefaultPatchEdgeColor',[0.5,0.5,0.5]);
set(0,'DefaultPatchFaceColor','none');
set(0,'DefaultPatchMarkerSize',markersize);
set(0,'DefaultFigureColormap',parula(64));
% =======================================================

if isunix
    CorreliPath = 'put_your_path_here';
else
    CorreliPath = 'D:\Imagerie\Correli\';
end

% this adds the correct folders to the matlab path
addpath(fullfile(CorreliPath,'ml_lib_core'));

% Read the image series
% ======================================================
imagepath = 'D:\Imagerie\Gabriele\';
images = dir([imagepath, filesep, 'image_*.bmp']);
images = {images.name}.';
N = numel(images);

% number of steps (inc 267 doesn't have c and d)
Ninc = 266;

% Pixel size (meter)
pix2m = 50e-3 / 1168;

% reference image
f = imread(fullfile(imagepath,images{1}));
[n, m] = size(f);

% Crack definition
% ======================================================

if false
    % draw a crack path by hand, note the orientation matters, the positive
    % crack direction is from the first point to the second.
    
    % displaying some results
    figure('Position',[150 150 1200 800])
    colormap(gray)
    imagesc(f)
    daspect([1 1 1]);
    colorbar
    
    A = [0.5*m, 0.5*n; 0.95*m, 0.5*n];
    title('Move the crack path, confirm with double click')
    [A, P] = selectarea(A,'Polyline');
    
    % copy past this below to save the crack path
    crackpath = crack_makepath(A,1,2);
    save('crackpath.mat','crackpath')
else
    load('crackpath.mat','crackpath');
end

% initial crack length
a0 = 5.0e-3 /pix2m;

% initialize the crack
crackpos = crack_position(crackpath,a0);
crackang = crack_angle(crackpath,a0);


% Material properties
% ======================================================

% Material parameters for titanium (change for 2024)
nu=0.33;              % Poisson's ratio
Young=70e3;          % Young's modulus (MPa)
G=Young/(1+nu)/2;     % Shear modulus (MPa)
kappa=(3-nu)/(1+nu);  % Bulk modulus for the Plane stress state


% Cracktip Mesh
% ======================================================

% create a cracktip mesh
Radii  = linspace(20,200,4); % 4 concentric radii 
ElSize = [5, 10, 25, 50];    % the element size at each radius
Mesh = crack_makemesh_gmsh(crackpos,crackang,Radii,ElSize,'CrackOpen',20);

% create a contour mesh for pretty plots
Contour = mesh_contour(Mesh);

figure('Position',[150 150 1200 800])
colormap(gray)
imagesc(f)
daspect([1 1 1]);
colorbar
hold on
mesh_plot(Mesh,'Wireframe',true,'EdgeColor','y');
%mesh_plot(Contour,'Wireframe',true,'EdgeColor','b','Linewidth',2);
plot(crackpos(:,1),crackpos(:,2),'.r');

caxis([0 255])
drawnow



% Williams parameters
% ======================================================

% Store some data in the Mesh, to be used by Williams_S
Mesh.kappa    = kappa;
Mesh.order    = [ -1, 3];
Mesh.crackpos = crackpos;
Mesh.crackang = crackang;
Mesh.Rnorm    = Radii(4); % normalize internal radial coordinate with Rnorm

% number of modes and degrees of freedom
Nmod = Mesh.order(2) - Mesh.order(1) + 1;
Ndof = 2*Nmod;

% mode list
modelist = repmat(Mesh.order(1):Mesh.order(2),2,1);

% some scaling for later use (to scale the williams amplitudes to physical
% dimensions)
scal  = 2*G*sqrt(2*pi)*(pix2m.^(1-modelist*.5)); % scale factor

% scaling correction when using a normilized radial coordinates
scal_R = Mesh.Rnorm.^(modelist./2);
scal = scal./scal_R;

% pointers to certain important parameters
Im1 = find(modelist == -1,2);  % amplitudes related to the crack tip shift
Ip1 = find(modelist ==  1,2);  % amplitudes related  SIF
Ip2 = find(modelist ==  2,2);  % amplitudes the T stress (and rotation)

% string for plotting
modstr{1,1} = 'Mode I';
modstr{2,1} = 'Mode II';

% Correlate_2D.m options
options.verbose = 0;
options.plotres = 0;
options.outputs = 'RN';
options.trustregion = 0;
options.initreg = 0;
options.threads = 1;
   
% readable convergence state
cstate{1,1} = '';
cstate{2,1} = ', max. iter.';
cstate{3,1} = ', div.';

% Put everything in a big structre
% ======================================================

% Store some data in D structure
D(1).Mesh = Mesh;
D(1).cor.U  = zeros(Ndof,1);
D(1).options = options;
D(1).Mesh = Mesh;
D(1).Ip1 = Ip1;
D(1).Ip2 = Ip2;
D(1).Im1 = Im1;
D(1).scal = scal;
D(1).cstate = cstate;
D(1).crackpath = crackpath;
D(1).a = a0;
D(1).init = zeros(Ndof,1);

% options for the crack-tip loop
D(1).convcrit = 1;
D(1).maxit = 20;
D(1).trustregion = 10;
D(1).dampen = 0.6;
D(1).alpha = 0;

% initialize some fields
D.cor = [];
D.da = [];
D.it = [];
D.K1 = [];
D.K2 = [];
D.TS = [];


% Loop over the increments
% ======================================================

figure('Position',[150 150 1200 800]);
colormap(gray)
hdl_img = imagesc(f);
daspect([1 1 1]);
hold on
hdl_cnt = mesh_plot(Contour,'Wireframe',true,'EdgeColor','b','Linewidth',2);
hdl_tip = plot(crackpos(:,1),crackpos(:,2),'.r');
caxis([0 255])
drawnow

hdl_fig = gcf;
hdl_ax = gca;

results = zeros(Ninc,5);
for inc = 1:Ninc
    
    fprintf('\nIncrement %d --------------------\n',inc);
    
    refim_name = sprintf('.*%03d%s.tif',inc,'c');
    defim_name = sprintf('.*%03d%s.tif',inc,'a');
    
    refim = find(not(cellfun(@isempty,regexp(images,refim_name))),1);
    defim = find(not(cellfun(@isempty,regexp(images,defim_name))),1);
    
    % load the reference image
    f = imread(fullfile(imagepath,images{refim}));
    set(hdl_img,'CData',f)
    [n,m] = size(f);
    
    % load the deformed image
    g = imread(fullfile(imagepath,images{defim}));
    
    % convert to single precision
    f = single(f);
    g = single(g);

    % load the images in the structure
    D(1).f = f;
    D(1).g = g;
    
    % correlate for one pacman
    D(1) = correlate_williams(D(1));
    
    % store some results
    results(inc,1) = D.a;
    results(inc,2) = D.da;
    results(inc,3) = D.K1;
    results(inc,4) = D.K2;
    results(inc,5) = D.TS;
    
    Contour = mesh_contour(D.Mesh);
    crackpos = crack_position(crackpath,D.a);
    
    % update the figure
    set(0,'CurrentFigure',hdl_fig);
    set(hdl_fig,'CurrentAxes',hdl_ax);
    delete(hdl_cnt);
    delete(hdl_tip);
    hdl_cnt = mesh_plot(Contour,'Wireframe',true,'EdgeColor','b','Linewidth',2);
    hdl_tip = plot(crackpos(:,1),crackpos(:,2),'.r');
    caxis([0 255])
    
    str{1} = sprintf('inc %d',inc);
    str{2} = sprintf('a = %9.2e [mm]',D(1).a*pix2m*1e3);
    str{3} = sprintf('K_{I} = %9.2e, K_{II} = %9.2e',D(1).K1,D(1).K2);
    title(str);
    drawnow
    
end





