addpath('Latent');
addpath('data');

clear all
close all

FALSE = (0 == 1);
TRUE = ~FALSE;

MAKE_AFF = FALSE;

Nid = [
    16 %1
    24 %2 
    32 %3
    48 %4
    63 %5
    64 %6
    65 %7
    80 %8
    100 %9
    120 %10
    124 %11
    127 %12
    128 %13
    129 %14
    130 %15
    160 %16
    255 %17
    256 %18
    257 %19
    300 %20
    362 %21
    511 %22
    512 %23
    513 %24
];


if TRUE %FALSE
  Tres = zeros(length(Nid),3);
  Tres(:,1) = Nid(:);
else
  load 'Tres'
  figure(1); clf; 
  % figure(1); hold on;
  plot(Tres(1:14,1).^2, Tres(1:14,2), '-*b');
end
T = zeros(length(Nid),2);
T(:,1) = Nid(:);


% for cid=1:15
% cid=15;
cid=4;

N = Nid(cid)
Tres(cid,1) = N;

dLogHalfAll = [];

%=== affty matrix parameters
afftyPar.sizeIm  = [N N];
afftyPar.dsThres = 1.1;
afftyPar.dsSupp  = 3.1; 
afftyPar.rho     = 1.5; 
beta0 = 40;
half0 = beta0/4;  
id0 = 2;

%=== create membrane. load from cache
if FALSE
  load(sprintf('im_%d',N));
else
  % or create a new one
  im = mkMembrane(N,0); % no foreground % save(sprintf('im_%d',N),'im');
  MAKE_AFF = TRUE;
  %im = mkMembrane(N,1); % with foreground

  % show the membrance
  %figure(2); clf; showIm(im); pause(.1);
end

sizeIm = [N N];%size(im);

%Pts = ones(prod(sizeIm),2);
%Pts(:,1) = im(:);
%figure(1);clf;showIm(im);

% affinity matrix. 
% load from the cache
if ~MAKE_AFF
  load(sprintf('affty_%d',N));
else
  % or create a new affinity matrix
  %profile on -detail 'builtin'
  %tic
  A = shiftAffty(im,afftyPar.rho); 
  %toc
  save(sprintf('../../svdms/data/affty_%d',N),'A');
end

nSVD = 51;
tic; [U,S] = svdms(A,nSVD,sizeIm); toc
