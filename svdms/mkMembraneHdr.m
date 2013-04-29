N = [];
while isempty(N)
  N = input('N: ');
end


%=== affty matrix parameters
afftyPar.sizeIm  = [N N];
afftyPar.dsThres = 1.1;
afftyPar.dsSupp  = 3.1; 
afftyPar.rho     = 1.5; 
beta0 = 40;
half0 = beta0/4;  
id0 = 2;

for kk = 1:10
  % or create a new one. no foreground 
  im = mkMembrane(N,0); 
  % show the membrance
  %figure(2); clf; showIm(im); pause(.1);
  
  fprintf('saving newdata/im_%d_%d.mat ...',N,kk);
  save(sprintf('newdata/im_%d_%d',N,kk),'im');
  fprintf('done\n');
  
  A = shiftAffty(im,afftyPar.rho); 
  
  fprintf('saving newdata/affty_%d_%d.mat ...',N,kk);
  save(sprintf('newdata/affty_%d_%d',N,kk),'A');
  fprintf('done\n');
end
