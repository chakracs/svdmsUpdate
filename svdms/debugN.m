%addpath('../speig');
addpath('Latent');
addpath('newdata');

clear all
close all

N = 16;
id = 4;

%=== affty matrix parameters
afftyPar.sizeIm  = [N N];
afftyPar.dsThres = 1.1;
afftyPar.dsSupp  = 3.1; 
afftyPar.rho     = 1.5; 
beta0 = 40;
half0 = beta0/4;  
id0 = 2;


load(sprintf('im_%d_%d',N,id));
load(sprintf('affty_%d_%d',N,id));

sizeIm = size(im);
Pts = ones(prod(sizeIm),2);
Pts(:,1) = im(:);
%figure(1);clf;showIm(im);

% scaling the affinity matrix
Dorig = full(sum(A,1))';
scl = median(full(Dorig));
A = A/scl;
Dorig = full(sum(A,1))';

clear hier
hier{1}.A = A;
hier{1}.sizeIm = sizeIm;

%tic
%profile on -detail 'builtin'
lev = 1;
fprintf(2, ' Latent:  lev %d, size %d\n', lev, prod(hier{lev}.sizeIm));
while prod(hier{lev}.sizeIm) > 180
  if lev ==1
    logpow = 1;%(original=1)
  else
    logpow = 2;
  end
  lev = lev+1;
  [hier{lev}.L hier{lev}.A hier{lev}.K hier{lev}.R hier{lev}.st ...
   hier{lev}.W hier{lev}.rbinNhbr sId1 sMp1] = ...
      buildLatent(hier{lev-1}.A, logpow, hier{lev-1}.sizeIm);
  hier{lev}.sizeIm = [size(hier{lev}.A, 1) 1];
  fprintf(2, ' Latent:  lev %d, size %d\n', lev, prod(hier{lev}.sizeIm));
  if hier{lev}.sizeIm(1) < 51
    clear hier{lev};
    lev = lev-1;
    break;
  end
end 
nLev = lev;

if nLev > 1
  [hier{nLev-1}.U, hier{nLev-1}.S,hier{nLev}.U,hier{nLev}.S,...
   hier{nLev-1}.itr] = coarseFineDebug(hier{nLev}.L,hier{nLev}.A,hier{nLev}.K,51);
  fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', [nLev, size(hier{nLev}.U)]);
  fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', ...
          [nLev-1, size(hier{nLev-1}.U)]);
  for lev = nLev-1:-1:2  
    [hier{lev-1}.U, hier{lev-1}.S,Pp,Qq,hier{lev-1}.itr] ...
        = coarseFineDebug(hier{lev}.L,hier{lev}.A,hier{lev}.K, ...
                     51,hier{lev}.U,hier{lev}.S);
    fprintf(2, ' Prolong: lev %d, size %d, nU %d\n', ...
            [lev-1, size(hier{lev-1}.U)]);
  end
end
%profile report /var/tmp/summary2
%T(cid,2) = toc;
%T(cid,:)


%end

ok = 0;
if ok
%% expts with aUold, aU
addpath('Latent');
load aUold25;
load aUold50;

aU25 = gramFixedModified(aUold25,0);
aU50 = gramFixedModified(aUold50,0);

for i = 1:size(aU25,2)
  nor(i,1) = norm(aUold25(:,i));
  nor(i,2) = norm(aUold50(:,i));  
end

figure; 
plot(nor(:,1),'linewidth',2); 
hold on;
plot(nor(:,2),'r','linewidth',2);
set(gca,'fontsize',15); grid on;

aUold = aUold25;
r = abs(aUold'*aUold);
r = r - diag(diag(r));
figure; showIm(r);
aU = aU25;
r = abs(aU'*aU);
r = r - diag(diag(r));
figure; showIm(r);



for k = 1:size(r,1)
end

end

ok = 0;
if ok

  TOTAL_POWER_ITR = input('pow iter: ');

  load trouble;

  % guess for the spectrum before power-iteration
  aSold = diag(aU'*sL*aU);
  [m1,n1] = max(aSold);
  aSold(n1) = 1;
  
  % power-iterations 
  aUold = aU;
  for k = 1:TOTAL_POWER_ITR
    aU = sL*aU;
  end
  figure(101); 
  if (TOTAL_POWER_ITR == 25)
    clf;
  end
  %plot(sum(aU.^2,1),'x-','linewidth',2); hold on;
  if (TOTAL_POWER_ITR == 25)
    semilogy(sum(aU.^2,1),'x-','linewidth',2); hold on;     
    %plot(aSold.^(TOTAL_POWER_ITR),'or-','linewidth',2);
    semilogy(aSold.^(TOTAL_POWER_ITR),'or-','linewidth',2);  
  else
    semilogy(sum(aU.^2,1),'x-g','linewidth',2); hold on;      
    %plot(aSold.^(TOTAL_POWER_ITR),'or-','linewidth',2);
    semilogy(aSold.^(TOTAL_POWER_ITR),'om-','linewidth',2);  
  end
  set(gca,'fontsize',15); grid on;  axis tight;
  
  figure; showIm(aU'*aU); 
  % orthgonalize
  aU = gramFixedModified(aU,0);          
  figure; showIm(aU'*aU);
  %fprintf('norm aU: %d %f %f\n',itr,max(sum(aU.*aU,1)),min(sum(aU.*aU,1)));
  aS = aU'*sL*aU;
  figure; showIm(aS);
  
  % rest of the code is not really necessary. just 
  % reproducing the pseudo-code in the paper.
  
  [us,ss,vs] = svd(aS);
  % update in space eigenvectors
  aU = aU*us;
  % grab in space eigenvalues
  aS = diag(ss);
  
  [aS,id] = sort(-aS);
  aS = -aS;
  aU = aU(:,id);

  
end
