%%%%%%%%%%%%%%% Updated debug.m:
TOTAL_POWER_ITR= 50; % input('pow iter: ');

tol = 1.0e-3;

load trouble;
aUold = aU;
nsvd = size(aU,2);

% Get true singular values...
[Ut St Vt] = svd(full(sL)); St = diag(St);

figure(1); clf; 
ht = semilogy(St, 's-r');

St(nsvd)
St(nsvd)^50
St(nsvd)^25

%% Properties of initial aU
sum(aU.^2,1)

Snsvd = min(St(1:nsvd));  %% CHANGE THIS TO BE THE COARSE ESTIMATE


%Chakra, THE FOLLOWING IS NOT THE SPECTRUM.  aU is not orthogonal.
% guess for the spectrum before power-iteration
%aSold = diag(aU'*sL*aU);
%[m1,n1] = max(aSold);
%aSold(n1) = 1; % this is my fix. not necessary.

% power-iterations 
%aUold = aU;  % I did this above, instead

TPI = min(TOTAL_POWER_ITR, floor(log(nsvd * eps/tol)/log(Snsvd)))


aU = aUold;
for k = 1:TPI  %% CHANGE from TOTAL_POWER_ITR to TPI
  aU = sL*aU;
end
figure(101); 
%if (TOTAL_POWER_ITR == 25)
%  clf;
%end
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

%% NEW Compare singular values
figure(1); clf; hold on;
plot(St(1:nsvd), 'x-r');  % True singular values
plot(aS(1:nsvd), 'o-g');  % Approximate

figure(2); clf;
plot(aS(1:nsvd)-St(1:nsvd))

