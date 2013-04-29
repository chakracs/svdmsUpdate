
TOTAL_POWER_ITR = input('pow iter: ');

load trouble;

% guess for the spectrum before power-iteration
aSold = diag(aU'*sL*aU);
[m1,n1] = max(aSold);
aSold(n1) = 1; % this is my fix. not necessary.

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
