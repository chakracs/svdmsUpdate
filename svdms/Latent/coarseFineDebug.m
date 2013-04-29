function [aU,aS,Ur,Sr,itr] = coarseFineDebug(sL,Ar,K,nSvd,Ur,Sr,st)
%function [aU,aS] = coarseFine(sL,Ar,K, nSVD)
%function [aU,aS] = coarseFine(sL,Ar,K, nSVD,Ur,Sr,st)
%
% using latent space affinities Ar, interpolation 
% kernels K and fine scale normalized affinities sL, 
% generate fine scale eigenvectors/eigenvalues.
% use power iteration on the fine scale eigenvectors.
% use only the top nSvd eigenvectors. the eigenvalues
% are computed by the rayleigh function.
%
%% calls: gramFixed.m
%
  USE_SPARSE_SVD = 0;  % 0 => use full svd on small problems, it is faster.
  POW_ITR_TR =  1e-4;%1e-3;
  TOTAL_POWER_ITR = 50;%25;%25 and 20 work.
  
  %======== STEP 1:  ===========
  if nargin<5 | prod(size(Ur) == 0)

    % scaling Ar. not necessary, done by setLatent.m
    %Ar = Ar/median((sum(Ar,1)));
    
    Dr = sum(Ar,2);
    
    sqrtD    = spdiags(Dr.^0.5, 0, length(Dr), length(Dr));
    sqrtDinv = spdiags(Dr.^ -0.5, 0, length(Dr), length(Dr));
    
    if issparse(Ar) & USE_SPARSE_SVD
      [Ur,Sr,Vr] = svds(sqrtDinv*Ar*sqrtDinv,...
                        min(nSvd,size(Ar,2)));
      Sr = diag(Sr);
    else
      [Ur,Sr,Vr] = svd(full(sqrtDinv * Ar * sqrtDinv));
      Sr = diag(Sr);
      % Crop to the first nSvd eigenvectors
      if nSvd < size(Ur,2)
        Ur = Ur(:,1:nSvd);
        Sr = Sr(1:nSvd);
      end
    end
  end

  %%% DEBUG MODE (May 04)
  %keyboard;
  
  %======== END STEP 1  ===========

  %======== STEP 2:  ===========  
  % power iteration to update eigenvectors.
  % orthogonalize with fixed gram procedure.
  %
  % eigen values are obtaineed by rayleigh coefficient
  % measure. can the rayleigh coefficients be ever negative?
  %
  % to do: stopping criterion for power iterations
  % calls: gramFixed (to orthogonalize power iterated basis)
  %
  
  % interpolate, that is extend Ur to fine scale
  aU = K*Ur;
  % normalize to be unit vectors
  for cp = 1:size(aU,2)
    aU(:,cp) = aU(:,cp)/norm(aU(:,cp));
  end
  aS = diag(aU'*sL*aU);
  figure(201); 
  plot(diag(aU'*aU),'-','linewidth',2); 
  set(gca,'fontsize',15); grid on;
  hold on; plot(aS,'r-','linewidth',2);  axis tight;
  title(sprintf('diag(U^T*U) iter 0'));
  pause;
  
  % orthogonalize right away (Apr'03, May 04)
  %[id,aU] = gramFixed(aU,0);    
  %aS = diag(aU'*sL*aU);
  %[aS,id] = sort(-aS);
  %aS = -aS;
  %aU = aU(:,id);

  % if one of the input arguments is the stationary 
  % distribution at the fine scale

  %if exist('st', 'var')
  %  aU(:,1) = sqrt(st(:));
  %end
  
  DEBUG = 0;
  itr = 1;
  if DEBUG
    fprintf('itr: ');  
  end

  keyboard;
  
  TOTAL_POWER_ITR = input('pow iter: ');
  converged = 0; 
  while (~converged)

    aUold = aU;

    for k = 1:TOTAL_POWER_ITR
      aU = sL*aU;

      % normalize to be unit vectors
      for cp = 1:size(aU,2)
	aU(:,cp) = aU(:,cp)/norm(aU(:,cp));
      end
      
      
      %[aU(:,1)'*aU(:,end) norm(aU(:,1)) norm(aU(:,end))]
      aUstats(k,1) = norm(aU(:,end));
      aUstats(k,2) = aU(:,1)'*aU(:,end) ;
      
      GS = 0;
      if GS
	% orthogonalize in space
	[id,aU] = gramFixed(aU,0);    
	%fprintf('norm aU: %d %f %f\n',itr,max(sum(aU.*aU,1)),min(sum(aU.*aU,1)));
	aS = aU'*sL*aU;
	% orthogonalize in space residuals 
	[us,ss,vs] = svd(aS);
	% update in space eigenvectors
	aU = aU*us;
	% grab in space eigenvalues
	aS = diag(ss);
	
	[aS,id] = sort(-aS);
	aS = -aS;
	aU = aU(:,id);
      end
      
    end
    %keyboard;
    figure(201); 
    plot(diag(aU'*aU),'-','linewidth',2); 
    set(gca,'fontsize',15); grid on;
    hold on; plot(aS,'r-','linewidth',2); axis tight;
    title(sprintf('diag(U^T*U) iter %d',itr));    
    %pause;
    keyboard;    
    if ~GS
      % orthogonalize in space
      %[id,aU] = gramFixed(aU,0);    
      aU = gramFixedModified(aU,0);          

      %fprintf('norm aU: %d %f %f\n',itr,max(sum(aU.*aU,1)),min(sum(aU.*aU,1)));
      aS = aU'*sL*aU;

      figure(202); 
      plot(diag(aU'*aU),'-','linewidth',2); 
      set(gca,'fontsize',15); grid on;  
      hold on; plot(-sort(-diag(aS)),'r-','linewidth',2); 
      title(sprintf('diag(U^T*U) iter %d after GS',itr));    
      %pause;
      
      %keyboard;
      % orthogonalize in space residuals 
      [us,ss,vs] = svd(aS);
      % update in space eigenvectors
      aU = aU*us;
      % grab in space eigenvalues
      aS = diag(ss);
      
      [aS,id] = sort(-aS);
      aS = -aS;
      aU = aU(:,id);
    end
    
    %figure(201); plot(diag(aU'*aU)); 
    figure(202); 
    plot(diag(aU'*aU),'-','linewidth',2); 
    set(gca,'fontsize',15); grid on;
    hold on; plot(aS,'r-','linewidth',2); 
    title(sprintf('diag(U^T*U) iter %d after GS + update',itr));    
    pause;

    ok = 0;
    if ok
      % reset values in aS that are > 1
      aS(aS > 1) = 1.0;
      % set lengths of eigenvectors to be 1
      lengthU = sum(aU.*aU,1);
      lidx    = lengthU > 1.0;
      if (lidx)
	aU(:,lidx) = aU(:,lidx) ./ (ones(size(aU,1),1)*sqrt(lengthU(lidx))) ;
      end
    end
    
    innerP = 1 - sum(aU .* aU,1)';
    %innerP = 1 - abs(sum(aU .* aUold,1))';    
    %fprintf('innerP aU: %f %f\n',max(sum(aU.*aUold,1)),min(sum(aU.*aUold,1)));
    converged = max(abs(innerP)) < POW_ITR_TR;
    
    if DEBUG 
      pause;
    end

    itr = itr + 1;
    
  end

  %======== END STEP 2  ===========  

  %fprintf('Coarse-Fine Iterations: %d\n',itr);
  
  
  return;
  

  %JUNK CODE

  %%ANGLE_TR = pi/36;
  %% the pow_itr_tr is roughly equal to (1 - cos(angle_tr))

  %GRAMS_ITR  = 50;
  %HALF_LIFE_TR = 50;

  
