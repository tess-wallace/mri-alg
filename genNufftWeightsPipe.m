function [wi, actmaxerr, psf_final] = genNufftWeightsPipe(ktr, N, maxerr)

% genNufftWeightsPipe generates a density compensation function for
% non-Cartesian reconstruction using the iterative algorithm proposed by
% Pipe and Menon
%
% Inputs
% ktr       radial k-space locations [nCol*nSpokes, nDim]
% N         reconstructed image size
% maxerr    maximum errror
% maxiter   maximum number of iterations
%
% Outputs
% wi        weights [nCol*nSpokes]
% actmaxerr computed error
%
% Onur Afacan

Nd = size(ktr,2);

kbw = repelem(4,Nd);
% kbw = repelem(3,Nd); % less memory...

nufftst = Gnufft({ktr, N, kbw, 2*N, N/2}); % generate st object

wi = ones(length(ktr), 1);
P = nufftst.arg.st.p; % scale factor?

goal = inf;
iter = 0;
saver = zeros(200,1);
maxav10e = 1;

while (max(abs(goal-1)) > maxerr && abs(maxav10e) > 1e-3 && iter<100)
  
    iter = iter + 1;
   
    goal = P * (P' * wi);
    
    errx(iter)=max(abs(goal-1));
    
    fprintf( 'Iteration %d: error = %.5f\n', iter, errx(iter));
        
    % psf1=nufftst'*wi;
    
    % psf2(:,:,:,iter)=reshape(psf1./max(abs(psf1(:))),N);

    wi = wi ./ real(goal);

    if iter > 1000
        warning('Iteration stuck?');
    end
    
    saver(iter) = max(abs(goal-1));
  
    % average of last 10 errs:
    if ( iter > 100 )
        av10e = mean(saver((iter-10):iter));
        maxav10e = max(abs(saver((iter-10):iter)-av10e));
    else
        maxav10e = 1;
    end
  
end

% Calculate final point spread function

psf_final = nufftst'*wi;
    
psf_final = reshape(psf_final./max(abs(psf_final(:))),N);


if ( max(abs(goal-1)) > maxerr && 0 )
    
    fprintf('NUFFT Pipe converged to %.5f error in %d iterations, instead of requested %.5f (maxav10e=%.6f)\n', max(abs(goal-1)), iter, maxerr, abs(maxav10e));
end

%scale = nufftst.arg.st.sn(end/2,end/2)^(-2) / fov^2 / prod(nufftst.arg.st.Kd) * prod(nufftst.arg.st.Nd);
%wi = wi * scale;

%wi = reshape(wi, size(ktr));

actmaxerr = max(abs(goal-1));


