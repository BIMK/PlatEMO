% Copyright (c) 2014, J.M. Hernandez-Lobato, M.W. Hoffman, Z. Ghahramani
% This function is adapted from the code for the paper
% Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
% Predictive Entropy Search for Efficient Global Optimization of Black-box
% Functions, In NIPS, 2014.
% https://bitbucket.org/jmh233/codepesnips2014
function [ l, sigma, sigma0 ] = sampleHypers(Xsamples, Ysamples, nSamples, fixhyp)
% Check if the hyper parameters are fixed
if isfield(fixhyp, 'l') && isfield(fixhyp, 'sigma') && isfield(fixhyp, 'sigma0')
    % Fix the hyper parameters
    l = repmat(fixhyp.l, [nSamples, 1]);
    sigma = ones(nSamples,1)*fixhyp.sigma;
    sigma0 = ones(nSamples,1)*fixhyp.sigma0;
else
    
    % We specify a gamma prior for the noise
    
    pNoise = prior_logunif();
    
    % We specify the likelihood
    
    lik = lik_gaussian('sigma2_prior', pNoise);
    
    % We specify a gamma prior for the lengthscales
    
    meanG = 0.5;
    varianceG = 0.1;
    alpha = meanG^2 / varianceG;
    beta = meanG / varianceG;
    pLengthscale = prior_gamma('sh', alpha, 'is', beta);
    
    % We specify a gamma prior for the amplitud
    
    pAmplitud = prior_logunif();
    
    % We specify the covariance function
    
    gpcf = gpcf_sexp('lengthScale', 0.3 * ones(1, size(Xsamples, 2)), ...
        'lengthScale_prior', pLengthscale, 'magnSigma2_prior', pAmplitud);
    
    % We infer the hyper-parameters
    
    gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-10);
    
    % We sample from the hyper-parameters
    
%     multiplier = 3;
%     burnin = 50;
%     [gp_rec, g, opt] = gp_mc(gp, Xsamples, Ysamples, 'nsamples', (nSamples + burnin) * multiplier, 'display', 0);
%     gp_rec = thin(gp_rec, multiplier * burnin + 1, multiplier);
    
    opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
    gp_rec=gp_optim(gp,Xsamples,Ysamples,'opt',opt);
    
    % We get the samples
    
    sigma0 = gp_rec.lik.sigma2;
    l = 1 ./ gp_rec.cf{ 1 }.lengthScale.^2;
    sigma = gp_rec.cf{ 1 }.magnSigma2;
end