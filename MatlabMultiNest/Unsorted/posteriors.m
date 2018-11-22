function posteriors(post_samples, wp, parnames)

% function posteriors(post_samples, wp)
%
% This function displays joint and/or marginalized posterior
% distributions, calculates summary statistics, etc.
%
% post_samples is an Nx(npars+2) array where N is the number of
% posterior samples, npars is the number of parameters, and the last two
% columns in this array contain the values of logL and logPosterior
% (=prior weighted likelihood/evidence).
%
% wp is a vector containing the parameter(s) for which you want to create
% the posterior e.g., wp=1 will provide a posterior only on the 1st
% parameter; wp=[1 3] will provide a 2D posterior on the 1st and 3rd
% parameters. For the moment wp can only contain a maximum of two params.
%
% parnames is a cell array of string names corresponding to the
% parameters specified by wp.  Eg., parnames = {'b', 'slope'}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npars = size(post_samples,2)-2;
lwp = length(wp);
lparnames = length(parnames);
nbins = 50; % number of bins for historgram plots

% check the length of the vector containing which parameters
if lwp > npars || max(wp) > npars || lwp > 2 || lwp < 1 || lwp~=lparnames
    error('Error... parameters not properly specified!');
end

if lwp == 1
    % make histogram plot
    figure()
    hist(post_samples(:,wp), nbins);
    xlabel(parnames{1},'fontsize',12);

    % calculate mean and std deviation
    par_mean = mean(post_samples(:,wp));
    par_std = std(post_samples(:,wp));
    fprintf('Parameter %s: Mean = %f, stddev = %f\n', ...
            parnames{1}, par_mean, par_std);

elseif lwp == 2
    % make 2-d histogram plot
    p1 = wp(1);
    p2 = wp(2);
    edges1 = linspace(min(post_samples(:,p1)), max(post_samples(:,p1)), nbins);
    edges2 = linspace(min(post_samples(:,p2)), max(post_samples(:,p2)), nbins);
    histmat = hist2(post_samples(:,p1), post_samples(:,p2), edges1, edges2);

    figure()
    imagesc(edges1, edges2, transpose(histmat));
    set(gca,'YDir','normal')
    xlabel(parnames{1},'fontsize',12);
    ylabel(parnames{2},'fontsize',12);

end

return

