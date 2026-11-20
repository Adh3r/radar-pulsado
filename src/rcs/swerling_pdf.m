function p = swerling_pdf(sigma, sigma_mean, model)
%SWERLING_PDF  Probability density of the RCS for the Swerling fluctuation models.
%   p = swerling_pdf(sigma, sigma_mean, model)
%     model 1 (or 2): many comparable scatterers -> exponential (chi-square, 2 dof)
%                     p = exp(-sigma/mean) / mean
%     model 3 (or 4): one dominant plus many small -> chi-square, 4 dof
%                     p = 4 sigma / mean^2 * exp(-2 sigma / mean)
%   Models 1/3 fluctuate scan to scan and 2/4 pulse to pulse; the pdf is the same.
switch model
    case {1, 2}
        p = exp(-sigma/sigma_mean) / sigma_mean;
    case {3, 4}
        p = 4*sigma/sigma_mean^2 .* exp(-2*sigma/sigma_mean);
    otherwise
        error('swerling_pdf:model', 'Model must be 1, 2, 3 or 4.');
end
end
