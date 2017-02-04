function [S,wl] = processOOIspectra(files, iron_window, tau)
%function [S,wl] = processOOIspectra(files, iron_window, tau)
%
% This function loads an pre-processes the spectra from Ocean Optics
% spectrometer. Preprocession is 
% 1) Correction of the wavelength so that it has equal step
% 2) spectral averaging (ironing) and
% 3) correction for different exposure tines.
%
% Input parameters: 
% files       - files to extract spectra from in the same format as for
%               loadOOIspectra function (type help loadOOIspectra for 
%               more help)
% tau         - (optional) array of exposure times for each spectrum. It should have the
%               length meeting number of files. Assumed to be equal for all the
%               spectra if omitted.
% iron_window - (optional) width of window for 'ironing'. Assumed to be 0
%               (no averaging) if omitted. 
% 
% Ouput:
% 
% S - matrix containing spectra as columns
% wl - Corrected wavelength with equal step
load OOIsens;
OOIsens = ones(2048,1);
if nargin == 0
    [S0, wavelen0] =loadOOIspectra;
else
    [S0,wavelen0] = loadOOIspectra(files);
end;
if (isempty(S0))
    error('No files found!');
end;
if nargin < 2 
    iron_window = 0;
end;

wl = min(wavelen0) : (max(wavelen0)-min(wavelen0))/(length(wavelen0)-1) : max(wavelen0);
S0 = S0 ./ (OOIsens*ones(1,size(S0,2)));
S = interp1(wavelen0, S0, wl,'spline');
if (iron_window ~= 0)
    S = iron_gauss(interp1(wavelen0, S0, wl,'spline'), iron_window);
end;

if nargin > 2
    tau = tau(:)';
    if length(tau) ~= size(S,2)
        warning('Size of tau array mismatches number uf spectra, ignoring tau');
    else
        S = S./(ones(size(S,1),1)*tau);
    end;
end;
