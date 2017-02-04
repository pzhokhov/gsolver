function [S, lambda] = creaate_lineshape(mask)

[spec, wavelen] = loadSPEfiles(mask);
[maxspec, maxspecindex] = max(spec);

[lambda, i] = sort(wavelen(macspecindex));
S = maxspec(i);