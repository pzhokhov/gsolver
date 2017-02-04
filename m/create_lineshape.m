function [S, lambda] = create_lineshape(mask)

[spec, wavelen] = loadSPEfiles(mask);
[maxspec, maxspecindex] = max(spec);

[lambda, i] = sort(wavelen(maxspecindex));
S = maxspec(i);