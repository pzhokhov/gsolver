function [dw, dwc, dwpa, E, jc, jpa] = load_absorbed_energy(Eg, w0, dt, N, sizes, name_temp, varargin)

const_SI; 

jc = readmatrices_text([name_temp '_jc.dat'], sizes, varargin);
jpa = readmatrices_text([name_temp '_jpa.dat'], sizes, varargin);
E = readmatrices_text([name_temp '_E.dat'], sizes, varargin);

if (sizes(1)==1)

dwc  = cumtrapz(squeeze(((jc) .*E))).*dt.*Eg.^2/SI.h_.*2*pi/w0.*N/1e9;
dwpa = cumtrapz(squeeze(((jpa).*E))).*dt.*Eg.^2/SI.h_.*2*pi/w0.*N/1e9;

else
    dwc  = cumtrapz(squeeze(sum((jc) .*E))).*dt.*Eg.^2/SI.h_.*2*pi/w0.*N/1e9;
    dwpa = cumtrapz(squeeze(sum((jpa).*E))).*dt.*Eg.^2/SI.h_.*2*pi/w0.*N/1e9;
end;

dw = dwc+dwpa;
