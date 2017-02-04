function [x,y] = read_plotdigitizer_xml(filename)

doc = xmlread(filename); 
P = doc.getElementsByTagName('point');

N = P.getLength; 
x = zeros(N,1); 
y = zeros(N,1);
for i=1:N; 
    cP = P.item(i-1); 
    cP_attr = cP.getAttributes;
    cx_text = char(cP_attr.item(0)); x(i)=str2num(cx_text(5:end-1));   
    cy_text = char(cP_attr.item(1)); y(i)=str2num(cy_text(5:end-1)); 
end

