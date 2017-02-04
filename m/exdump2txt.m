function exdump2txt(filemask)

files = dir(filemask);
if isempty(files)
    error('No files match requred name mask.');
end;

for i = 1:length(files)
    filename = files(i).name;
    outname = [filename '.txt'];
    [A, xnet, znet] = load_exdump(filename);
    eqistepsave(A,xnet,znet, outname);
    disp(sprintf('%s --> %s ... Done', filename, outname));
end;