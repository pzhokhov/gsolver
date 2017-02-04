function f = animate_we(X, maxn, imagefn, title_templ)

h =figure;
if (nargin==2)
    imagefn = @imagesc;
    title_templ = 'Picture %d    \\leftarrow previous, \\rightarrow next';
end;
if (nargin == 3)
    if (isa(imagefn,'function_handle'))
        title_templ = 'Picture %d    \\leftarrow previous, \\rightarrow next';
    else 
        title_templ = imagefn;
    end;
end;
               
imagefn(X(1));

set(h,'KeyPressFcn', @(hObject, eventdata)keypress_fcn(hObject, eventdata, imagefn, X, maxn, title_templ));
set(h,'Tag', 'animatefig_n1');
title(sprintf(title_templ, 1));

return;

function keypress_fcn(hObject, eventdata, imagefn, X, maxN, title_templ)

if strcmp(eventdata.Key,'leftarrow')
 t = get(hObject, 'Tag');
 n = str2double(t(13:end));
 if (n>1)
    n=n-1;
    animate3d_redraw(hObject, imagefn, X, n, title_templ);
 end;
end;
if strcmp(eventdata.Key,'rightarrow')
 t = get(hObject, 'Tag');
 n = str2double(t(13:end));
 if (n<maxN)
    n=n+1;
    animate3d_redraw(hObject, imagefn, X, n, title_templ);
 end;
end;
return;

function animate3d_redraw(hObject, imagefn, X, n, title_templ)
    figure(hObject);
%     xlim = get(gca, 'xlim');
%     ylim = get(gca, 'ylim');
    ch = get(hObject, 'children');
    xlim = zeros(length(ch),2); ylim=xlim;
    for i=1:length(ch);
     xlim(i,:) = get(ch(i),'xlim');
     ylim(i,:) = get(ch(i),'ylim');
    end;
     imagefn(X(n));
%      set(gca,'xlim', xlim);
%      set(gca,'ylim', ylim);
    ch_ = get(hObject, 'children'); 
    for i=1:length(ch);
     set(ch_(ch_==ch(i)),'xlim',xlim(i,:));
     set(ch_(ch_==ch(i)),'ylim',ylim(i,:));
    end;

    
    set(hObject, 'Tag', sprintf('animatefog_n%d',n));
    title(sprintf(title_templ, n));
return;