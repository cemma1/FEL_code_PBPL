function enhance_plot(fontname,fontsize,linewidth)

%  enhance_plot([fontname,fontsize,linewidth]);
%
%  Function to enhance MATLAB's lousy text choices on plots.  Sets the
%  current figure's Xlabel, Ylabel, Title, and all Text on plots, plus
%  the axes-labels to the "fontname" and "fontsize" input here where
%  the defaults have been set to 'times' and 16.
%
%  INPUTS:  fontname:   (Optional,DEF='TIMES') FontName string to use
%                       MATLAB's ugly default is 'Helvetica'
%           fontsize:   (Optional,DEF=16) FontSize integer to use
%                       MATLAB's tiny default is 10
%           linewidth:  (Optional,DEF=1.0) LineWidth value to use
%                       MATLAB's spidery default is 0.5
%======================================================================

if ~exist('fontname')
  fontname='times';
end
if ~exist('fontsize')
  fontsize=16;
end
if ~exist('linewidth')
  linewidth=1.0;
end

Hf=gcf;

% loop over children of figure ... this gets the legends, too

Hnf=get(Hf,'Children');
for n=1:length(Hnf)
  Ha=Hnf(n);
  set(Ha,'FontName',fontname);
  set(Ha,'FontSize',fontsize);
  set(Ha,'LineWidth',linewidth);
%{
  Hx=get(Ha,'XLabel');
  set(Hx,'FontName',fontname);
  set(Hx,'FontSize',fontsize);
  set(Hx,'VerticalAlignment','cap');

  Hy=get(Ha,'YLabel');
  set(Hy,'FontName',fontname);
  set(Hy,'FontSize',fontsize);
  set(Hy,'VerticalAlignment','bottom');

  Hy=get(Ha,'ZLabel');
  set(Hy,'FontName',fontname);
  set(Hy,'FontSize',fontsize);
 %set(Hy,'VerticalAlignment','bottom');

  Ht=get(Ha,'Title');
  set(Ht,'FontName',fontname);
  set(Ht,'FontSize',fontsize);
  set(Ht,'VerticalAlignment','baseline');
%}
  Hn=get(Ha,'Children');
  for m=1:length(Hn)
    t=get(Hn(m),'Type');
    if (strcmp(t,'text'))
      set(Hn(m),'FontName',fontname);
      set(Hn(m),'FontSize',fontsize);
    elseif (strcmp(t,'line'))
      set(Hn(m),'LineWidth',linewidth);
    end
  end
end

figure(Hf);
