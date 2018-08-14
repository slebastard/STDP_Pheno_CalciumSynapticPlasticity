function [ax,h]=subtitle(text)
%
%Centers a title over a group of subplots.
%Returns a handle to the title and the handle to an axis.
% [ax,h]=subtitle(text)
%           returns handles to both the axis and the title.
% ax=subtitle(text)
%           returns a handle to the axis only.
ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title(text,'FontSize',15);
if (nargout < 2)
    return
end
% t = get(ax,'Title');
% p = get(t, 'Position');
% set(t,'pos',p+[0 0.04 0 0.04]);
h=get(ax,'Title');

