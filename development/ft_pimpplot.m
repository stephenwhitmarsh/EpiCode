function varargout = ft_pimpplot(hfig,cmap,pimp,set_printing_properties)
%--------------------------------------------------------------------------
% FT_PIMPPLOT can pimp single channel time frequency representations and
% topographical representations created by FieldTrip. Basically, it
% replots the plots created by imagesc and surface using contourf.
%
% See also FT_PLOT
%
% This file is a FieldTrip fork as part of the FieldTripWrapper toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% check input
if nargin < 4,      set_printing_properties = true; end 
if nargin < 3,      pimp = true;    end
if isempty(pimp),	pimp = true;	end
if nargin < 2,      cmap = @jet;    end
if nargin < 1,      hfig = gcf;     end

% set output
varargout = {};
if nargout > 0
    varargout = {hfig};
end
    
% % return if requested
% if isfalse(pimp)
%     return
% end

% specify number of steps in colormap
stps = 2^9;

% create colormap if needed
if isa(cmap,' function_handle')
    cmap = cmap(stps);
end
    

% get images axes handles
hs = findobj(hfig,'Type','image','-and','Tag','cip');
ht = findobj(hfig,'Type','surface');
ha = findobj(hfig,'Type','axes');

% pimp the single channel time-frequency representation
for i = 1:length(hs)
    
    % get title and data
    %title_str = get(get(get(hs(i),'Parent'),'Title'),'String');
    X = get(hs(i),'XData');
    Y = get(hs(i),'YData');
    Z = get(hs(i),'CData');
    A = get(hs(i),'AlphaData');
    %c = caxis(get(hs(i),'Parent'));
    
    % expand matrix. This is necessary because the original plot was made
    % using image
    x_step = mode(abs(diff(X)));
    y_step = mode(abs(diff(Y)));
    X = [X(1)-x_step X X(end)+x_step];
    Y = [Y(1)-y_step Y Y(end)+y_step];
    Z = [Z(:,1) Z Z(:,end)];
    Z = [Z(1,:); Z; Z(end,:)];
    A = [A(:,1) A A(:,end)];
    A = [A(1,:); A; A(end,:)];
    
    % create new figure if requested
    if pimp < 0
        if ~exist('hnfig','var')
            hnfig = copyobj(hfig,0);
        end
        figure(hfig);
    end
    
    % re-plot data using contourf
    set(hfig,'CurrentAxes',get(hs(i),'Parent'));
    delete(hs(i));
    hold on;
    [C,hn] = contourf(X,Y,Z,stps,'LineStyle','none');
    colormap(cmap);
    set(gca,'Box','off');
    axis([min(X)+x_step/2 max(X)-x_step/2 min(Y)+y_step/2 max(Y)-y_step/2]);
    % axis tight;
    
    % add the alpha-map from the original image
    if ~all(A(:))
        whitemat = ones([size(A) 3]);
        A(A>=1) = 0;    % A(A> 0 & A<=0.5) = 0.75;
        image(X,Y,whitemat,'AlphaData',A);
    end
    
    % restack the objects
    hl = findobj(hfig,'Type','line');
    uistack(hl,'top');
    %uistack(hn,'bottom');
    
    % % add line at time = 0
    % hold on; line([0 0],[min(Y) max(Y)],'Color','k','LineWidth',2);
    hold off;
    
end

% pimp the topographical plot
for i = 1:length(ht)
    
    % get data
    X = get(ht(i),'XData');
    Y = get(ht(i),'YData');
    Z = get(ht(i),'CData');
    
    % adjust X and Y to have centre at zero. This is necessary because the
    % original plot was made using surface
    X = X + 0.5/size(X,1);
    Y = Y + 0.5/size(Y,1);
    
    % create new figure if requested
    if pimp < 0
        if ~exist('hnfig','var')
            hnfig = copyobj(hfig,0);
        end
        figure(hfig);
    end
        
    % re-plot data using contourf
    set(hfig,'CurrentAxes',get(ht(i),'Parent'));
    delete(ht(i));
    hold on; %hold(ha(i),'on');
    [C,hn] = contourf(X,Y,Z,stps,'LineStyle','none');
    colormap(cmap);
    
    % restack the objects
    uistack(hn,'bottom');
        
end

% % continue working on new figure if requested
% if pimp < 0
%     hfig = hnfig;
% end
    
% % make background white
% set(hfig,'Color','w');

if set_printing_properties
    % set correct printing properties
    set(hfig,'PaperType','<custom>');
    %set(h,'PaperUnits','centimeters');
    if strcmpi(get(hfig,'PaperUnits'),'centimeters')
        set(hfig,'PaperSize',[20 15]);
        set(hfig,'PaperPosition',[0 0 20 15]);
    elseif strcmpi(get(hfig,'PaperUnits'),'inches')
        set(hfig,'PaperSize',[8 6]);
        set(hfig,'PaperPosition',[0 0 8 6]);
    end
end

% set output
if nargout > 0
    varargout = {hfig};
end
