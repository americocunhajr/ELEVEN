
% -----------------------------------------------------------------
%  graph_contourf2_pnt.m
%
%  This functions plots the contour map of the scalar functions
%  F: R^2 -> R and G: R^2 -> R
%
%  input:
%  x      - x mesh vector
%  y      - y mesh vector
%  F      - scalar field 1
%  G      - scalar field 2
%  gtitle - graph title
%  xlab   - x axis label
%  ylab   - y axis label
%  xmin   - x axis minimum value
%  xmax   - x axis maximum value
%  ymin   - y axis minimum value
%  ymax   - y axis maximum value
%  gname  - graph name
%  flag   - output file format (optional)
%
%  output:
%  gname.eps - output file in eps format (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 2, 2018
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_contourf2_pnt(x,y,F,G,Xmax,Ymax,gtitle,xlab,ylab,...
                                  xmin,xmax,ymin,ymax,gname,flag)
                                    
    % check number of arguments
    if nargin < 14
        error('Too few inputs.')
    elseif nargin > 15
        error('Too many inputs.')
    elseif nargin == 14
        flag = 'none';
    end
    
    % generate 2D mesh grid
    %xi      = linspace(xmin,xmax,3*length(x));
    %yi      = linspace(ymin,ymax,3*length(y));
    %[XI,YI] = meshgrid(xi,yi);
    %FI      = griddata(x,y,F,XI,YI,'nearest'); % discont
    %FI      = griddata(x,y,F,XI,YI,'linear');    % C0
    %FI      = griddata(x,y,F,XI,YI,'natural'); % C1
    %FI      = griddata(x,y,F,XI,YI,'cubic');   % C2
    %FI      = griddata(x,y,F,XI,YI,'v4');      % C2

    fig = figure('Name',gname,'NumberTitle','off');
    
    [XI,YI] = meshgrid(x,y);
    fh = contourf(XI,YI,F);
    
    
    %fh = contourf(XI,YI,FI);
    %colormap jet;
    %colormap cool
    colormap pink
    %colormap hot
    colorbar
    hold on
    fh1 = contour(XI,YI,G,'--');
    fh2 = plot(Xmax,Ymax,'xm','LineWidth',5);
    hold off
    set(gcf,'color','white');
    set(gca,'position',[0.15 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    set(fh2,'LineWidth',5.0);
    set(fh2,'MarkerSize',20.0);
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    if ( strcmp(xmin,'auto') || strcmp(xmax,'auto') )
        xlim('auto');
    else
        xlim([xmin xmax]);
    end
    if ( strcmp(ymin,'auto') || strcmp(ymax,'auto') )
        ylim('auto');
    else
        ylim([ymin ymax]);
    end
    labX = xlabel(xlab,'FontSize',18,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',18,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    

    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        gname = [gname, '.eps'];
        %graph_fixPSlinestyle(gname,gname);
    end

return
% -----------------------------------------------------------------