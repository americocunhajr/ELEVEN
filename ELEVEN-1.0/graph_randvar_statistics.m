
% -----------------------------------------------------------------
%  graph_randvar_statistics.m
%
%  This functions plots in the same figure the PDF curve and
%  other statistics of a given random variable.
%
%  input:
%  X_supp - pdf x data vector
%  X_ksd  - pdf y data vector
%  X_mean - mean value
%  X_std  - standard deviation
%  X_low  - lower quantile
%  X_upp  - upper quantile
%  gtitle - graph title
%  xlab   - x axis label
%  ylab   - y axis label
%  leg1   - legend 1
%  leg2   - legend 2
%  leg3   - legend 3
%  leg4   - legend 4
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
%  last update: Dec 26, 2017
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_randvar_statistics(X_supp,X_ksd,X_mean,...
                                        X_std,X_low,X_upp,...
                                        bins1,freq1,...
                                        gtitle,xlab,ylab,...
                                        leg1,leg2,leg3,leg4,...
                                        xmin,xmax,ymin,ymax,gname,flag)
	
    % check number of arguments
    if nargin < 20
        error('Too few inputs.')
    elseif nargin > 21
        error('Too many inputs.')
    elseif nargin == 20
        flag = 'none';
    end

    % check arguments
    if length(X_supp) ~= length(X_ksd)
        error('X_supp and X_ksd vectors must be same length')
    end
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh7 = bar(bins1,freq1,0.8);
    hold all
    fh1 = plot(X_supp,X_ksd,'-b','DisplayName',leg1);
    fh2 =  line([X_mean X_mean],[ymin ymax],'DisplayName',leg2);
    fh3 =  line([X_mean-X_std X_mean-X_std],[ymin ymax],'DisplayName',leg3);
    fh4 =  line([X_low X_low],[ymin ymax],'DisplayName',leg4);
    fh5 =  line([X_mean+X_std X_mean+X_std],[ymin ymax]);
    fh6 =  line([X_upp X_upp],[ymin ymax]);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    %leg = legend(leg1,leg2,leg3,leg4,'Location','northwest');
    leg = legend([fh1; fh2; fh3; fh4],'Location','northwest');
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    
    set(fh1,'Color','b');
    set(fh1,'LineStyle','-');
    set(fh1,'LineWidth',2.0);
    set(fh2,'Color','r');
    set(fh2,'LineStyle','-');
    set(fh2,'LineWidth',3.0);
    set(fh3,'Color','m');
    set(fh3,'LineStyle','-.');
    set(fh3,'LineWidth',2.0);
    set(fh4,'Color','g');
    set(fh4,'LineStyle','--');
    set(fh4,'LineWidth',2.0);
    set(fh5,'Color','m');
    set(fh5,'LineStyle','-.');
    set(fh5,'LineWidth',2.0);
    set(fh6,'Color','g');
    set(fh6,'LineStyle','--');
    set(fh6,'LineWidth',2.0);
    set(fh7,'FaceColor','none');
    set(fh7,'EdgeColor','k');
    set(fh7,'LineStyle','-');
    labX = xlabel(xlab,'FontSize',20,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',20,'FontName','Helvetica');
    %set(Xlab,'interpreter','latex');
    %set(Ylab,'interpreter','latex');
%     uistack(fh6,'top');
% 	uistack(fh5,'top');
%     uistack(fh4,'top');
% 	uistack(fh3,'top');
%     uistack(fh2,'top');
%     uistack(fh1,'top');
    
    hold off
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    %if ( strcmp(flag,'eps') )
        %saveas(gcf,gname,'png');
        saveas(gcf,gname,'epsc2');
    %    gname = [gname, '.eps'];
    %end

return
% -----------------------------------------------------------------
