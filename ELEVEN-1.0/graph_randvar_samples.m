
% -----------------------------------------------------------------
%  graph_pdf_statistics.m
%
%  This functions plots in the same figure the samples and
%  other statistics of a given random variable.
%
%  input:
%  data_indx    - random variable sample indices
%  data_samples - random variable sample values
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
function fig = graph_randvar_samples(data_indx,data_samples,X_mean,...
                                    X_std,X_low,X_upp,...
                                    gtitle,xlab,ylab,...
                                    leg1,leg2,leg3,leg4,...
                                    xmin,xmax,ymin,ymax,gname,flag)
	
    % check number of arguments
    if nargin < 18
        error('Too few inputs.')
    elseif nargin > 19
        error('Too many inputs.')
    elseif nargin == 18
        flag = 'none';
    end

    % check arguments
    if length(data_indx) ~= length(data_samples)
        error('data_indx and data_samples vectors must be same length')
    end
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh1 = plot(data_indx,data_samples,'xb');
    hold all
    fh2 =  line([xmin xmax],[X_mean X_mean]);
    fh3 =  line([xmin xmax],[X_mean-X_std X_mean-X_std]);
    fh4 =  line([xmin xmax],[X_low X_low]);
    fh5 =  line([xmin xmax],[X_mean+X_std X_mean+X_std]);
    fh6 =  line([xmin xmax],[X_upp X_upp]);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    leg = legend(leg1,leg2,leg3,leg4,'Location','northwest');
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    
    set(fh1,'Color','b');
    %set(fh1,'LineStyle','-');
    set(fh1,'LineWidth',2.0);
    set(fh2,'Color','r');
    set(fh2,'LineStyle','-');
    set(fh2,'LineWidth',2.0);
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
    labX = xlabel(xlab,'FontSize',20,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',20,'FontName','Helvetica');
    %set(Xlab,'interpreter','latex');
    %set(Ylab,'interpreter','latex');
    
    hold off
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        %gname = [gname, '.eps'];
    end

return
% -----------------------------------------------------------------
