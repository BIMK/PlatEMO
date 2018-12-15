function Draw(Data,varargin)
%Draw - Display data.
%
%   Draw(P) displays all the points in matrix P, where each row of P
%   indicates one point, and each column of P indicates one dimension.
%
%   Draw(P,Name,Value,...) displays the points with the specified style.
%
%   If each point contains 2 or 3 dimensions, the points will be displayed
%   in 2-D or 3-D space. And if the number of dimensions of each point is
%   more than 3, display the points by parallel coordinates, where each
%   value of abscissa indicates one dimension, and each value of ordinate
%   indicates the value of the corresponding dimension.
%
%   Example:
%       Draw(rand(100,3),'sk')
%       Draw(rand(100,10))

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(Data);
    
    %% The size of the figure
    set(gca,'Unit','pixels');
    if get(gca,'Position') <= [inf inf 400 300]
        Size = [3 5 .8 8];
    else
        Size = [6 8 2 13];
    end
    
    %% The styple of the figure
    if nargin < 2
        if M == 2
            varargin = {'ok','MarkerSize',Size(1),'Marker','o','Markerfacecolor',[.7 .7 .7],'Markeredgecolor',[.4 .4 .4]};
        elseif M == 3
            varargin = {'ok','MarkerSize',Size(2),'Marker','o','Markerfacecolor',[.7 .7 .7],'Markeredgecolor',[.4 .4 .4]};
        elseif M > 3
            varargin = {'Color',[.5 .5 .5],'LineWidth',Size(3)};
        end
    end
    
    %% Draw the figure
    set(gca,'NextPlot','add','Box','on','Fontname','Times New Roman','FontSize',Size(4));
    if M == 2
        plot(Data(:,1),Data(:,2),varargin{:});
        xlabel('\itf\rm_1'); ylabel('\itf\rm_2');
        set(gca,'XTickMode','auto','YTickMode','auto','View',[0 90]);
        axis tight;
    elseif M == 3
        plot3(Data(:,1),Data(:,2),Data(:,3),varargin{:});
        xlabel('\itf\rm_1'); ylabel('\itf\rm_2'); zlabel('\itf\rm_3');
        set(gca,'XTickMode','auto','YTickMode','auto','ZTickMode','auto','View',[135 30]);
        axis tight;
    elseif M > 3
        Label = repmat([0.99,2:M-1,M+0.01],N,1);
        Data(2:2:end,:)  = fliplr(Data(2:2:end,:));
        Label(2:2:end,:) = fliplr(Label(2:2:end,:));
        Data  = Data';
        Label = Label';
        plot(Label(:),Data(:),varargin{:});
        xlabel('Dimension No.'); ylabel('Value');
        set(gca,'XTick',1:ceil(M/10):M,'XLim',[1,M],'YTickMode','auto','View',[0 90]);
    end
end