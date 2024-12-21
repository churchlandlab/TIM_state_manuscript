
function stdshade(amatrix,alpha,acolor,F,smth,naninds,varargin)
% usage: stdshade(amatrix,acolor,F,alpha,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - amatrix is a data matrix (observations,datalength)
% - acolor defines the used color (default is red)
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23
%
% Edit: Optionally added possibility to define amean and astd as additional
% input: stdshade(amatrix,acolor,F,alpha,smth,amean,astd)
%

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r';
end

if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1;
end

if ne(size(F,1),1)
    F=F';
end

if length(varargin)==2
    amean=smooth(varargin{1},smth)';
    astd=smooth(varargin{2},smth)';
    
else
    if size(amatrix,1) == 1 || size(amatrix,2) == 1
        amean(1,:) = amatrix;
        astd = zeros(1,length(amean));
        F = 1 : length(amean);
    else
%         amean=smooth(nanmean(amatrix,1),smth)';
%         astd=nanstd(amatrix)/sqrt(sum(~isnan(amatrix(:,1)))); % to get sem shading
%         astd(isnan(astd)) = 0;

        amean=smoothdata(nanmean(amatrix,1),'gaussian',smth);
        astd=nanstd(amatrix)/sqrt(sum(~isnan(amatrix(:,1)))); % to get sem shading
        astd(isnan(astd)) = 0;
        
    end
    %    astd=nanstd(amatrix); % to get std shading
end

lineMean = amean;
amean = fillmissing(amean, 'previous');
fillx = [F fliplr(F)];
filly = [amean+astd fliplr(amean-astd)];
ignoreinds = naninds;

x = cell(length(ignoreinds)+1,1);
y = cell(length(ignoreinds)+1,1);
%create separate boxes to fill around nans
for i = 1:length(ignoreinds)+1
    if i==1 %first box
        inds = 1:ignoreinds(1)-1;
        fillx = [F(inds) fliplr(F(inds))];
        filly = [amean(inds)+astd(inds) fliplr(amean(inds)-astd(inds))];
        x{i} = fillx;
        y{i} = filly;
    elseif i==length(ignoreinds)+1 %last box
        inds = ignoreinds(i-1)+1:length(F);
        fillx = [F(inds) fliplr(F(inds))];
        filly = [amean(inds)+astd(inds) fliplr(amean(inds)-astd(inds))];
        x{i} = fillx;
        y{i} = filly;
    else %all other boxes
        inds = ignoreinds(i-1)+1:ignoreinds(i)-1;
        fillx = [F(inds) fliplr(F(inds))];
        filly = [amean(inds)+astd(inds) fliplr(amean(inds)-astd(inds))];
        x{i} = fillx;
        y{i} = filly;
    end
end

hold on;
for i = 1:length(x) %iterate thru boxes
    if exist('alpha','var')==0 || isempty(alpha)
        fill(x{i},y{i},acolor,'linestyle','none');
        acolor='k';
    else
        fill(x{i},y{i},acolor, 'FaceAlpha', alpha,'linestyle','none');
    end
    
    if ishold==0
        check=true; else check=false;
    end
end
hold on;
lineMean(naninds) = NaN;
cLine = plot(F,lineMean,'color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end



