function [lnw,varargout] = readLNW(varargin)
%READLNW -read Hypack line file(.lnw)
%
%   lnw = readLNW reads the formatted text file that
%   HYPACK uses for planning lines.  The function
%   requires no inputs.
%
%   OUTPUT:
%       line data are read into a structure
%       with fields name (line number), x and y.
%
%   USAGE OPTIONS:
%       lnw=readLNW('filename',example.lnw)- reads
%           in specified file rather than prompting
%           the user.
%
%       lnw=readLNW('plotit')- plots lines with
%           labels.
%
%       [lnw, sumstats]=readLNW - outputs a structure
%           summarizing line lengths, total distance 
%           for planning purposes.
%
% A. Stevens 3/22/2007
% astevens@usgs.gov

%check for input file
fpath=[];
if nargin>=1
    ind=strcmpi(varargin,'filename');
    if any(ind)==1
        fpath=varargin{ind+1};
        if exist(fpath,'file')==0;
            erstr=['File does not exist in MATLAB path. ',...
                'Check your file and path and try again...'];
            error(erstr)
        end
    end
end

if isempty(fpath)==1;
    [filename, pathname] = ...
        uigetfile('*.lnw', 'Pick an LNW-file');
    fpath=[pathname,filename];
end

%read file
fid=fopen(fpath,'r');
l1=fgetl(fid);
[~,numlines]=strread(l1,'%s %n');

lind=1;
while lind<=numlines;
    ln=fgetl(fid);
    [nmr,numpoints]=strread(ln,'%s %n');
    if strcmpi(nmr,'LIN')==1;
        for i=1:numpoints;
            lp=fgetl(fid);
            [~,x(i),y(i)]=strread(lp,'%s %f %f');
        end
        lend=fgetl(fid);
        [~,lname]=strread(lend,'%s %s');
        [base,ext]=strtok(lname,'_');
        if numel(base{:})<3
            
            lname=cellstr(sprintf('%03.3s%s',base{:},ext{:}));
        end


        lnw(lind).name=char(lname);
        lnw(lind).x=x;
        lnw(lind).y=y;

        lind=lind+1;
        clear x y lname
    end
end
fclose(fid);

%make a simple plot if specified
if any(strcmpi(varargin,'plotit'))==1;
    figure
    hold on

    %make ~50 labels per plot
    int=ceil(numlines/50);

    for i=1:numlines;
        plot(lnw(i).x,lnw(i).y);
        if rem(i,int)==0
            text(lnw(i).x(1),lnw(i).y(1),...
                num2str(lnw(i).name),'fontsize',10)
        end
    end
    axis equal
end

if nargout>1

    stats.numlines=numlines;
    for i=1:numlines
        x1=lnw(i).x(1:end-1);
        x2=lnw(i).x(2:end);
        y1=lnw(i).y(1:end-1);
        y2=lnw(i).y(2:end);

        xd = x2 - x1;
        yd = y2 - y1;
        distr(i)=sum(sqrt(xd.*xd +yd.*yd))';
    end
    stats.totalLength=sum(distr);
    stats.lineLens=distr;
    varargout{1}=stats;
end






