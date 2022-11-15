function [pdata,varargout] = xyz2lines(xyz,lnw,varargin)
% XYZ2LINES - Extract elevation transects from scattered data
% 
% PDATA = XYZ2LINES(XYZ,LNW) - extracts points from a scattered
%   dataset, XYZ (m x 3), along transects specified in a structure,
%   LNW. The structure of LNW is derived from a Hypack linefile (.lnw)
%   and can be read using READLNW. Extracted profiles are supplied in
%   a structure array, PDATA, with the following fields:
%
%       'line_number' - Transect number in LNW
%       'x'           - Easting
%       'y'           - Northing
%       'z'           - Elevation
%       'dist'        - Distance along transect
%       'offset'      - Distance from transect
%
% [PDATA,GDATA] = XYZ2LINES(XYZ,LNW) - points that did not fall on a 
%   transect are included in the optional output, GDATA. 
%   
%   For each transect in LNW, XYZ2LINES finds all points that fall close 
%   to the line (default = 5 m).  The line is then broken into segments 
%   (default = 1 m) and the closest point to the line in that segment is 
%   selected. Transects that do not contain a minimum number of points 
%   are discarded (default = 20 points). If there is a large gap along the
%   transect, a NaN is inserted.  Optional input parameters control 
%   the search radius, segment length, minimum number of points in a 
%   transect, and maximum number of points in each segment. 
%
% OPTIONAL INPUT
%       'max_offset'   - maximum search radius from the specified line
%       'max_gap'      - inserts a NaN for gaps larger than this value
%       'interval'     - window length (survey units)
%       'min_points'   - minimum number of points to be considered a transect
%       'extend_lines' - extend search beyond line endpoints. Extends both
%                        start and end of lines equally.
%       'points_per_seg' - maximum number of points in each segment. The
%                       closest n points will be used. 
%       'clipxy'       - specify a polygon (m x 2 matrix) to include in
%                        transects.  Points outside of polygon will not be
%                        included in profiles
%       
% SEE ALSO READLNW


%process inputs
p=inputParser;
opts={'max_offset',   3,    {'numeric'}, {'nonnegative'};...
    'max_gap',        5,    {'numeric'}, {'scalar';'nonnegative'};...
    'interval',       1,    {'numeric'}, {'nonnegative'};...
    'min_points',     20,   {'numeric'}, {'scalar'};...
    'extend_lines',   0,    {'numeric'},  {'scalar'};...
    'points_per_seg', 1,    {'numeric'},  {'scalar'};...
    'clipxy',         [],   {'numeric'},  {'ncols',2}};

addRequired(p,'xyz',@(x)validateattributes(x,{'numeric'},{'ncols',3}));
addRequired(p,'lnw',@(x)validateattributes(x,{'struct'},{'nonempty'}));
%optional input
cellfun(@(x)(p.addOptional(x{1},x{2},...
    @(y)(validateattributes(y, x{3},x{4})))),...
    num2cell(opts,2));
parse(p,xyz,lnw,varargin{:})
opt=p.Results;

point_id=(1:size(xyz,1))';

pdata=repmat(struct(...
    'x',[],'y',[],'z',[],...
    'dist',[],'offset',[],...
    'idx',[]),...
    length(lnw),1);
pidx=cell(length(lnw),1);
lidx=ones(length(lnw),1);

%find points in polygon if specified
if ~isempty(opt.clipxy)
    in=inpolygon(xyz(:,1),xyz(:,2),opt.clipxy(:,1),opt.clipxy(:,2));
else
    in=true(size(xyz,1),1);
end


for i = 1:length(lnw)
    if numel(lnw(i).x)==2
        
        m=(lnw(i).y(end)-lnw(i).y(1))/...
            (lnw(i).x(end)-lnw(i).x(1));
        theta=atan(m);
        xt = xyz(:,1).*cos(theta) + xyz(:,2).*sin(theta);
        yt = -xyz(:,1).*sin(theta) + xyz(:,2).*cos(theta);
        
        xl=lnw(i).x.*cos(theta) + ...
            lnw(i).y.*sin(theta);
        yl=-lnw(i).x.*sin(theta) + ...
            lnw(i).y.*cos(theta);
        
        
        %calculate along line distance
        xlt=sort(xl,'descend');
        dist=xt-xlt(1);
        
        %determine if the points are along the the line
        if ~isinf(opt.extend_lines)
            lrange=[0 diff(xlt)];
            xrange=[min(lrange)-opt.extend_lines,...
                max(lrange)+opt.extend_lines];
            fun=@(x,y)(x>=y(1) & x<=y(2));
            didx=fun(dist,xrange);
        end
            
            
        
        %determine distance from line to points
        offset=yt-yl(1);
        
        
       
        oidx=(abs(offset)<opt.max_offset & didx & in);
        
        xi=[min(dist(oidx)):opt.interval:max(dist(oidx)),max(dist(oidx))];
        [n,bin]=histc(dist(oidx),xi);
        
        obin=accumarray(bin,offset(oidx),[numel(n) 1],@(x){x},{NaN});
        ibin=accumarray(bin,point_id(oidx),[numel(n) 1],@(x){x},{NaN});
        
        if opt.points_per_seg==1
            [~,bidx]=cellfun(@(x)(min(abs(x))),obin,'un',0);
            tidx=cellfun(@(x,y)(x(y)),ibin,bidx);
            tidx(isnan(tidx))=[];
        else
            [~,bidx]=cellfun(@(x)(sort(abs(x))),obin,'un',0);
            
            if isfinite(opt.points_per_seg)
                np=cellfun(@numel,bidx);
                nflag=np>opt.points_per_seg;
                nidx=cellfun(@(x)(x(1:opt.points_per_seg)),...
                    bidx(nflag),'un',0);
                bidx(nflag)=nidx;
            end
            tidx=cell2mat(cellfun(@(x,y)(x(y)),ibin,bidx,'un',0));
            tidx(isnan(tidx))=[];
        end
        
        if numel(tidx)>opt.min_points
            [dsort,di]=sort(dist(tidx));
            pdata(i).idx=tidx(di);
            pdata(i).x=xyz(tidx(di),1);
            pdata(i).y=xyz(tidx(di),2);
            pdata(i).z=xyz(tidx(di),3);
            pdata(i).dist=dsort;
            pdata(i).offset=offset(tidx(di));
            
            %insert nans between large gaps
            dd=abs(diff(pdata(i).dist));
            gaps=find(dd>opt.max_gap);
            if ~isempty(gaps)
                pdata(i)=structfun(@(x)(insertrows(x,NaN,gaps)),...
                    pdata(i),'un',0);
            end
            
            pidx{i}=tidx;
        else
            lidx(i)=0;
        end 
        
    end
end
    

        
        
pdata(~lidx)=[];

%insert the line number
lnum=arrayfun(@(x)(x.name),lnw,'un',0);
xline=arrayfun(@(x)(x.x),lnw,'un',0);
yline=arrayfun(@(x)(x.y),lnw,'un',0);
lnum(~lidx)=[];
xline(~lidx)=[];
yline(~lidx)=[];
[pdata(:).line_num]=deal(lnum{:});
[pdata(:).line_x]=deal(xline{:});
[pdata(:).line_y]=deal(yline{:});

if nargout>1
    pidx(~lidx)=[];
    
    gidx=ismember(point_id,cell2mat(pidx));
    gdata.x=xyz(~gidx,1);
    gdata.y=xyz(~gidx,2);
    gdata.z=xyz(~gidx,3);
    varargout={gdata}';
end

%%%----------------------------------------------------------------------
function [C,RA,RB] = insertrows(A,B,ind)
% INSERTROWS - Insert rows into a matrix at specific locations
%   C = INSERTROWS(A,B,IND) inserts the rows of matrix B into the matrix A at
%   the positions IND. Row k of matrix B will be inserted after position IND(k)
%   in the matrix A. If A is a N-by-X matrix and B is a M-by-X matrix, C will
%   be a (N+M)-by-X matrix. IND can contain non-integers.
%
%   If B is a 1-by-N matrix, B will be inserted for each insertion position
%   specified by IND. If IND is a single value, the whole matrix B will be
%   inserted at that position. If B is a single value, B is expanded to a row
%   vector. In all other cases, the number of elements in IND should be equal to
%   the number of rows in B, and the number of columns, planes etc should be the
%   same for both matrices A and B. 
%
%   Values of IND smaller than one will cause the corresponding rows to be
%   inserted in front of A. C = INSERTROWS(A,B) will simply append B to A.
%
%   If any of the inputs are empty, C will return A. If A is sparse, C will
%   be sparse as well. 
%
%   [C, RA, RB] = INSERTROWS(...) will return the row indices RA and RB for
%   which C corresponds to the rows of either A and B.
%
%   Examples:
%     % the size of A,B, and IND all match
%        C = insertrows(rand(5,2),zeros(2,2),[1.5 3]) 
%     % the row vector B is inserted twice
%        C = insertrows(ones(4,3),1:3,[1 Inf]) 
%     % matrix B is expanded to a row vector and inserted twice (as in 2)
%        C = insertrows(ones(5,3),999,[2 4])
%     % the whole matrix B is inserted once
%        C = insertrows(ones(5,3),zeros(2,3),2)
%     % additional output arguments
%        [c,ra,rb] = insertrows([1:4].',99,[0 3]) 
%        c.'     % -> [99 1 2 3 99 4] 
%        c(ra).' % -> [1 2 3 4] 
%        c(rb).' % -> [99 99] 
%
%   Using permute (or transpose) INSERTROWS can easily function to insert
%   columns, planes, etc:
%
%     % inserting columns, by using the transpose operator:
%        A = zeros(2,3) ; B = ones(2,4) ;
%        c = insertrows(A.', B.',[0 2 3 3]).'  % insert columns
%     % inserting other dimensions, by using permute:
%        A = ones(4,3,3) ; B = zeros(4,3,1) ; 
%        % set the dimension on which to operate in front
%        C = insertrows(permute(A,[3 1 2]), permute(B,[3 1 2]),1) ;
%        C = ipermute(C,[3 1 2]) 
%
%  See also HORZCAT, RESHAPE, CAT

% for Matlab R13
% version 2.0 (may 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% 1.0, feb 2006 - created
% 2.0, may 2008 - incorporated some improvements after being selected as
% "Pick of the Week" by Jiro Doke, and reviews by Tim Davis & Brett:
%  - horizontal concatenation when two arguments are provided
%  - added example of how to insert columns
%  - mention behavior of sparse inputs
%  - changed "if nargout" to "if nargout>1" so that additional outputs are
%    only calculated when requested for

error(nargchk(2,3,nargin)) ;

if nargin==2,
    % just horizontal concatenation, suggested by Tim Davis
    ind = size(A,1) ;
end

% shortcut when any of the inputs are empty
if isempty(B) || isempty(ind),    
    C = A ;     
    if nargout > 1,
        RA = 1:size(A,1) ;
        RB = [] ;
    end
    return
end

sa = size(A) ;

% match the sizes of A, B
if numel(B)==1,
    % B has a single argument, expand to match A
    sb = [1 sa(2:end)] ;
    B = repmat(B,sb) ;
else
    % otherwise check for dimension errors
    if ndims(A) ~= ndims(B),
        error('insertrows:DimensionMismatch', ...
            'Both input matrices should have the same number of dimensions.') ;
    end
    sb = size(B) ;
    if ~all(sa(2:end) == sb(2:end)),
        error('insertrows:DimensionMismatch', ...
            'Both input matrices should have the same number of columns (and planes, etc).') ;
    end
end

ind = ind(:) ; % make as row vector
ni = length(ind) ;

% match the sizes of B and IND
if ni ~= sb(1),
    if ni==1 && sb(1) > 1,
        % expand IND
        ind = repmat(ind,sb(1),1) ;
    elseif (ni > 1) && (sb(1)==1),
        % expand B
        B = repmat(B,ni,1) ;
    else
        error('insertrows:InputMismatch',...
            'The number of rows to insert should equal the number of insertion positions.') ;
    end
end

sb = size(B) ;

% the actual work
% 1. concatenate matrices
C = [A ; B] ;
% 2. sort the respective indices, the first output of sort is ignored (by
% giving it the same name as the second output, one avoids an extra 
% large variable in memory)
[abi,abi] = sort([[1:sa(1)].' ; ind(:)]) ;
% 3. reshuffle the large matrix
C = C(abi,:) ;
% 4. reshape as A for nd matrices (nd>2)
if ndims(A) > 2,
    sc = sa ;
    sc(1) = sc(1)+sb(1) ;
    C = reshape(C,sc) ;
end

if nargout > 1,
    % additional outputs required
    R = [zeros(sa(1),1) ; ones(sb(1),1)] ;
    R = R(abi) ;
    RA = find(R==0) ;
    RB = find(R==1) ;
end









        