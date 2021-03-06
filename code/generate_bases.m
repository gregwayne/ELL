function [celltypes,bfns] = generate_bases(type,pth,varargin)
%function [celltypes,bfns] = generate_bases(type,varargin)
% let's make simplified basis functions and see how changing them affect
% the kernel!
% basis types:
% raw
% nopad (raw but don't fill in stretches with no data)
% thr (pass threshold value)
% 2thr (pass thresholds: PCA, everyone else)

% I am bad at parameter-fitting, so as a first go let's just make everyone
% binary. Next I'll try accounting for cell response variability (still
% assume NC's are zero, because they probably are!)

rcell=loader_means(pth);
ncells=length(rcell);
tmax=4500;
bfns=zeros(ncells,tmax);
celltypes=rcell(:,1);

for i=1:ncells
    celltime=length(rcell{i,2});
    offset=0;
    temp= [rcell{i,2}(offset+1:min(tmax+offset,celltime)); rcell{i,2}(end)*ones(max(tmax-celltime+offset,0),1)];

    switch type
        case 'raw'
            bfns(i,:)=temp-mean(temp);
        case 'nopad'
            bfns(i,:)=nan(size(temp));
            bfns(i,1:min(tmax+offset,celltime)) = rcell{i,2}(offset+1:min(tmax+offset,celltime));
        case 'thr'
            vmin=min(rcell{i,2});
            vmax=max(rcell{i,2});
            thrval=.5; %safety default value
            if(~length(varargin))
                if(i==1)
                    disp('warning: no threshold value provided; using 0.5 default');
                end
            else
                thrval = varargin{1};
            end
            thr=(vmax-vmin)*thrval + vmin;
            bfns(i,:) = double((temp>thr));
        case '2thr'
            vmin=min(rcell{i,2});
            vmax=max(rcell{i,2});
            thrval=.5; %safety default value
            if(~length(varargin))
                if(i==1)
                    disp('warning: no threshold value provided; using 0.5 default');
                end
            else
                if(~isempty(strfind(rcell{i,1},'PCA')))
                    thrval = varargin{1};
                else
                    thrval = varargin{2};
                end
            end
            thr=(vmax-vmin)*thrval + vmin;
            bfns(i,:) = double(temp>thr);
        otherwise
            disp 'error: define your basis type!';
    end
end

