function hfig = orthogui(varargin)
% hfig = orthogui(voldata,varargin)

% KWJ 2013

%%%%%% can call orthogui(fig,[x y z]) to move cursor externally
if(numel(varargin) > 1 && numel(varargin{1}) == 1 && ishghandle(varargin{1}) && ~isempty(regexp(get(varargin{1},'tag'),'^orthogui')))
    hfig = varargin{1};
    switch lower(varargin{2})
        case {'getlocation','getloc'}
            M = getappdata(hfig,'guidata');
            hfig =  [M.cx M.cy M.cz];
        case {'setdata'}
            M = getappdata(hfig,'guidata');
            D = getappdata(hfig,'data');
            D.V = varargin{3};
            setappdata(hfig,'data',D);
            updatefig(M,false);
        case {'setalphadata'}
            M = getappdata(hfig,'guidata');
            D = getappdata(hfig,'data');
            D.Valpha = varargin{3};
            setappdata(hfig,'data',D);
            updatefig(M,false);
        case {'location','loc'}
            [x, y, z] = splitvars(varargin{3});
            
            M = getappdata(hfig,'guidata');
            M.cx = x;
            M.cy = y;
            M.cz = z;
            setappdata(hfig,'guidata',M);
            updatefig(M);
        case {'location_nocb','loc_nocb'}
            [x, y, z] = splitvars(varargin{3});
            
            M = getappdata(hfig,'guidata');
            M.cx = x;
            M.cy = y;
            M.cz = z;
            setappdata(hfig,'guidata',M);
            updatefig(M,false);
        case {'callback'}
            M = getappdata(hfig,'guidata');
            cb = varargin{3};
            if(~iscell(cb))
                cb = {cb};
            end
            M.cbfunc = cb;
            setappdata(hfig,'guidata',M);
        case {'addcallback'}
            M = getappdata(hfig,'guidata');
            cb = varargin{3};
            if(iscell(cb))
                M.cbfunc = [M.cbfunc cb];
            else
                M.cbfunc = [M.cbfunc {cb}];
            end
            
            setappdata(hfig,'guidata',M);            
    end
   
    return;
    
%%%%%% main function: hfig = orthogui(V, parameters)
else
    V = double(varargin{1});
    varargin = varargin(2:end);
    
    %if(~iscell(V))
    %    V = {V};
    %end
end

%hack to avoid annoying mrVista function shadowing
if(regexpimatch(which('annotation'),'mrloadret'))
    d = justdir(which('annotation'));
    rmpath(d);
    addpath(d,'-end');
end

spmdir = justdir(which('spm'));
defaultstruct_colin = [spmdir '/canonical/single_subj_T1.nii'];
defaultstruct= [spmdir '/templates/T1.nii'];
%structimage = [spmdir '/templates/EPI.nii'];

p = inputParser;
p.addParamValue('location',[]);
p.addParamValue('loc',[]);
p.addParamValue('background',[]);
p.addParamValue('bg',[0]);
p.addParamValue('callback',{});
p.addParamValue('title','Orthoviews');
p.addParamValue('colormap','jet');
p.addParamValue('alpha',[]);
p.addParamValue('maxalpha',[1]);
p.addParamValue('mip',true);
p.addParamValue('dim',[]);
p.addParamValue('bgdim',[]);
p.addParamValue('link',[]);
p.addParamValue('clim',[]);

p.parse(varargin{:});
r = p.Results;
loc = r.location;
cbfunc = r.callback;
titlestr = r.title;
bgfile = r.background;
if(isempty(bgfile))
    bgfile = r.bg;
end
if(isempty(loc))
    loc = r.loc;
end
cmap = r.colormap;
alpha = r.alpha;
maxalpha = r.maxalpha;
showmip = r.mip;
dimpermute = r.dim;
bgpermute = r.bgdim;
figlink = r.link;
colormap_clim = r.clim;

if(~iscell(cbfunc))
    cbfunc = {cbfunc};
end

Vbg = [];

is_bg_default = false;
if(isempty(bgfile))
    is_bg_default = true;
    bgstruct = spm_vol(defaultstruct);
    [Vbg xyz_bg] = spm_read_vols(bgstruct);
elseif(ischar(bgfile))
    if(strcmpi(bgfile,'default') || strcmpi(bgfile,'def'))
        is_bg_default = true;
        bgfile = defaultstruct;
    elseif(strcmpi(bgfile,'colin'))
        is_bg_default = true;
        bgfile = defaultstruct_colin;
    end
    bgstruct = spm_vol(bgfile);
    [Vbg xyz_bg] = spm_read_vols(bgstruct);
elseif(isstruct(bgfile))
    [Vbg xyz_bg] = spm_read_vols(bgfile);
elseif(isnumeric(bgfile))
    if(numel(bgfile) == 1)
        Vbg = bgfile*ones(size(V));
    else
        Vbg = bgfile;
    end
else
end

if(isempty(bgpermute))
    bgpermute = dimpermute;
end


if(isempty(maxalpha) && numel(alpha) == 1)
    maxalpha = alpha;
end

if(numel(maxalpha) == 0)
    if(isempty(bgfile))
        maxalpha = 1;
    else
        maxalpha = .5;
    end
end

Valpha = [];
if(numel(alpha) > 1)
    Valpha = abs(alpha);
    Valpha = Valpha/max(Valpha(:));    
end

if(~isempty(dimpermute))
    %V = cellfun(@(x)(permute(x,dimpermute)),V,'uniformoutput',false);
    if(any(dimpermute < 0))
        df = find(dimpermute < 0);
        for i = 1:numel(df)
            V = flipdim(V,df(i));
            if(~isempty(Valpha))
                Valpha = flipdim(Valpha,df(i));
            end
        end
    end
    V = permute(V,abs(dimpermute));
    if(~isempty(Valpha))
        Valpha = permute(Valpha,abs(dimpermute));
    end
end

if(~isempty(bgpermute))
    if(~isempty(Vbg) && ~is_bg_default)
        
        if(any(bgpermute < 0))
            df = find(bgpermute < 0);
            for i = 1:numel(df)
                Vbg = flipdim(Vbg,df(i));
            end
        end
        Vbg = permute(Vbg,abs(bgpermute));
    end
end

do_updatelinked = false;
sz = size(V);
if(isempty(loc))
    loc = round(sz/2);
end
if(ischar(loc) && strcmpi(loc,'max'))
    do_updatelinked = true;
    [~,midx] = nanmax(abs(V(:)));
    [mx my mz] = ind2sub(size(V),midx);
    loc = [mx my mz];
end

[cx cy cz] = splitvars(loc);

if(ischar(cmap))
    cmap = eval(sprintf('%s(%d)',cmap,256));
end

if(isempty(figlink))
    figtag = 'orthogui';
else
    figtag = sprintf('orthogui.%d',figlink);
end
hfig = figure('name',titlestr,'NumberTitle','off','WindowButtonMotionFcn',@fig_mousemove,...
    'WindowButtonUpFcn',@ax_mouseup,'tag',figtag,'WindowKeyPressFcn',@fig_keypress);

Vbg = Vbg./nanmax(Vbg(:));

if(isempty(bgfile))
    Vbg = 0*Vbg;
end

[mipx mipy mipz] = splitvars(round(size(V)/2));
[mipbg1 mipbg2 mipbg3] = bgslice(Vbg,[mipx mipy mipz],sz);
%if(~isempty(bgfile))
    [bg1 bg2 bg3] = bgslice(Vbg,[cx cy cz],sz);
%end


cursorgap = 1;

ax = [];
img = [];
bgimg = [];
hcurH = {};
hcurV = {};
hcurVmip = {};
hcurHmip = {};
mipsize = [];
mipidx = [];
imgmipbg = [];
imgmip = [];
axmip = [];

ax(1) = axes('position',[0 .5 .5 .5]);
%if(~isempty(bgfile))
    bgimg(1) = image([1 sz(1)],[1 sz(3)],bg1);
    hold on;
%end

img(1) = imagesc(squeeze(V(:,cy,:)).');
axis equal xy;
colormap(cmap);
hold on;
%[hcurV(1) hcurH(1)] = splitvars(plot([cx cx; 0 sz(3)]',[0 sz(1); cz cz]','w'));
[hcurV{1} hcurH{1}] = plotcursor(ax(1),cx,cz,[0 sz(1)],[0 sz(3)],cursorgap,'w');

ax(2) = axes('position',[.5 .5 .5 .5]);
%if(~isempty(bgfile))
    bgimg(2) = image([1 sz(2)],[1 sz(3)],bg2);
    hold on;
%end
img(2) = imagesc(squeeze(V(cx,:,:)).');
axis equal xy;
colormap(cmap);
hold on;
%[hcurV(2) hcurH(2)] = splitvars(plot([cy cy; 0 sz(3)]',[0 sz(2); cz cz]','w'));
[hcurV{2} hcurH{2}] = plotcursor(ax(2),cy,cz,[0 sz(3)],[0 sz(2)],cursorgap,'w');

ax(3) = axes('position',[0 0 .5 .5]);
%if(~isempty(bgfile))
    bgimg(3) = image([1 sz(1)],[1 sz(2)],bg3);
    hold on;
%end
img(3) = imagesc(squeeze(V(:,:,cz)).');
axis equal xy;
colormap(cmap);
hold on;
%[hcurV(3) hcurH(3)] = splitvars(plot([cx cx; 0 sz(2)]',[0 sz(1); cy cy]','w'));
[hcurV{3} hcurH{3}] = plotcursor(ax(3),cx,cy,[0 sz(2)],[0 sz(1)],cursorgap,'w');

hpanel = uipanel('position',[.5 0 .5 .5]);
if(showmip)
    %axhist_pos = [.5 .25 .5 .2];
    %axmip_pos = [.5 0 .5 .25];
    axhist_pos = [0 .5 1 .45];
    axmip_pos = [0 0 1 .5]; 
    axmip = axes('parent',hpanel,'outerposition',axmip_pos);
    mipbgsz = max([size(mipbg1); size(mipbg2); size(mipbg3)],[],1);
    mipbg1 = padimageto(mipbg1,mipbgsz,[],[0 0 0]);
    mipbg2 = padimageto(mipbg2,mipbgsz,[],[0 0 0]);
    mipbg3 = padimageto(mipbg3,mipbgsz,[],[0 0 0]);
    mipbg = [mipbg1 mipbg2 mipbg3];


    %unsigned MIP (mip(abs))
    [umip1 umipidx1] = nanmax(abs(V),[],2);
    [umip2 umipidx2] = nanmax(abs(V),[],1);
    [umip3 umipidx3] = nanmax(abs(V),[],3);
    
    [mip1 mipidx1] = nanmax(V,[],2);
    [mip2 mipidx2] = nanmax(V,[],1);
    [mip3 mipidx3] = nanmax(V,[],3);
    
    %%%%%%
    mip1(mip1 < umip1) = -umip1(mip1 < umip1);
    mip2(mip2 < umip2) = -umip2(mip2 < umip2);
    mip3(mip3 < umip3) = -umip3(mip3 < umip3);
    
    mipidx1(mip1 < umip1) = umipidx1(mip1 < umip1);
    mipidx2(mip2 < umip2) = umipidx2(mip2 < umip2);
    mipidx3(mip3 < umip3) = umipidx3(mip3 < umip3);
    
    mip1 = squeeze(mip1).';
    mip2 = squeeze(mip2).';
    mip3 = squeeze(mip3).';
    mipidx1 = squeeze(mipidx1).';
    mipidx2 = squeeze(mipidx2).';
    mipidx3 = squeeze(mipidx3).';
    mipidx = {mipidx1, mipidx2, mipidx3};
    
    mipsize = [size(mip1); size(mip2); size(mip3)];
    padsz = max(mipsize,[],1);
    mip1 = padimageto(mip1,padsz,[],nan);
    mip2 = padimageto(mip2,padsz,[],nan);
    mip3 = padimageto(mip3,padsz,[],nan);
    mip = [mip1 mip2 mip3];
    
    imgmipbg = image([1 size(mip,2)],[1 size(mip,1)],mipbg);
    hold on;
    imgmip = imagesc(mip);
    colormap(cmap);
    set(imgmip,'alphadata',maxalpha.*(mip ~= 0 & ~isnan(mip)));
    axis equal xy tight;
    set(axmip,'xtick',[],'ytick',[],'color',[0 0 0]);
    
    [hcurVmip{1} hcurHmip{1}] = plotcursor(axmip,0,0,[0 1],[0 1],0,'w');
    [hcurVmip{2} hcurHmip{2}] = plotcursor(axmip,0,0,[0 1],[0 1],0,'w');
    [hcurVmip{3} hcurHmip{3}] = plotcursor(axmip,0,0,[0 1],[0 1],0,'w');
else
   
    %axhist_pos = [.5 .05 .5 .4];
    axhist_pos = [0 0 1 .9];
end

axhist = axes('parent',hpanel,'outerposition',axhist_pos,'color',[1 1 1]);
[hx xi] = hist(V(V(:) ~= 0 & ~isnan(V(:))),100);
hhistfill = fill(xi([1 end end 1]),[0 0 max(hx) max(hx)],'k','linestyle','none','facealpha',.1);
hold on;
hbar = bar(xi,hx,1,'k');


%set(axhist, 'yscale','log');
yl = get(gca,'ylim');

hcurhist = plot([0 0],yl,'r');
set(hhistfill,'ydata',yl([1 1 2 2]));
title('Intensity hist (click to find voxel)');

if(~isempty(colormap_clim) || all(isnan(V(:))))
    cl = colormap_clim;
else
    cl = [nanmin(V(:)) nanmax(V(:))];
    if(all(cl >= 0))
        cl = [0 max(cl)];
    else
        cl = max(abs(cl))*[-1 1];
    end
end
set(hfig,'color',[0 0 0]);
set(img,'alphadata',maxalpha);
set([hcurH{:} hcurV{:} hcurHmip{:} hcurVmip{:} img bgimg hbar hpanel imgmip imgmipbg hhistfill],'hittest','off');
set([ax axmip],'xtick',[],'ytick',[],'color',[0 0 0],'clim',cl);

%%%% make colorbar and fix its position
set(axhist,'clim',cl);
cb = colorbar('peer',axhist,'location','eastoutside','color',[0 0 0],'xcolor',[0 0 0],'ycolor',[0 0 0]);
pc = get(cb,'position');
ph = get(axhist,'position');
pc(2) = ph(2);
set(cb,'position',pc);
set(axhist,'position',ph);

for i = 1:numel(ax)
    set(ax(i),'ButtonDownFcn',{@ax_mousedown,i});
end

set(axhist,'ButtonDownFcn',@hist_mousedown);
set(axmip,'ButtonDownFcn',@axmip_mousedown);

curax = 1;
mouseax = [];
curax_rect = annotation(hfig,'rectangle',[0 0 1 1],'color',[1 0 0],'hittest','off');

D = fillstruct(V,Vbg,Valpha);
M = fillstruct(hfig,ax,img,bgimg,bgfile,maxalpha,cx,cy,cz,...
    hcurV,hcurH,hpanel,axmip,mouseax,curax,curax_rect,hcurVmip,hcurHmip,...
    imgmip,imgmipbg,axhist,hcurhist,titlestr,cbfunc,figlink,cursorgap,mipsize,mipidx);

setappdata(hfig,'guidata',M);
setappdata(hfig,'data',D);
if(do_updatelinked)
    update_linked_figures(M);
end
updatefig(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatefig(M, run_callback)
if(nargin == 1)
    run_callback = true;
end
dumpstruct(M);
D = getappdata(hfig,'data');
sz = size(D.V);
cx = min(sz(1),max(1,round(cx)));
cy = min(sz(2),max(1,round(cy)));
cz = min(sz(3),max(1,round(cz)));


set(hfig,'name',sprintf('%s [%d %d %d] %f',titlestr,[cx cy cz],D.V(cx,cy,cz)));
set(hcurhist,'xdata',D.V(cx,cy,cz)*[1 1]);
s1 = squeeze(D.V(:,cy,:)).';
s2 = squeeze(D.V(cx,:,:)).';
s3 = squeeze(D.V(:,:,cz)).';
mask1 = (s1 ~= 0 & ~isnan(s1));
mask2 = (s2 ~= 0 & ~isnan(s2));
mask3 = (s3 ~= 0 & ~isnan(s3));
if(~isempty(D.Valpha))
    a1 = squeeze(D.Valpha(:,cy,:)).';
    a2 = squeeze(D.Valpha(cx,:,:)).';
    a3 = squeeze(D.Valpha(:,:,cz)).';    
    mask1 = mask1.*a1;
    mask2 = mask2.*a2;
    mask3 = mask3.*a3;
end
set(img(1),'CData',s1,'alphadata',maxalpha.*mask1);
set(img(2),'CData',s2,'alphadata',maxalpha.*mask2);
set(img(3),'CData',s3,'alphadata',maxalpha.*mask3);
if(~isempty(bgfile))
    [bg1 bg2 bg3] = bgslice(D.Vbg,[cx cy cz],sz);
    set(bgimg(1),'CData',bg1);
    set(bgimg(2),'CData',bg2);
    set(bgimg(3),'CData',bg3);
end

if(~isempty(curax))
    axrect = get(ax(curax),'position');
    axrect(1:2) = axrect(1:2)+.005*axrect(3:4);
    axrect(3:4) = .99*axrect(3:4);
    set(curax_rect,'units','normalized','position',axrect,'visible','on');
else
    set(curax_rect,'visible','off');
end

if(~isempty(axmip))
    [mipx mipy xl yl] = ax2mip(1,cx,cz,mipsize);
    updatecursor(hcurVmip{1},hcurHmip{1},mipx,mipy,xl,yl,0);

    [mipx mipy xl yl] = ax2mip(2,cy,cz,mipsize);
    updatecursor(hcurVmip{2},hcurHmip{2},mipx,mipy,xl,yl,0);

    [mipx mipy xl yl] = ax2mip(3,cx,cy,mipsize);
    updatecursor(hcurVmip{3},hcurHmip{3},mipx,mipy,xl,yl,0);
end

%plot(axmip,[1 size(s1,1)],[1 size(s1,2)],'w');
updatecursor(hcurV{1},hcurH{1},cx,cz,[0 sz(1)],[0 sz(3)],cursorgap);
% set(hcurV{1},'XData',[cx cx],'YData',[0 sz(1)]);
% set(hcurH{1},'XData',[0 sz(1)],'YData',[cz cz]);

updatecursor(hcurV{2},hcurH{2},cy,cz,[0 sz(2)],[0 sz(3)],cursorgap);
%set(hcurV{2},'XData',[cy cy],'YData',[0 sz(3)]);
%set(hcurH{2},'XData',[0 sz(2)],'YData',[cz cz]);

updatecursor(hcurV{3},hcurH{3},cx,cy,[0 sz(1)],[0 sz(2)],cursorgap);
%set(hcurV{3},'XData',[cx cx],'YData',[0 sz(2)]);
%set(hcurH{3},'XData',[0 sz(1)],'YData',[cy cy]);
if(run_callback && ~isempty(cbfunc))
    for i = 1:numel(cbfunc)
        if(isempty(cbfunc{i}))
            continue;
        end
        cbfunc{i}([cx cy cz]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hist_mousedown(src,event)
p = get(src,'CurrentPoint');
p = p(1);
M = getappdata(gcbf,'guidata');
D = getappdata(gcbf,'data');
idx1 = find(D.V(:) <= p);
if(isempty(idx1))
    [~,idx1] = min(D.V(:));
end
[~, idx2] = max(D.V(idx1));
idx = idx1(idx2);
[x y z] = ind2sub(size(D.V),idx);
M.cx = x;
M.cy = y;
M.cz = z;
update_linked_figures(M);
setappdata(gcbf,'guidata',M);
updatefig(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function axmip_mousedown(src,event)
M = getappdata(gcbf,'guidata');
M.mouseax = -1;
setappdata(gcbf,'guidata',M);
fig_mousemove(src,event);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ax_mousedown(src,event, idx)
rightclick = strcmpi(get(gcbf,'SelectionType'),'alt');
M = getappdata(gcbf,'guidata');

if(rightclick)
    D = getappdata(gcbf,'data');
    sz = size(D.V);
    cx = min(sz(1),max(1,round(M.cx)));
    cy = min(sz(2),max(1,round(M.cy)));
    cz = min(sz(3),max(1,round(M.cz)));
    
    
    fprintf('[%d %d %d] %f \n',[cx cy cz],D.V(cx,cy,cz));
    return;
end
M.mouseax = idx;
M.curax = idx;

setappdata(gcbf,'guidata',M);
fig_mousemove(src,event);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ax_mouseup(src,event)
M = getappdata(gcbf,'guidata');
M.mouseax = [];
setappdata(gcbf,'guidata',M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig_mousemove(src,event)
M = getappdata(gcbf,'guidata');
if(~isempty(M.mouseax))
    idx = M.mouseax;
    rightclick = strcmpi(get(gcbf,'SelectionType'),'alt');
    
    if(idx == -1)
        p = get(M.axmip,'CurrentPoint');
        p = round(p(1,1:2));
        [idx px py] = mip2ax(p(1),p(2),M.mipsize);
        
        px = min(M.mipsize(idx,2),max(1,round(px)));
        py = min(M.mipsize(idx,1),max(1,round(py)));
        pz = M.mipidx{idx}(py,px);
        
        %if right click, move to slice where max voxel was found
        if(rightclick)
            p = [px py pz];
        else
            p = [px py];
        end
    else
        p = get(M.ax(idx),'CurrentPoint');
        p = round(p(1,1:2));
    end
    
    if(idx == 1)
        M.cx = p(1);
        M.cz = p(2);
        if(numel(p) == 3)
            M.cy = p(3);
        end
    elseif(idx == 2)
        M.cy = p(1);
        M.cz = p(2);
        if(numel(p) == 3)
            M.cx = p(3);
        end        
    elseif(idx == 3)
        M.cx = p(1);
        M.cy = p(2);
        if(numel(p) == 3)
            M.cz = p(3);
        end        
    end
    
    update_linked_figures(M);
    setappdata(gcbf,'guidata',M);
    updatefig(M);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s1 s2 s3] = bgslice(Vbg, xyz, imgsz)
szratio = size(Vbg)./imgsz;
bxyz = round([xyz(:).'-1].*szratio + 1);
bxyz = min(size(Vbg),max(1,bxyz));

s1 = squeeze(Vbg(:,bxyz(2),:)).';
s2 = squeeze(Vbg(bxyz(1),:,:)).';
s3 = squeeze(Vbg(:,:,bxyz(3))).';

cmap = gray(256);
s1 = ind2rgb(ceil(size(cmap,1)*s1),cmap);
s2 = ind2rgb(ceil(size(cmap,1)*s2),cmap);
s3 = ind2rgb(ceil(size(cmap,1)*s3),cmap);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig_keypress(src,event)
M = getappdata(gcbf,'guidata');
movecursor = [];
if(strcmpi(event.Key,'h'))
    if(strcmpi(get(M.img(1),'visible'),'on'))
        set(M.img,'visible','off');
    else
        set(M.img,'visible','on');
    end
elseif(strcmpi(event.Key,'leftarrow'))
    movecursor = [-1 0];
elseif(strcmpi(event.Key,'rightarrow'))
    movecursor = [1 0];
elseif(strcmpi(event.Key,'uparrow'))
    movecursor = [0 1];
elseif(strcmpi(event.Key,'downarrow'))
    movecursor = [0 -1];
end

if(~isempty(movecursor) && ~isempty(M.curax))
    idx = M.curax;
    if(idx == 1)
        M.cx = M.cx + movecursor(1);
        M.cz = M.cz + movecursor(2);
    elseif(idx == 2)
        M.cy = M.cy + movecursor(1);
        M.cz = M.cz + movecursor(2);
    elseif(idx == 3)
        M.cx = M.cx + movecursor(1);
        M.cy = M.cy + movecursor(2);
    end
    update_linked_figures(M);
    setappdata(gcbf,'guidata',M);
    updatefig(M);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_linked_figures(M)
if(~isempty(M.figlink))
    allfigs = findobj(0,'type','figure');
    figtags = get(allfigs,'tag');
    linkedfigs = allfigs(regexpmatch(figtags,sprintf('^orthogui\\.%d$',M.figlink)));
    linkedfigs = setdiff(linkedfigs,M.hfig);
    if(~isempty(linkedfigs))
        for i = 1:numel(linkedfigs)
            orthogui(linkedfigs(i),'location',[M.cx M.cy M.cz]);
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hv hh] = plotcursor(ax,x,y,xl,yl,gap,varargin)
if(gap > 0)
    hv = plot(ax,[x x; x x]',[yl(1) y-gap; y+gap yl(2)]',varargin{:});
    hh = plot(ax,[xl(1) x-gap; x+gap xl(2)]',[y y; y y]',varargin{:});
else
    hv = plot(ax,[x x]',[yl(1) yl(2)]',varargin{:});
    hh = plot(ax,[xl(1) xl(2)]',[y y]',varargin{:});    
end
hv = reshape(hv,1,[]);
hh = reshape(hh,1,[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatecursor(hv,hh,x,y,xl,yl,gap,varargin)

if(numel(hv) > 1)
    set(hv(1),'XData',[x x],'YData',[yl(1) y-gap]);
    set(hv(2),'XData',[x x],'YData',[y+gap yl(2)]);

    set(hh(1),'XData',[xl(1) x-gap],'YData',[y y]);
    set(hh(2),'XData',[x+gap xl(2)],'YData',[y y]);
else
    set(hv,'XData',[x x],'YData',[yl(1) yl(2)]);
    set(hh,'XData',[xl(1) xl(2)],'YData',[y y]);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mipx mipy xl yl] = ax2mip(axidx,px,py,mipsize)

padsz = max(mipsize,[],1);
padgap = (repmat(padsz,size(mipsize,1),1)-mipsize)/2;
ox = (axidx-1)*padsz(2)+padgap(axidx,2);
oy = padgap(axidx,1);

mipx = px + ox;
mipy = py + oy;
xl = [ox ox+mipsize(axidx,2)];
yl = [oy oy+mipsize(axidx,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axidx px py] = mip2ax(mipx,mipy,mipsize)


padsz = max(mipsize,[],1);
padgap = (repmat(padsz,size(mipsize,1),1)-mipsize)/2;
axidx = floor(mipx/padsz(2)) + 1;
px = mipx - (axidx-1)*padsz(2) - padgap(axidx,2);
py = mipy -  padgap(axidx,1);
