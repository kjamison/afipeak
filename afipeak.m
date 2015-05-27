function afipeak_kj(varargin)

% afipeak_kj
%
% make sure this is in the Matlab path.
% (ie: open in matlab editor, run it and click "Add to path" when it
%   prompts you
afi_target = 52;

roisize = 3;
%% get dcm

%dcmstart = [getuserdir '/MR-SE005-AFI']; %test
%dcmstart = '/export/dicom'; %shim7tas dicom
%dcmstart = '/export/raid1/scratch'; %shim7tas scratch
dcmstart = getuserdir;
xyz_start = [];

if(nargin >= 1 && exist(varargin{1},'file'))
    dcmstart = varargin{1};
end

if(nargin >= 2)
    xyz_start = varargin{2};
end

dcmstart_dir = '';
if(exist(dcmstart,'dir'))
    dcmstart_dir = dcmstart;
elseif(exist(dcmstart,'file'))
    [dcmstart_dir,~,~] = fileparts(dcmstart);
end

dcmfile = '';
if(~isempty(dcmstart_dir))
    D = dir(dcmstart_dir);
    D = D(~arrayfun(@(d)(all(d.name=='.')),D));
    if(all(arrayfun(@(d)(~d.isdir && regexpimatch(d.name,'\.dcm$')),D)))
        dcmfile = fullfile(dcmstart_dir,D(1).name);
    end
end

if(isempty(dcmfile))
    [dcmfilename,dcmdir] = uigetfile({'*.dcm';'*.*'},'Select AFI DICOM file',dcmstart);
    if(~isnumeric(dcmfilename))
        dcmfile = fullfile(dcmdir,dcmfilename);
    end
end

if(isempty(dcmfile) || ~exist(dcmfile,'file'))
    if(~usejava('desktop'))
        exit;
        
    end
    return;
end
%

%%
fprintf('\n####################################################\n');
fprintf('Reading AFI: %s\n\n',dcmfile);

%%

Astruct = afi2flipangle(dcmfile);
A = Astruct.A;
ref_voltage = Astruct.Vref;
ref_flipangle = Astruct.ref_flipangle;

BG = (Astruct.D1+Astruct.D2)/2;

%% 
sz3d = size(A);

% A reasonable starting guess for whole-brain AFI
xyz_expected = sz3d.*[.66 .5 .5];

%% remove SOME noisy voxels
Ahpf = A-smooth3(A,'gaussian',[5 5 5]);
dmask = roimask(sz3d,xyz_expected,20);
Aclean = A.*(abs(Ahpf - median(Ahpf(dmask(:)))) < 1*std(Ahpf(dmask(:))));
%Aclean(Aclean == 0) = nan;

%neighbor_cutoff = 4/26;
%neighbor_cutoff = .75;
%Aclean(smooth3(Aclean>0,'gaussian',[3 3 3]) < neighbor_cutoff) = 0;
%Aclean(isnan(smooth3(Aclean,'gaussian',[3 3 3]))) = nan;

%Aclean = smooth3(Aclean,'gaussian',[3 3 3]);


Aclean(Aclean == 0) = nan;

% find center of mass on smoothed volume
Asmooth = smooth3(Aclean,'box',[5 5 5]);
[my, mx, mz] = meshgrid(1:sz3d(2),1:sz3d(1),1:sz3d(3));
cx = nansum(mx(:).*Asmooth(:))/nansum(Asmooth(:));
cy = nansum(my(:).*Asmooth(:))/nansum(Asmooth(:));
cz = nansum(mz(:).*Asmooth(:))/nansum(Asmooth(:));
xyz_cmass = [cx cy cz];

%% find peak 

%if(isempty(xyz_start))
%    xyz_start = xyz_expected;
%end

if(isempty(xyz_start))
    xyz_peak0 = xyz_cmass;
else
    xyz_peak0 = roipeak(Aclean,xyz_start,20);
end

xyz_peak = findpeak(Aclean,xyz_peak0);

dmask = roimask(sz3d,xyz_expected,20);
Arange = [0 nanmax(Aclean(dmask(:)))];

%% initialize figures
orthfig = orthogui(Aclean,'bg',BG,'loc',round(xyz_peak),'clim',Arange);
%set(orthfig,'toolbar','none','menubar','none');

figinfo = hgload('info_fig.fig');
bgcol = get(figinfo,'color');
p_info = getpixelposition(figinfo);

bgcol = [1 1 1];
controls = uistruct(figinfo);
controlnames = fieldnames(controls);
for i = 1:numel(controlnames)
    ctrl = controls.(controlnames{i});
    switch get(ctrl,'Type')
        case 'uicontrol'
            set(ctrl,'backgroundcolor',bgcol);
        case 'figure'
            set(ctrl,'color',bgcol);
        otherwise
    end
end

gdata = guidata(orthfig);
gdata.dcmfile = dcmfile;
gdata.position3d = xyz_peak;
gdata.sz3d = sz3d;
gdata.Aclean = Aclean;
gdata.ref_voltage = ref_voltage;
gdata.afi_target = afi_target;
gdata.ref_flipangle = ref_flipangle;
gdata.figinfo = figinfo;
gdata.roisize = roisize;
gdata.controls = controls;
guidata(orthfig,gdata);



printfcn = @(xyz,ref_voltage,ref_flipangle,afi_target,roisize)(print_shim_info(xyz,A,ref_voltage,ref_flipangle,afi_target,roisize,controls));
ortho_cb = @(xyz)(datacursor_updatefcn(xyz,orthfig,printfcn));
orthogui(orthfig,'callback',ortho_cb);
p_orth = getpixelposition(orthfig);
p_info_new = [p_orth(1)-p_info(3) p_orth(2) p_info(3) p_info(4)];
setpixelposition(figinfo,p_info_new);

set(controls.butFindPeak,'callback',@(~,~)(findpeak_callback(orthfig,ortho_cb)));
set(controls.butNewAFI,'callback',@(~,~)(load_new_afi(orthfig,figinfo)));
set(controls.butSetFlip,'callback',@(~,~)(set_afi_target(orthfig,ortho_cb)));
set(controls.butSetROI,'callback',@(~,~)(set_roisize(orthfig,ortho_cb)));
set(figinfo,'CloseRequestFcn',@(~,~)(closewindow([orthfig figinfo])));
set(orthfig,'CloseRequestFcn',@(~,~)(closewindow([orthfig figinfo])));

ortho_cb(round(xyz_peak));

%%
function datacursor_updatefcn(xyz,orthfig,printfcn)

gdata = guidata(orthfig);
gdata.position3d = xyz;
if(gdata.roisize > 1)
    dmask = roimask(gdata.sz3d,round(xyz),gdata.roisize);
    orthogui(orthfig,'setalphadata',dmask*.5+.5);
else
    orthogui(orthfig,'setalphadata',[]);
end
guidata(orthfig,gdata);

printfcn(xyz,gdata.ref_voltage,gdata.ref_flipangle,gdata.afi_target,gdata.roisize);

%%
function print_shim_info(xyz,A,ref_voltage,ref_flipangle,afi_target,roisize,controls)

xyz_round = round(xyz);

peakheight = A(xyz_round(1),xyz_round(2),xyz_round(3));
new_ref = ref_voltage*afi_target/peakheight;

if(roisize > 1)
    dmask = roimask(size(A),xyz_round,roisize);
    peakheight_roi = nanmean(A(dmask(:)));
    new_ref_roi = ref_voltage*afi_target/peakheight_roi;
else
    peakheight_roi = [];
    new_ref_roi = [];
end

flip_pct = 100*(afi_target-ref_flipangle)/ref_flipangle;

set(controls.txtX,'String',sprintf('%.0f',xyz_round(1)));
set(controls.txtY,'String',sprintf('%.0f',xyz_round(2)));
set(controls.txtZ,'String',sprintf('%.0f',xyz_round(3)));
set(controls.txtFlip,'String',sprintf('%.3f',peakheight));
set(controls.txtNewRef,'String',sprintf('%.0f',new_ref));
set(controls.txtFlip_roi,'String',sprintf('%.3f',peakheight_roi));
set(controls.txtNewRef_roi,'String',sprintf('%.0f',new_ref_roi));
set(controls.txtInitRef,'String',sprintf('%.0f',ref_voltage));
set(controls.txtFlipTarget,'String',sprintf('%.0f (%.0f+%.0f %%)', afi_target, ref_flipangle, flip_pct));
set(controls.txtROI,'String',sprintf('%.0f',roisize));

%%
function findpeak_callback(orthfig,ortho_cb)
gdata = guidata(orthfig);


xyz_peak = findpeak(gdata.Aclean,gdata.position3d);
ortho_cb(round(xyz_peak));
orthogui(orthfig,'loc',round(xyz_peak));

%%
function set_afi_target(orthfig, ortho_cb)
gdata = guidata(orthfig);

tmp_str = inputdlg('Enter target flip angle:','Set Target',1,{num2str(gdata.afi_target)});
if(isempty(tmp_str) || isempty(str2num(tmp_str{1})))
    return;
end

gdata.afi_target = str2num(tmp_str{1});
guidata(orthfig,gdata);

xyz_peak = gdata.position3d;
ortho_cb(round(xyz_peak));
orthogui(orthfig,'loc',round(xyz_peak));


%%
function set_roisize(orthfig, ortho_cb)
gdata = guidata(orthfig);

tmp_str = inputdlg('Enter ROI radius:','ROI Size',1,{num2str(gdata.roisize)});
if(isempty(tmp_str) || isempty(str2num(tmp_str{1})))
    return;
end
roisize = str2num(tmp_str{1});
gdata.roisize = max(1,roisize);
guidata(orthfig,gdata);

xyz_peak = gdata.position3d;
ortho_cb(round(xyz_peak));
orthogui(orthfig,'loc',round(xyz_peak));

%%
function closewindow(figs)
if(usejava('desktop'))
    delete(figs(ishandle(figs)));
else
    exit;
end

%%
function load_new_afi(orthfig,figinfo)
gdata = guidata(orthfig);
dcmstart = gdata.dcmfile;
figs = [orthfig figinfo];
eval('delete(figs(ishandle(figs))); afipeak(dcmstart)');
