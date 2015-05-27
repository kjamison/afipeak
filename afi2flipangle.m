function Astruct = afi2flipangle(dcmfile)
% Flip angle computation from:
% Yarnykh, V. L. (2007), Actual flip-angle imaging in the pulsed steady 
% state: A method for rapid three-dimensional mapping of the transmitted 
% radiofrequency field. Magn Reson Med, 57: 192–200. doi: 10.1002/mrm.21120

Astruct = [];

%dcmdir = '/home/range1-raid1/kjamison/shim_scratch/20140712-ST001-251833/MR-SE005-AFI';
[dcmdir, dcmfilename, ext] = fileparts(dcmfile);

if(~strcmpi(ext,'.dcm'))
    dcmdir = dcmfile;
end

dcmfiles = dir([dcmdir '/*.dcm']);

%TR1 = 25ms, TR2=105ms
D = [];
dcm = {};
for i = 1:numel(dcmfiles)
    dcmfile = [dcmdir '/' dcmfiles(i).name];
    dcm{i} = dicominfo(dcmfile);
    dcmimg = dicomread(dcmfile);
    if(isempty(D))
        D = zeros([size(dcmimg) numel(dcmfiles)]);
    end
    D(:,:,i) = dcmimg;
end
D1 = D(:,:,1:end/2);
D2 = D(:,:,end/2+1:end);



dcm_tr = dcm{1}.RepetitionTime;
ref_flipangle = dcm{1}.FlipAngle;
%dcm_csa = read_csa(dcm{1}.Private_0029_1010);
dcm_csa2 = read_csa(dcm{1}.Private_0029_1020);
numslices = parse_phoenix(dcm_csa2.MrPhoenixProtocol,'sKSpace.lImagesPerSlab');
Vref = parse_phoenix(dcm_csa2.MrPhoenixProtocol,'sTXSPEC.asNucleusInfo[0].flReferenceAmplitude');
try
	tr_offset = parse_phoenix(dcm_csa2.MrPhoenixProtocol,'sWiPMemBlock.alFree[11]');
catch err
	tr_offset = parse_phoenix(dcm_csa2.MrPhoenixProtocol,'sWipMemBlock.alFree[2]');
end

tr1 = dcm_tr-tr_offset/1000;
tr2 = dcm_tr+tr_offset/1000;

if(isempty(numslices) || isempty(Vref) || ~isnumeric(numslices) || ~isnumeric(Vref))
    return;
end

r = D2./D1;
n = tr2/tr1;

r(isnan(r)) = 0;
r(~isfinite(r)) = 0;


A = real(acos((n*r-1)./(n-r))*180/pi);
a_max = acos(-1/n)*180/pi;


A(~isfinite(A)) = 0;
A(A >= floor(a_max)) = 0;

Astruct = fillstruct(Vref,A,D1,D2,tr1,tr2,ref_flipangle);
