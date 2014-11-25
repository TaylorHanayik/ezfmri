function nii_batch (fmriname, t1name, TRsec, slice_order)
%preprocess fMRI data using standard settings
%  fmriname : name of 4D fMRI volumes
%  t1name   : name of anatomical scan
%  TRsec    : [optional] TR for fMRI data
%  slice_order : EPI slice order 0=auto,1=[1234],2=[4321],3=[1324],4=[4231], 5=[2413],6=[3142] 
% Motion correct, slice time corrects, coregisters, normalizes, smooths data
%  Final data has 'swa' prefix
%
%Examples 
%   nii_batch('fmri.nii','T1.nii');  %fMRI, T1
%   nii_batch('fmriblocks009.nii','T1s005.nii', 1.92, 1);  %fMRI, T1, 2 second TR, ascending
%   nii_batch(strvcat('fmriA.nii','fmriB.nii'),'T1.nii'); %fMRI (2 sessions), T1
% nii_batch(strvcat('/Users/rorden/Downloads/nii_process/2ses/afmri.nii','/Users/rorden/Downloads/nii_process/2ses/bfmri.nii'),'/Users/rorden/Downloads/nii_process/2ses/T1s005.nii', 1.92, 1);
resliceMM = 3; %resolution for reslicing data
if ~exist('fmriname','var'), fmriname = ''; end; %fMRI not specified
if ~exist('t1name','var'), t1name = ''; end; %T1 not specified
if ~exist('TRsec','var'), TRsec = 0; end; %TR not specified - autodetect
if ~exist('slice_order','var'), slice_order = 0; end; %TR not specified - autodetect
[fmriname, t1name, TRsec, slice_order, prefix] = validateInputsSub(fmriname, t1name, TRsec, slice_order); 
setOriginSub(strvcat(t1name, fmriname), 1); %#ok<REMFF1> %align image to anterior commissure
meanname = mocoSub(prefix, fmriname); %motion correct
prefix = slicetimeSub(prefix, fmriname, TRsec, slice_order); %slice-time correct
coregEstSub(t1name, meanname, prefix, fmriname); %make sure fMRI is aligned with T1
prefix = normSegSub(t1name, meanname, prefix, fmriname, resliceMM);
% % %ALTERNATE, with new segment...
% % % prefix = normNewSeg(t1name, meanname, prefix, fmriname);
prefix = smoothSub(8, prefix, fmriname); %smooth images
deleteImagesSub(prefix, fmriname); %delete intermediate images
stat_1st_block (prefix, fmriname,TRsec); %do first level statistics 1 SESSION
%ALTERNATE WITH TWO SESSIONS
%stat_1st_block_2ses (prefix, fmriname,TRsec); %do first level statistics 2 SESSIONS

%end nii_batch()

%---------- LOCAL FUNCTIONS FOLLOW

function deleteImagesSub(prefix, fmriname)
%delete intermediate images, e.g. if prefix is 'swa' then delete 'wa' and 'a' images
if length(prefix) < 2, return; end;
for s = 1 : length(fmriname(:,1))
    for i = 2 : length(prefix)
        [pth,nam,ext,vol] = spm_fileparts( deblank (fmriname(s,:))); %#ok<NASGU>
        nam = fullfile(pth,[prefix(i:length(prefix)), nam, ext]);
        if exist(nam, 'file')
            fprintf('Deleting  %s\n',nam );
            delete(nam);
        end;
    end;
end;
%end deleteImagesSub()

function prefix = smoothSub(FWHMmm, prefix, fmriname)
%blur images
fprintf('Smoothing with a %fmm FWHM Gaussian kernel\n',FWHMmm);
matlabbatch{1}.spm.spatial.smooth.data = getsesvolsSubFlat(prefix, fmriname);
matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHMmm FWHMmm FWHMmm];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);
prefix = ['s' prefix];
%smoothSub()

function prefix = normSegSub(t1name, meanname, prefix, fmriname, resliceMM)
%warp data to standard space with new segment
segnormSub(t1name);
segnormwriteSub(t1name, {t1name}, 1);
segnormwriteSub(t1name, {meanname}, resliceMM);
segnormwriteSub(t1name, getsesvolsSubFlat(prefix, fmriname), resliceMM);
prefix = ['w' prefix];
%end normSegSub()

function segnormSub(t1)
%estimate normalization based on unified segmentation normalization of T1
fprintf('Unified segmentation normalization of %s\n',t1);
matlabbatch{1}.spm.spatial.preproc.data = {t1};
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 0];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 0;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.opts.tpm = {fullfile(spm('Dir'),'tpm','grey.nii');fullfile(spm('Dir'),'tpm','white.nii');fullfile(spm('Dir'),'tpm','csf.nii')};
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2;2;2;4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
spm_jobman('run',matlabbatch);
%end segnormSub()

function segnormwriteSub(t1,mod, mm)
%reslice fMRI data based on previous segmentation-normalization
[pth,nam,ext,vol] = spm_fileparts( deblank (t1)); %#ok<NASGU,ASGLU>
fprintf('Applying unified segmentation normalization parameters from %s to %d image[s], resliced to %gmm\n',t1,length(mod(:,1)),mm);
matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {fullfile(pth,[ nam, '_seg_sn.mat'])};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = mod;
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [mm mm mm];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%end segnormwriteSub()

function coregEstSub(t1, meanfmri, prefix, fmriname)
%coregister fmri data to match T1 image
fprintf('Coregistering %s to match %s\n',meanfmri,t1);
fMRIses = getsesvolsSubHier(prefix, fmriname);
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {t1};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfmri};
matlabbatch{1}.spm.spatial.coreg.estimate.other = fMRIses;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregEstSub()

function [fmriname, t1name, TRsec, slice_order, prefix] = validateInputsSub(fmriname, t1name, TRsec, slice_order)
%check all inputs
if exist('spm','file')~=2; fprintf('%s requires SPM\n',which(mfilename)); return; end;
%X spm('Defaults','fMRI');
%X spm_jobman('initcfg'); % useful in SPM8 only
clear matlabbatch
prefix = ''; %originally '', if 'swa' it means we have smoothed, warped, aligned data
if ~exist('fmriname','var') || isempty(fmriname) %fMRI not specified
 fmriname = spm_select(2,'image','Select 4D fMRI volumes[s]');
end
fmriCell = getsesvolsSubHier(prefix, fmriname);
nSessions = numel(fmriCell);
nVol = sum(cellfun('prodofsize',fmriCell));
fprintf('fMRI has %d sessions for a total of %d volumes\n',nSessions, nVol);
if (nVol < 10) 
    error('Too few volumes: this script expects 4D fMRI images');
end
if (nSessions > 5) 
    error('Too many sessions: provide the FIRST volume from each 4D image');
end
if ~exist('t1name','var')  || isempty(t1name)
  t1name = spm_select(1,'image','Select T1 image volume');
end
%determine TR for fMRI data (in secounds)
if ~exist('TRsec','var')  || isempty(TRsec) || (TRsec == 0)
    TRsec = getTRSub(deblank (fmriname(1,:)));
end
if (TRsec ==0)
    answer = inputdlg('TR (sec)', 'Input required',1,{'2'});
    TRsec = str2double(answer{1});
end
%determine slice order
if ~exist('slice_order','var')  || isempty(slice_order) || (slice_order == 0)
    slice_order = getSliceOrderSub(deblank (fmriname(1,:)));
end
if (slice_order ==0)
    answer = inputdlg('slice order (1=[1234],2=[4321],3=[1324],4=[4231], 5=[2413],6=[3142])', 'Input required',1,{'1'});
    slice_order = str2double(answer{1});
end
%end validateInputsSub()

function slice_order =  getSliceOrderSub(fmriname)
%detect whether slices were acquired ascending, descending, interleaved
[pth,nam,ext,vol] = spm_fileparts( deblank(fmriname(1,:))); %#ok<NASGU>
fMRIname1 = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
fid = fopen(fMRIname1);
fseek(fid,122,'bof');
slice_order = fread(fid,1,'uint8');
fclose(fid);
%end getSliceOrderSub()

function tr =  getTRSub(fmriname)
%returns Repeat Time in seconds for volume fMRIname 
% n.b. for original images from dcm2nii - SPM will strip this information
hdr = spm_vol(fmriname);
if isfield(hdr(1,1).private.timing,'tspace')
  tr = hdr(1,1).private.timing.tspace;
else
  fprintf('%s error: unable to determine TR for image %s (perhaps SPM stripped this information)\n',mfilename,fmriname); 
  tr = 0;
end
%end getTRSub()

function coivox = setOriginSub(vols, modality)
%Use first image to estimate location of anterior commissure, apply to all
coivox = ones(4,1);
if ~exist('vols','var') %no files specified
 vols = spm_select(inf,'image','Reset origin for selected image(s) (estimated from 1st)');
end
vols = vol1OnlySub(vols); %only process first volume of 4D datasets...
if ~exist('modality','var') %no files specified
 modality = 1;
 fprintf('%s Modality not specified, assuming T1\n', mfilename);
end
%extract filename 
[pth,nam,ext, ~] = spm_fileparts(deblank(vols(1,:)));
fname = fullfile(pth,[nam ext]); %strip volume label
%report if filename does not exist...
if (exist(fname, 'file') ~= 2) 
 	fprintf('%s set origin error: unable to find image %s.\n',mfilename,fname);
	return;  
end;
hdr = spm_vol([fname,',1']); %load header 
img = spm_read_vols(hdr); %load image data
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
fprintf('%s center of brightness differs from current origin by %.0fx%.0fx%.0fmm in X Y Z dimensions\n',fname,XYZ_mm(1),XYZ_mm(2),XYZ_mm(3)); 
for v = 1:   size(vols,1) 
    fname = deblank(vols(v,:));
    if ~isempty(fname)
        [pth,nam,ext, ~] = spm_fileparts(fname);
        fname = fullfile(pth,[nam ext]); 
        hdr = spm_vol([fname ',1']); %load header of first volume 
        fname = fullfile(pth,[nam '.mat']);
        if exist(fname,'file')
            destname = fullfile(pth,[nam '_old.mat']);
            copyfile(fname,destname);
            fprintf('%s is renaming %s to %s\n',mfilename,fname,destname);
        end
        hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
        hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
        hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);
        spm_create_vol(hdr);
        if exist(fname,'file')
            delete(fname);
        end
    end
end%for each volume
coregSub(vols, modality);
for v = 1:   size(vols,1) 
    [pth, nam, ~, ~] = spm_fileparts(deblank(vols(v,:)));
    fname = fullfile(pth,[nam '.mat']);
    if exist(fname,'file')
        delete(fname);
    end
end%for each volume
%end setOriginSub()

function meanname = mocoSub(prefix, fmriname)
%motion correct fMRIdata
fMRIses = getsesvolsSubHier(prefix, fmriname);
fprintf('Motion correction\n');
matlabbatch{1}.spm.spatial.realign.estwrite.data = fMRIses;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; %0 <- do not reslice data!
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
[pth,nam,ext] = fileparts(deblank(fmriname(1,:)));
meanname = fullfile(pth,['mean', nam, ext]); %moco creates mean image, used for subsequent processing
%end mocoSub()

function prefix = slicetimeSub(prefix, fmriname, TRsec, slice_order) 
%slice time correct data
kNIFTI_SLICE_SEQ_INC = 1; %1,2,3,4
kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
%kNIFTI_SLICE_ALT_INC = 3; %1,3,2,4 Siemens: interleaved with odd number of slices, interleaved for other vendors
%kNIFTI_SLICE_ALT_DEC = 4; %4,2,3,1 descending interleaved
kNIFTI_SLICE_ALT_INC2 = 5; %2,4,1,3 Siemens interleaved with even number of slices 
kNIFTI_SLICE_ALT_DEC2 = 6; %3,1,4,2 Siemens interleaved descending with even number of 
[pth,nam,ext,vol] = spm_fileparts( deblank(fmriname(1,:))); %#ok<NASGU>
fMRIname1 = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
if slice_order == 0 %attempt to autodetect slice order
    fid = fopen(fMRIname1);
    fseek(fid,122,'bof');
    slice_order = fread(fid,1,'uint8');
    fclose(fid);
    if (slice_order > kNIFTI_SLICE_UNKNOWN) && (slice_order <= kNIFTI_SLICE_ALT_DEC2)
        fprintf('Auto-detected slice order as %d\n',slice_order);
    else
        error('Error: unable to auto-detect slice order. Please manually specify slice order or use recent versions of dcm2nii.\n');
    end;
end
hdr = spm_vol([fMRIname1 ',1']);
if TRsec == 0
    TRsec = hdr.private.timing.tspace;
    if TRsec == 0
        error('%s error: unable to auto-detect slice timing. Please manually specify slice order or use recent versions of dcm2nii.\n');
    end; 
end
nslices = hdr.dim(3);
if nslices <= 1 %automatically detect TR
    error('Fatal Error: image %s does not have multiple slices per volume - slice time correction inappropriate. Please edit m-file named %s\n',fMRIname,which(mfilename));
end;
if (slice_order == kNIFTI_SLICE_ALT_INC2) || (slice_order == kNIFTI_SLICE_ALT_DEC2) %sequential
    isSiemens = true;
end;
if (slice_order == kNIFTI_SLICE_SEQ_INC) || (slice_order == kNIFTI_SLICE_SEQ_DEC) %sequential
	so = 1:1:nslices;
else % if sequential else Interleaved
	if (mod(nslices,2) == 0) && (isSiemens) %even number of slices, Siemens
		so =[2:2:nslices 1:2:nslices ];
	else
		so =[1:2:nslices 2:2:nslices];
	end
end
if (mod(slice_order,2) == 0) %isDescending
	so = (nslices+1)-so;
end; %isDescending
TA = (TRsec/nslices)*(nslices-1);
fprintf('  Slice order=%d, slices=%d, TR= %0.3fsec, TA= %fsec, referenced to 1st slice.\n', slice_order,nslices, TRsec,TA);
if (TRsec < 0.1) || (TRsec > 5.0) 
    fprintf('  Aborting: strange Repeat Time (TR). Please edit the m-file.\n');
    if  (TRsec > 5.0) 
          fprintf('  Long TR often used with sparse imaging: if this is a sparse design please set the TA manually.\n');
    end; 
    if  (TRsec < 0.1) 
          fprintf('  Short TR may be due to DICOM-to-NIfTI conversion. Perhaps use dcm2nii.\n');
    end; 
    return;
end; %unusual TR
fMRIses = getsesvolsSubHier(prefix, fmriname);
matlabbatch{1}.spm.temporal.st.scans = fMRIses;                               
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TRsec;
matlabbatch{1}.spm.temporal.st.ta = TA;
matlabbatch{1}.spm.temporal.st.so = so;
matlabbatch{1}.spm.temporal.st.refslice = so(1); %set slice order to the first slice http://www.alivelearn.net/?p=1037
%so(1) is first acquired slice, set it as the refernece
%  to see how this works add the line "fprintf('spm_slice_timing slice %d shift %f\n',k, shift amount);' if the 'for k =' loop of spm_slice_timing.m
fprintf('Setting reference slice as %d\n',so(1));
matlabbatch{1}.spm.temporal.st.prefix = 'a';
spm_jobman('run',matlabbatch);
prefix = ['a' prefix]; %new files will have 'a' appended to filename
%end slicetimeSub()

function coregSub(vols, modality)
if modality == 2
    template = fullfile(spm('Dir'),'templates','T2.nii');
elseif modality == 3
    template  = fullfile(spm('Dir'),'toolbox','Clinical','scct.nii');
else
    template = fullfile(spm('Dir'),'templates','T1.nii');
end
if ~exist(template,'file')
    error('Unable to find template named %s\n', template);
end
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols(1,:)),',1']};%{'/Users/rorden/Desktop/3D.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(vols(2:end,:));% {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregSub()

function vols = vol1OnlySub(vols)
%only select first volume of multivolume images '/dir/img.nii' -> '/dir/img.nii,1', '/dir/img.nii,33' -> '/dir/img.nii,1'
oldvols = vols;
vols = [];
for v = 1:   size(oldvols,1) 
    [pth,nam,ext, ~] = spm_fileparts(deblank(oldvols(v,:)));
    vols = strvcat(vols, fullfile(pth, [ nam ext ',1']) ); %#ok<REMFF1>
end
%end vol1OnlySub()

function [fMRIses] = getsesvolsSubHier(prefix, fmriname)
%load all images from all sessions... AS MULTIPLE SESSIONS (e.g. moco)
nsessions = length(fmriname(:,1));
fMRIses = cell(nsessions,1);
for s = 1 : nsessions 
	[pth,nam,ext,vol] = spm_fileparts( deblank (fmriname(s,:))); %#ok<NASGU>
	sesname = fullfile(pth,[prefix, nam, ext]);
    fMRIses(s,1) = {getsesvolsSub(sesname)};
end;
%end getsesvolsSubHier()

function [fMRIses] = getsesvolsSubFlat(prefix, fmriname)
%load all images from all sessions... AS SINGLE SESSION (e.g. norm writing)
nsessions = length(fmriname(:,1));
fMRIses = '';
for s = 1 : nsessions 
	[pth,nam,ext,vol] = spm_fileparts( deblank (fmriname(s,:))); %#ok<NASGU>
	sesname = fullfile(pth,[prefix, nam, ext]);
    fMRIses = [fMRIses; getsesvolsSub(sesname)]; %#ok<AGROW>
end;
% end getsesvolsSubFlat()

function [sesvols] = getsesvolsSub(sesvol1)
% input: single volume from 4D volume, output: volume list
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1';'img.nii,2';'img.nii,3'}
[pth,nam,ext,vol] = spm_fileparts( deblank (sesvol1)); %#ok<NASGU>
sesname = fullfile(pth,[nam, ext]);
hdr = spm_vol(sesname);
nvol = length(hdr);
sesvols = cellstr([sesname,',1']);
if nvol < 2, return; end;
for vol = 2 : nvol
    sesvols = [sesvols; [sesname,',',int2str(vol)]  ]; %#ok<AGROW>
end;
%end getsesvolsSub()
