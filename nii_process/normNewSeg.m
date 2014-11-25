function prefix = normNewSeg(t1name, meanname, prefix, fmriname)
%warp data to standard space with new segment
% t1name : filename of T1-weighted image
% meanname : name of mean fMRI data
% prefix   : prefix appended to name of fMRI data
% fmriname : base name(s) of 4D fMRI images
%Example
% normNewSeg('T1.nii','meanfmri.nii','fmri.nii','');
% normNewSeg('T1.nii','meanfmri.nii','fmri.nii','a'); %normalize 'afmri.nii'
newSeg(t1name); %normalize images
newSegWriteSub(t1name, t1name, '');
newSegWriteSub(t1name, meanname, '');
prefix = newSegWriteSub(t1name, fmriname, prefix);
%end normNewSeg()

function  prefix = newSegWriteSub(t1name, warpname, prefix)
%reslice img using pre-existing new-segmentation deformation field
[pth,nam,ext, vol] = spm_fileparts(t1name); %#ok<NASGU>
defname = fullfile(pth,['y_' nam ext]);
if ~exist(defname,'file')
    error('Unable to find new-segment deformation image %s',defname);
end
matlabbatch{1}.spm.util.defs.comp{1}.def = {defname};
matlabbatch{1}.spm.util.defs.ofname = '';
warpses = getsesvolsSubFlat(warpname, prefix);
matlabbatch{1}.spm.util.defs.fnames = warpses;
matlabbatch{1}.spm.util.defs.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.interp = 1;
spm_jobman('run',matlabbatch);
prefix = ['w' prefix];
%end for subfunction newsegwritesub

function newSeg(t1name)
%apply new segment - return name of warping matrix
template = fullfile(spm('Dir'),'toolbox','Seg','TPM.nii');
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
matlabbatch{1}.spm.tools.preproc8.channel.vols = {t1name};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[template ',1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[template ',2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[template ',3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[template ',4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[template ',5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[template ',6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];
spm_jobman('run',matlabbatch);
%end newSeg()

function [fMRIses] = getsesvolsSubFlat(fmriname, prefix)
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
