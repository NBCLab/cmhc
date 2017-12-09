addpath /home/applications/spm12

maskimg = spm_vol('Tal_2mm_mask.nii.gz');
mask = spm_read_vols(maskimg);
thelocs = find(mask);

for a = 1:numel(ExpImg)
    temp = ExpImg(a).ModActs;
    voxval(:,a) = temp(thelocs);
    clear temp
end

corrmat = corr(voxval);
thelink = linkage(corrmat, 'average', 'correlation');
dendrogram(thelink, 0)
