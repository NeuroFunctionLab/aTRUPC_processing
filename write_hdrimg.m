function write_hdrimg(imgarray, imgfname, voxsize, dt, ss)
% write 3D/4D images into one analyze format (hdr/img)
%
% 	imgarray: 3D/4D mat
% 	imgfname: full directory of the file intended to write
% 	voxsize: 1*3 mat storing the voxel size in mm
% 	dt: datatype
% 		{'uint8','int16','int32','float32','float64','int8','uint16','uint32'};
% 		[     2       4       8        16        64    256      512      768];
% 	ss: 1 - no rescale

imgsize = [size(imgarray) 1 1];
matsize = imgsize(1:3);
numvol  = imgsize(4);

outVol.fname    = imgfname;
outVol.dim      = matsize;
outVol.dt       = [dt,0];
outVol.pinfo    = [ss;0;0];
outVol.descrip  = 'Analyze file';
outVol.private  = '';
outVol.mat      = [ -voxsize(1) 0 0 (floor(matsize(1)/2)+0.5)*voxsize(1);
                    0 voxsize(2) 0 -(floor(matsize(2)/2)+0.5)*voxsize(2);
                    0 0 voxsize(3) -(floor(matsize(3)/2)+0.5)*voxsize(3);
                    0 0 0 1];
                
for ii = 1:numvol
    outVol.n = [ii,1];
    spm_write_vol(outVol,imgarray(:,:,:,ii));
end

