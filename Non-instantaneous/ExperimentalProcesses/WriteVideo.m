function [] = WriteVideo(outName, vid)

warning('off','MATLAB:audiovideo:VideoWriter:mp4FramePadded');
v = VideoWriter(outName,'MPEG-4');
open(v)
writeVideo(v, vid);
close(v)
implay(strcat(outName,'.mp4'));

end