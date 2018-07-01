function record()
check = exist ('ambience.wav');
fs = 200000;
b = 24;
coder.extrinsic('audiorecorder');
recObj = audiorecorder(fs,b,1);
disp('Start recording.')
coder.extrinsic('recordblocking');
recordblocking(recObj, 3);
disp('End of recording.');
coder.extrinsic('getaudiodata');
data = getaudiodata(recObj);
filename = 'ambience.wav';
if check == 0
    coder.extrinsic('audiowrite');
    audiowrite(filename,data,fs);
else
    oldFilename = 'ambience.wav';
    check = 2;
    i = 1;
    while check >1,
    newFilename = ['ambience' num2str(i) '.wav'];
    check = exist(newFilename);
    i = i+1;
    end
    movefile(oldFilename,newFilename);
    coder.extrinsic('audiowrite');
    audiowrite(filename,data,fs);
end    

end