%Convert string to binary
function n=text2bin(textstring)
bintext = dec2bin(double(textstring),7);%text into binary..
[rp,cp] = size(bintext);
n = str2num(reshape(bintext',1,rp*cp)')';