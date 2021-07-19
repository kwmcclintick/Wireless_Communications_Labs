% Function convert binary to ASCII string...
function textstring = bin2text(n)
rp = floor(length(n)/7);
rez = num2str(n(1:7*rp)')';
textstring = char(bin2dec(reshape(rez,7,rp)'))';
