function SLMdata = func_quantSLM10bit(phaseMap,Gmax,sizex,sizey)

phase = uint16(phaseMap.*Gmax);

% transform to 10bit bmp
red1 = uint16(2^9+2^8+2^7);
red2 = bitand(phase,red1);
R = uint8(bitshift(red2,-2));

green1 = uint16(2^6+2^5+2^4);
green2 = bitand(phase,green1);
G = uint8(bitshift(green2,1));

blue1 = uint16(2^3+2^2+2^1+2^0);
blue2 = bitand(phase,blue1);
B = uint8(bitshift(blue2,4));

SLMdata=uint8(zeros(sizey,sizex,3));
SLMdata(:,:,1) = R;
SLMdata(:,:,2) = G;
SLMdata(:,:,3) = B;
end
