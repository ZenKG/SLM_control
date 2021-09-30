function SLMdata = func_quantSLM10bit(phaseMap,Gmax,sizex,sizey)

SLMdata = uint8(phaseMap.*Gmax);

end
