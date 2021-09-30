function CGHdata = func_roundCGH(W,iter,sizex,sizey)

  
N=2048;
i = sqrt(-1);
A = zeros(N,N);

for u = 1:N
    for v = 1:N
        
       if (u-N/2)^2+(v-N/2)^2<W^2
           A(u,v)=1;
       end
    end
end

PH=rand(N,N);
I=A.*PH;
n=iter;
%Error=zeros(n);

for p=1:n
    B=fftshift(fft2(fftshift(I)));
    B1=angle(B);
    C=A.*exp(i*B1); 
    
    D=ifftshift(ifft2(fftshift(C)));
    D1=angle(D);
    I=exp(i*D1);  
end

CGHdata = I(N/2-sizey/2:N/2+sizey/2-1,N/2-sizex/2:N/2+sizex/2-1);

end