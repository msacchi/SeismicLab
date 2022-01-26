function [out]=Operator_Radon_Stolt(in,PARAM,transform)

%PARAM.time;        %time vector  
%PARAM.receivers;   %x-axis vector 
%PARAM.v;           %velocities vector 
%PARAM.shift         %Shifted apexes location vector 
%PARAM.tpad;        %zero padding ratio for time 
%PARAM.xpad         %zero padding ratio for x-axis 
%PARAM.freq_cut     %cutt off frequency 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Input Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nt=length(PARAM.time);                  %time axis     
dt=(PARAM.time(nt)-PARAM.time(1))/(nt-1); 
nx=length(PARAM.shift_x);             %horizontal axis (apexes)
nR=length(PARAM.receivers);         %horizontal axis (receivers)
dx=(PARAM.shift_x(nx)-PARAM.shift_x(1))/(nx-1); 
nv=length(PARAM.v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Cross-Correlated with wavelet (if adjoint) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(transform,'a') %Adj operator m(tau,v,x) =L' d(t,h) 
    if isequal(PARAM.use_w,'yes')
        in= conv_xcorr(in,PARAM.w,transform); 
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Extend Apexes Range 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(transform,'a') %Adjoint 
    [in]=truncating_padding(in,PARAM,transform);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Zero Padding                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmax=PARAM.time(nt)-PARAM.time(1);  
tpad=PARAM.tpad*tmax;     nt_pad = round((tmax+tpad)/dt+1);
xmax=PARAM.shift_x(nx)-PARAM.shift_x(1);  
xpad=PARAM.tpad*xmax;     nx_pad = round((xmax+xpad)/dx+1);  
if isequal(transform,'a')
    in = [in; zeros(nt_pad-nt,nx)];
    in = [in  zeros(nt_pad,nx_pad-nx)];
elseif isequal(transform,'f') 
    tmp=zeros(nt_pad,nv,nx_pad);
    tmp(1:nt,:,1:nx)=in;
    in=tmp; 
end 
nf=2^(nextpow2(nt_pad));    
nkx=2^(nextpow2(nx_pad));   
%%Compute the frequency and wavenumber axes 
nf2 = nf/2+2;                           nkx2=nkx/2+2;
ifreq=1:1:nf;                           ikx=1:1:nkx; 
ifreq2=ifreq-1-nf*floor(ifreq/nf2);     ikx2=ikx-1-nkx*floor(ikx/nkx2);
freq=ifreq2/nf/dt;                      kx=ikx2/nkx/dx;   
%cut off freq 
dfreq=1/nf/dt;
ifreq_cut=round(PARAM.freq_cut/dfreq)+1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Forward FTT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(transform,'a')            %Adj operator m(tau,v,x) =L' d(t,h)
    %kx_shift=exp(-2*pi*1i*kx(ikx)*PARAM.shift(1));
    in_fk=fft2(in,nf,nkx); %2D FFT+1 d(t,h) ----> d(w,kx)
    out_fk=zeros(nf/2+1,nv,nkx);    %Memory for m(w,v,kx) 
    out_fx=zeros(nf,nv,nkx);        %Memory for m(w,v,x)
    out=zeros(nt,nv,nx);            %Memory for m(tau,v,x) 
elseif isequal(transform,'f')         %Forward d(t,h) = L m(tau,v,x)
    %kx_shift=exp(2*pi*1i*kx(ikx)*PARAM.shift(1));
    in_fk=zeros(nf,nv,nkx);         %Memory for m(w,v,kz)
    for iv=1:nv;                    %2D FFT+1 m(t,v,x) ----> m(w,v,kx)
        in_fk(:,iv,:)=fft2(squeeze(in(:,iv,:)),nf,nkx); 
    end 
    out_fk=zeros(nf/2+1,nkx);       %Memory for d(w,kx)
    out_fx=zeros(nf,nkx);           %Memory for d(w,x)
    out=zeros(nt,nx);               %Memory for d(t,h)
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Stolt mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for iv=1:nv  
    for ikx=1:nkx 
        for ifreq_tau=1:ifreq_cut  
            freq_map2=kx(ikx)^2*PARAM.v(iv)^2+freq(ifreq_tau)^2; 
            freq_map=sqrt(freq_map2); 
            ifreq_map=ceil(freq_map/dfreq)+1;     % mapped freq index  
            if(freq_map2 ~=0)
                scale=freq(ifreq_tau)/freq_map; 
            else 
                scale=1;
            end 
            if(0<ifreq_map && ifreq_map<ifreq_cut);
                if isequal(transform,'a') %adj m =L' d
                    out_fk(ifreq_tau,iv,ikx)=out_fk(ifreq_tau,iv,ikx)+in_fk(ifreq_map,ikx)*scale;
                elseif isequal(transform,'f') %forward d =L m
                    out_fk(ifreq_map,ikx)=out_fk(ifreq_map,ikx)+in_fk(ifreq_tau,iv,ikx)*scale;  
                end
            end 
        end 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Inverse FTT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compute the -ve frequencies using Symmetry 
if isequal(transform,'a')
    for iv=1:nv
        out_fx(1:nf/2+1,iv,:)=ifft(squeeze(out_fk(:,iv,:)),[],2);
        out_fx(nf/2+2:nf,iv,:)=conj(flipud(squeeze(out_fx(2:nf/2,iv,:))));
    end 
elseif isequal(transform,'f')
    out_fx(1:nf/2+1,:)=ifft(out_fk,[],2);
    out_fx(nf/2+2:nf,:)=conj((flipud(out_fx(2:nf/2,:)))); 
end
%%Inverse Fourier Transform and Remove Zero Pads 
if isequal(transform,'a') %adj m =L' d
    for iv=1:nv
        temp=(ifft(squeeze(out_fx(:,iv,:)),[],1));%2D FFT-1  m(tau,v,x)
        out(:,iv,:)=real(temp(1:nt,1:nx));
    end 
elseif isequal(transform,'f') %forward d=L m
    temp=(ifft(out_fx,[],1));%2D FFT-1  d(t,x)
    out=real(temp(1:nt,1:nx));   
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Remove extended offsets (Forward)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(transform,'f') 
    [out]=truncating_padding(out,PARAM,transform);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Convolve with wavelet (Forward)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(transform,'f') 
    if isequal(PARAM.use_w,'yes')
        out=conv_xcorr(out,PARAM.w,transform);  %Convolve with wavelet
    end 
end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Functions used in the Stolt2 Operator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Function to add and remove zero padding to extend apexes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]=truncating_padding(input,PARAM,transform)
% Forward -->truncate 
if isequal(transform,'f')  
    [nt,~]=size(input); 
    output=zeros(nt,length(PARAM.receivers));
    for iR=1:length(PARAM.receivers)
        [index]=Find_element_index(PARAM.shift_x, PARAM.receivers(iR));
        if index > 0
            output(:,iR)=input(:,index); 
        end 
    end 
end 
%Adjoint -->padd 
if isequal(transform,'a') %adjoint
    [nt,~]=size(input); 
    output=zeros(nt,length(PARAM.shift_x));
    for iSh=1:length(PARAM.shift_x)
        [index]=Find_element_index(PARAM.receivers,PARAM.shift_x(iSh));
        if index > 0
            output(:,iSh)=input(:,index); 
        end 
    end 
end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Function to find element location for zero padding  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [index]=Find_element_index(input, element)
for i=1:length(input)
   if input(i)== element
       index=i; 
       break; 
   else 
       index=0; 
   end 
end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Function to cross-correlate/convolve with wavelet 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out]=conv_xcorr(x,w,transform)
[nt,nx]=size(x); 
out=zeros(nt,nx);
if isequal(transform,'f')
    for ix=1:nx
        tr=conv(x(:,ix),w);
        out(1:nt,ix)=tr(1:nt);
    end 
end 
if isequal(transform,'a')
    for ix=1:nx
        tr=xcorr(x(:,ix),w);
        out(1:nt,ix)=tr(nt:end);
    end 
end 
end

