function timeshift_coef(obj,t0,kk)
          freq=get_freq(obj,kk);
          n=numel(obj.coef{kk}.cos);
          
           
           
           for i=1:1:n
               timeshifttotal(i)=(obj.coef{kk}.cos(i)/2+j*obj.coef{kk}.sin(i)/(-2)) * exp(j*2*pi*freq*i*t0);
           end
                     
           obj.coef{kk}.timeshiftcos = real(timeshifttotal) *2;
           obj.coef{kk}.timeshiftsin = imag(timeshifttotal) * -2;
     
      end %timeshift func