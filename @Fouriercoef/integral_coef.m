function integral_coef(obj,kk)
           freq=get_freq(obj,kk);
          n=numel(obj.coef{kk}.cos);
          
         
           
           for i=1:1:n
               integraltotal(i)=(obj.coef{kk}.cos(i)/2+j*obj.coef{kk}.sin(i)/(-2))/(j*2*pi*freq*i);
           end
                     
           obj.coef{kk}.integralcos = real(integraltotal) *2;
           obj.coef{kk}.integralsin = imag(integraltotal) * -2;
           obj.coef{kk}.integraldc  = obj.coef{kk}.dc * pi;
     
      end %integral func