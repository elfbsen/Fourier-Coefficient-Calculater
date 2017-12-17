function derivative_coef(obj,kk)
          freq=get_freq(obj,kk);
          n=numel(obj.coef{kk}.cos);
          
          
           
           for i=1:1:n
               derivativetotal(i)=(obj.coef{kk}.cos(i)/2+j*obj.coef{kk}.sin(i)/(-2))*(j*2*pi*freq*i);
           end
                     
           obj.coef{kk}.derivativecos = real(derivativetotal) *2;
           obj.coef{kk}.derivativesin = imag(derivativetotal) * -2;
          
      end %derivative func