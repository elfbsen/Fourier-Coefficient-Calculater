classdef Substairs < Absinterpol
    
   
   methods 
       function obj = Substairs(d,t)
            if nargin > 0                         
                if isnumeric(d)
                     if isnumeric(t)
                     
                        obj.ppval=interp1(t,d,'nearest','pp');    
                         
                     end
                end 

            end
       end %constructor        
       function intervec=interpolation_vector(obj,pointnum);
          
           intervec=ppval(obj.ppval,...
               linspace(obj.ppval.breaks(1),...
               obj.ppval.breaks(obj.ppval.pieces +1),pointnum));
           
       end %func interpolation_vector
       function result=integral_for_dc(obj)
            period=(obj.ppval.breaks(obj.ppval.pieces +1)-obj.ppval.breaks(1));
            result=2/period * [ obj.ppval.breaks(:).^1 ];
            result1=result(2:(obj.ppval.pieces+1),:);
            result=result1-result(1:obj.ppval.pieces,:);
        end %integral_for_dc func
       function result=integral_for_harmonics(obj,howmanyn)
            period=(obj.ppval.breaks(obj.ppval.pieces +1)-obj.ppval.breaks(1));
            w=2*pi*(1/period);
            for n=1:howmanyn
                for i=1:(obj.ppval.pieces +1)
                    t=obj.ppval.breaks(i);
                    
                   result(n).cos(i,:)=2/period * [ sin(n*t*w)/(n*w) ];                                                                    
                   result(n).sin(i,:)=2/period * [ -cos(n*t*w)/(n*w)];                                                
                
                end
            end
            
            for n=1:howmanyn
                
            resultcos1=result(n).cos(2:(obj.ppval.pieces +1),:);
            result(n).cos=resultcos1-result(n).cos(1:obj.ppval.pieces,:);
            
            resultsin1=result(n).sin(2:(obj.ppval.pieces +1),:);
            result(n).sin=resultsin1-result(n).sin(1:obj.ppval.pieces,:);
            
            end

        end% integral_for_harmonics func  
       function coef=getcoef(obj,howmany_n) 
           coef_of_int=obj.ppval.coefs; %general coefs for integration(stairs)
                %sorting c0
              integralresult=integral_for_harmonics(obj,howmany_n);

            for n=1:howmany_n
               res(n).cos=integralresult(n).cos.*coef_of_int;
               res(n).sin=integralresult(n).sin.*coef_of_int;   
            end

            integralfordc=integral_for_dc(obj);
            coef.dc=sum(sum(integralfordc.*coef_of_int))/2;
            coef.complex = zeros(1,2*howmany_n +1);
            coef.complex(howmany_n+ 1)= coef.dc;

            for n=1:howmany_n
                coef.cos(n) =sum(sum(res(n).cos));
                coef.sin(n) =sum(sum(res(n).sin)); 
                coef.complex(howmany_n + n+1)= coef.cos(n)*2 + coef.sin(n)*-2j;
                coef.complex(howmany_n - n +1 ) = coef.cos(n)*2 - coef.sin(n)*-2j;
            end
      end %getcoefstairs func 
       
   end %methods
    methods (Static)
        function [symintcos, symintsin] = get_symint ()
           syms t;
           syms n;
           syms w;
           myfunc_cos= [ cos(n*w*t)];
          myfunc_sin= [ sin(n*w*t)];
          symintcos = int( myfunc_cos,t);
          symintsin=int(myfunc_sin,t);
       end% get_symint func
    end%static methods
end %Substairs class