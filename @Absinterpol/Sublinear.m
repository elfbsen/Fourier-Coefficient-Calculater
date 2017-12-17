classdef Sublinear < Absinterpol

   
   methods
       function obj = Sublinear(d,t)
            if nargin > 0          
                if isnumeric(d)
                     if isnumeric(t)                     
                        obj.ppval=interp1(t,d,'linear','pp');                             
                     end
                end     
            end
       end %constructor       
       function intervec=interpolation_vector(obj,pointnum)
              
           intervec=ppval(obj.ppval,...
               linspace(obj.ppval.breaks(1),...
               obj.ppval.breaks(end),pointnum));
           
       end %interpolation func
       function result=integral_for_dc(obj)
            period=(obj.ppval.breaks(end)-obj.ppval.breaks(1));
            result=2/period * [obj.ppval.breaks(:).^2/2 ...
                        obj.ppval.breaks(:).^1 ];
            result1=result(2:end,:);
            result=result1-result(1:end-1,:);
        end %integral_for_dc func
       function result=integral_for_harmonics(obj,howmanyn)
            period=(obj.ppval.breaks(end)-obj.ppval.breaks(1));
            w=2*pi*(1/period);
            for n=1:howmanyn
                for i=1:(obj.ppval.pieces +1)
                    t=obj.ppval.breaks(i);
                    
                   result(n).cos(i,:)=2/period * [ (cos(n*t*w) + n*t*w*sin(n*t*w))/(n^2*w^2) ...                                           
                                                                                   sin(n*t*w)/(n*w) ];                                                                    
                   result(n).sin(i,:)=2/period * [ (sin(n*t*w) - n*t*w*cos(n*t*w))/(n^2*w^2) ...
                                                  -cos(n*t*w)/(n*w)];                                                
                
                end
            end
            
            for n=1:howmanyn       
            result(n).cos=result(n).cos(2:end,:)...
                -result(n).cos(1:end-1,:);  
            result(n).sin=result(n).sin(2:end,:)...
                -result(n).sin(1:end-1,:);
            end
        end% integral_for_harmonics func
       function  coef=getcoef(obj,howmany_n)

      coefs.t0=-(obj.ppval.breaks(1:obj.ppval.pieces))'.*obj.ppval.coefs(:,1)...
          +obj.ppval.coefs(:,2);
      coefs.t1=obj.ppval.coefs(:,1);
 
      coef_of_int=[coefs.t1  coefs.t0];
       
      
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

    coef.cos(n)=sum(sum(res(n).cos));
    coef.sin(n)=sum(sum(res(n).sin));
    coef.complex(howmany_n + n+1)= coef.cos(n)/2 + coef.sin(n)*j/-2;
    coef.complex(howmany_n - n +1 ) = coef.cos(n)/2 - coef.sin(n)*j/-2;
    
end


end %getcoeflinear func   
 
   end %methods
    methods(Static)
        function [symintcos, symintsin] = get_symint ()
           syms t;
           syms n;
           syms w;
           myfunc_cos= [           cos(n*w*t)*t;
                                     cos(n*w*t)];
          myfunc_sin= [          sin(n*w*t)*t;
                                    sin(n*w*t)];
          symintcos = int( myfunc_cos,t);
          symintsin=int(myfunc_sin,t);
          
       end% get_symint func
    end%static methods
end %Sublinear classdef