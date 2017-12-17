classdef Subspline < Absinterpol   
    methods 
       function obj = Subspline(d,t)
            if nargin > 0                 
                if isnumeric(d)
                     if isnumeric(t)
                            obj.ppval=spline(t,d);
                     end
                end
            end
       end %constructor       
       function intervec=interpolation_vector(obj,pointnum);
          
           intervec=ppval(obj.ppval,...
               linspace(obj.ppval.breaks(1),...
               obj.ppval.breaks(end),pointnum));
           
       end %func interpolation_vector       
       function result=integral_for_dc(obj)
            period=(obj.ppval.breaks(end)-obj.ppval.breaks(1));
            result=2/period * [obj.ppval.breaks(:).^4/4 ...
                        obj.ppval.breaks(:).^3/3  ...
                        obj.ppval.breaks(:).^2/2 ...
                        obj.ppval.breaks(:).^1 ];
            result1=result(2:(obj.ppval.pieces+1),:);
            result=result1-result(1:obj.ppval.pieces,:);
        end %integral_for_dc func
       function result=integral_for_harmonics(obj,howmanyn)
           period=(obj.ppval.breaks(obj.ppval.pieces +1)-obj.ppval.breaks(1));
            w=2*pi*(1/period);
            
            for n=1:howmanyn
                for i=1:(obj.ppval.pieces +1)
                    t=obj.ppval.breaks(i);
            
                    
                   result(n).cos(i,:)=2/period * [ -(6*cos(n*t*w) - n^3*t^3*w^3*sin(n*t*w) + 6*n*t*w*sin(n*t*w) - 3*n^2*t^2*w^2*cos(n*t*w))/(n^4*w^4) ... 
                                   (n^2*t^2*w^2*sin(n*t*w) - 2*sin(n*t*w) + 2*n*t*w*cos(n*t*w))/(n^3*w^3) ...                   
                                                        (cos(n*t*w) + n*t*w*sin(n*t*w))/(n^2*w^2) ...                                           
                                                                                   sin(n*t*w)/(n*w) ];                                                                    
                   result(n).sin(i,:)=2/period * [ cos(n*t*w)*((6*t)/(n^3*w^3) - t^3/(n*w)) - sin(n*t*w)*(6/(n^4*w^4) - (3*t^2)/(n^2*w^2)) ...
                                                 cos(n*t*w)*(2/(n^3*w^3) - t^2/(n*w)) + (2*t*sin(n*t*w))/(n^2*w^2) ...
                                                  (sin(n*t*w) - n*t*w*cos(n*t*w))/(n^2*w^2) ...
                                                  -cos(n*t*w)/(n*w)];
                                                
                
                end
            end
            
            for n=1:howmanyn
                
            resultcos1=result(n).cos(2:(obj.ppval.pieces +1),:);
            result(n).cos=resultcos1-result(n).cos(1:obj.ppval.pieces,:);
            
            resultsin1=result(n).sin(2:(obj.ppval.pieces +1),:);
            result(n).sin=resultsin1-result(n).sin(1:obj.ppval.pieces,:);
            end
            
            

        end% integral_for_harmonics func
       function coef=getcoef(obj, howmany_n)
                    binom=[1 0 0 0;1 -1 0 0; 1 -2 1 0; 1 -3 3 -1];
  
          %%% COEFs for SPLINE %%%%
            for i=1:(obj.ppval.pieces)
                a=[(obj.ppval.breaks(i)).^0 ...
                (obj.ppval.breaks(i)).^1 (obj.ppval.breaks(i)).^2 ...
                (obj.ppval.breaks(i)).^3];
            
                %t^0
                coefs.t0(i,:)=a.*binom(1,:).*obj.ppval.coefs(i,4);
                %t^1
                coefs.t1(i,:)=a.*binom(2,:).*obj.ppval.coefs(i,3);
                %t^2
                coefs.t2(i,:)=a.*binom(3,:).*obj.ppval.coefs(i,2);
                %t^3
                coefs.t3(i,:)=a.*binom(4,:).*obj.ppval.coefs(i,1);
          
            end
                zeromat=zeros(obj.ppval.pieces, 1);                 
                coefs.t0=[ zeromat zeromat zeromat coefs.t0(:,1)];
                coefs.t1=[zeromat zeromat coefs.t1(:,1) coefs.t1(:,2)];
                coefs.t2=[zeromat coefs.t2(:,1) coefs.t2(:,2) coefs.t2(:,3)];
    
     coef_of_int=coefs.t3+coefs.t2+coefs.t1+coefs.t0;      %general coefs for integration(spline)
     %sorting c3 c2 c1 c0
     
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
       end % getcoefspline function
    end %methods
    methods (Static)
        function [symintcos, symintsin] = get_symint ()
           syms t;
           syms n;
           syms w;
           myfunc_cos= [ cos(n*w*t)*t^3;
                                     cos(n*w*t)*t^2;
                                     cos(n*w*t)*t;
                                     cos(n*w*t)];
          myfunc_sin= [ sin(n*w*t)*t^3;
                                    sin(n*w*t)*t^2;
                                    sin(n*w*t)*t;
                                    sin(n*w*t)];
          symintcos = int( myfunc_cos,t);
          symintsin=int(myfunc_sin,t);
       end% get_symint func
    end%static methods
end % Subspline class

