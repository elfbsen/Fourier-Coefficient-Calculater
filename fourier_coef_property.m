classdef fourier_coef_property < handle
    properties
        flag_isuniform_ 
        flag_ispointnumenough_ 
        str_interpoltype_ 
    end%properties
    
     methods
            function obj=fourier_coef_property (cdata,cstr_interpoltype,cflag_isuniform, cflag_ispointenough)
                   
                if nargin<4
                    warning('number of input arguments are not enough, they re taking values automatically');
                    obj.str_interpoltype_='linear';
                    obj.flag_ispointnumenough_=1;
                    obj.flag_isuniform_=0;
                end
             if ~isnumeric(cdata)
                    warning('data must be numeric');
                end
                
                if (cflag_ispointenough==1)|(cflag_ispointenough==0)
                    obj.flag_ispointnumenough_=cflag_ispointenough;
                else
                    warning('flag_ispointenough must be 0 or 1, this flag is taken 1 automatically');
                    obj.flag_ispointnumenough_=1;
                end
                
                  if (cflag_isuniform==1)|(cflag_isuniform==0)
                    obj.flag_isuniform_=cflag_isuniform;
                else
                    warning('flag_isuniform must be 0 or 1, this flag is taken 0 automatically');
                    obj.flag_isuniform_=0;
                  end
                
                if isstr(cstr_interpoltype)
                    obj.str_interpoltype_=cstr_interpoltype;
                else
                    warning('str_interpoltype must be string, this str choosen "linear" automatically');
                    obj.str_interpoltype_='linear';
                end                            
            end %const
        end%methods
end %fcoefprop 
            
