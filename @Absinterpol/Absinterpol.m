classdef (Abstract) Absinterpol < handle
    
     properties (GetAccess=public, SetAccess=protected)
        ppval
                
     end
     
   methods (Abstract)
       interpolation_vector
       integral_for_dc
       integral_for_harmonics
       getcoef
   end
   
   methods(Static)
       get_symint
   end
   

end