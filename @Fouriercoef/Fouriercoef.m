classdef Fouriercoef < handle
properties (GetAccess=public, SetAccess=protected) 
      interclass
      data
      time
      old_time = [];
      no_inputs
      fp %fouriercoef prop class
      coef
end    %prop    

methods (Static)
    test_all_functions_matrices; 
end% methods static
methods             
      function obj=Fouriercoef(d,fcoef_prop) %d numeric, data and time together,    

        obj.fp=fcoef_prop;
                       
            if nargin>0 %nargout varargin         
              obj.no_inputs=size(d,2)-1; 
              obj.time = d(:,1);
             %%"For" statement for datas in the given matris%%
             for kk=1:1:obj.no_inputs
                 
             obj.data{kk}= d(:,kk+1); 
             
              if obj.fp.flag_isuniform_==1 && obj.fp.flag_ispointnumenough_==1                       
                        
              else
                                  
                  switch lower(obj.fp.str_interpoltype_)
                     
                      case 'spline'
                          obj.interclass{kk}=Subspline(d(:,kk+1),d(:,1));
                      case 'linear'
                          obj.interclass{kk}=Sublinear(d(:,kk+1),d(:,1));
                      case 'stairs'
                          obj.interclass{kk}=Substairs(d(:,kk+1),d(:,1));
                      otherwise
                          warning('Interpolation type is not given properly, it is taken linear automatically');
                          obj.interclass{kk}=Sublinear(d(:,kk+1),d(:,1));
                          obj.fp.str_interpoltype_='linear';
                  end
              end
             end
            end
      end%constructor
      function getcoef(obj, pointnum, howmany_n, flag_plotter)
        for kk=1:1:obj.no_inputs
        switch obj.fp.flag_isuniform_             
            case 1 %flag_isuniform=1 
                           
            switch obj.fp.flag_ispointnumenough_               
               case 0      %f2=0 point number is not enough, interpolation is needed 
                tmp_vec = interpolation_vector(obj.interclass{kk}, pointnum +1);
                xfft= fft(tmp_vec(1:end-1))/pointnum ;    %fft
                xfft= [ xfft(((pointnum-1)/2 + 2) : pointnum)  xfft(1: ((pointnum - 1)/2 +1))] ;%fftshift
                obj.coef{kk}.complex = xfft;
                if kk==1
                obj.old_time = obj.time;
                obj.time = linspace(obj.interclass{kk}.ppval.breaks(1),obj.interclass{kk}.ppval.breaks(end),pointnum+1);
                obj.time = obj.time(1:end-1);
                end
                if flag_plotter
                  figure;                     
                              subplot(2,1,1), stem(-(pointnum - 1)/2 : (pointnum - 1)/2 , real(xfft)*2,'bo-','filled');title(sprintf ('Real Part of Shifted Fourier Coefficients - node %d  ',kk)),grid on; 
                              subplot(2,1,2), stem(-(pointnum - 1)/2 : (pointnum - 1)/2 ,imag(xfft)*-2,'bo-','filled');title(sprintf ('Imaginary Part of Shifted Fourier Coefficients - node %d  ',kk)), grid on;
                               elif_plot_set(22,3);
                         grid on;
                end
                xfft=xfft(1,(pointnum-1)/2 +1 : pointnum ); %positive freqs

               obj.coef{kk}.dc = real(xfft(1));
               obj.coef{kk}.cos= real(xfft(2:end)) *2;
               obj.coef{kk}.sin = imag(xfft(2:end)) *-2;
               
           
              
               case 1 %f2=1, point number is enough,  interpolation is not needed.  
                xfft= fft( obj.data{kk}')/(numel( obj.data{kk}));    %fft    
                pointnum_in=size(xfft, 2)
                xfft=[ xfft(((pointnum_in-1)/2 +2):pointnum_in) xfft(1:((pointnum_in-1)/2 +1))]; %fft shift  
                obj.coef{kk}.complex = xfft;
                if flag_plotter
                  figure;   
                   elif_plot_set(22,3);
                              subplot(2,1,1), stem(-(pointnum_in-1)/2 :(pointnum_in-1)/2 ,real(xfft)*2,'bo-','filled','LineWidth', 4);title(sprintf ('Real Part of Shifted Fourier Coefficients - node %d  ',kk),'FontSize',22),grid on; 
                              elif_plot_set(22,3);
                              subplot(2,1,2), stem(-(pointnum_in-1)/2 :(pointnum_in-1)/2 ,imag(xfft)*-2,'bo-','filled','LineWidth',4);title(sprintf ('Imaginary Part of Shifted Fourier Coefficients - node %d  ',kk),'FontSize',22),grid on;
                              elif_plot_set(22,3);
                         grid on;
                end   
                xfft=xfft((pointnum_in-1)/2 +1 : pointnum_in ); %positive freqs

               obj.coef{kk}.dc = real(xfft(1));
               obj.coef{kk}.cos= real(xfft(2:end)) *2;
               obj.coef{kk}.sin = imag(xfft(2:end)) *-2;
               
             
            end   
            
            case 0  %f1=0, data is nonuniform, interpolation is exact                
               switch obj.fp.flag_ispointnumenough_              
                 case 0      %f2=0 point number is not enough, interpolation is needed
                     if kk==1
                 obj.old_time= obj.time;
                 obj.time = linspace(obj.interclass{kk}.ppval.breaks(1),obj.interclass{kk}.ppval.breaks(end),pointnum);
                     end
             switch lower(obj.fp.str_interpoltype_)
                 case 'spline'
                     obj.interclass{kk}=Subspline(interpolation_vector(obj.interclass{kk},pointnum),...
                          obj.time);
                          obj.coef{kk}=getcoef(obj.interclass{kk},howmany_n);
                 case 'linear'
                        obj.interclass{kk}=Sublinear(interpolation_vector(obj.interclass{kk},pointnum),...
                          obj.time);
                          obj.coef{kk}=getcoef(obj.interclass{kk},howmany_n);
                 case 'stairs'
                          obj.interclass{kk}=Substairs(interpolation_vector(obj.interclass{kk},pointnum),...
                          obj.time);
                          obj.coef{kk}=getcoef(obj.interclass{kk},howmany_n);
             end    
      
                          
        case 1 %f2=1, point number is enough,  interpolation is not needed.
             if kk==1
                 obj.old_time= obj.time;
                 obj.time = linspace(obj.interclass{kk}.ppval.breaks(1),obj.interclass{kk}.ppval.breaks(end),howmany_n*2 +1);
                     end
                switch lower(obj.fp.str_interpoltype_)                
                        
                 case 'spline'
                     obj.coef{kk}=getcoef(obj.interclass{kk},howmany_n);
                 case 'linear'
                     obj.coef{kk}=getcoef(obj.interclass{kk},howmany_n);
                 case 'stairs'
                     obj.coef{kk}=getcoef(obj.interclass{kk},howmany_n);
                end
                  %plotting coefs  
        if flag_plotter
            %bradaki katsayılar organik değil bişi yapmak lazım??
            plotting_cos=[ obj.coef{kk}.cos(end:-1:1) obj.coef{kk}.dc obj.coef{kk}.cos ];
            plotting_sin=[ -obj.coef{kk}.sin(end:-1:1) 0  obj.coef{kk}.sin ];
                          figure;                     
                              subplot(2,1,1), stem(-howmany_n:1:howmany_n,plotting_cos,'bo-','filled');title(sprintf ('cos coef - data %d  ',kk)),grid on; 
                              subplot(2,1,2), stem(-howmany_n:1:howmany_n,plotting_sin,'bo-','filled');title(sprintf ('sin coef - data %d  ',kk)),grid on;
                              set(gca,'Box','on','FontName','Arial',...
                         'FontSize',22,'FontWeight','bold','LineWidth',3)
                         set( gcf , 'units' , 'normalized' , 'Position' , [0 0 1 1] );
                          style = hgexport('factorystyle');
                            style.Bounds = 'tight';
                          hgexport( gcf , '.tmpmatlab' , style , 'applystyle' , true );
                         grid on;
        end    
               end
               
        end
        
      
        end
      end   %getcoef function
      function dback=getdataback(obj,data_no, coef_type,plot_list) 
          if numel(obj.old_time)==0
              obj.old_time=obj.time
          end
          dback= zeros(data_no, numel(obj.time));
          kk=data_no;
               freq=get_freq(obj,kk);
               coef_str=coef_type;
               mycoef.dc=obj.coef{kk}.dc;        
               switch lower(coef_str)
                   case 'original'
                       mycoef.cos=obj.coef{kk}.cos;
                       mycoef.sin = obj.coef{kk}.sin;
                       mycoef.complex = obj.coef{kk}.complex;
                   case 'integral'  
                       mycoef.cos=obj.coef{kk}.integralcos;
                       mycoef.sin = obj.coef{kk}.integralsin;
                       mycoef.dc =obj.coef{kk}.integraldc;
                   case 'derivative'
                       mycoef.cos=obj.coef{kk}.derivativecos;
                       mycoef.sin = obj.coef{kk}.derivativesin;
                   case 'timeshift'
                       mycoef.cos=obj.coef{kk}.timeshiftcos;
                       mycoef.sin = obj.coef{kk}.timeshiftsin;
                        case 'low_pass_filtered'
                       mycoef.cos=obj.coef{kk}.lpfcos;
                       mycoef.sin = obj.coef{kk}.lpfsin;
                   otherwise
                       warning...
                          ('coef_type must be "integral", "derivative", "timeshift","low_pass_filtered" or "original". Given coef type is not one of them and coef_type chosen original automatically.');
                       mycoef.cos=obj.coef{kk}.cos;
                       mycoef.sin = obj.coef{kk}.sin;
               end
                for n=1:numel(mycoef.cos)
                         cosx(n,:)=cos (2*pi*n*freq*obj.time);
                         sinx(n,:)=sin(2*pi*n*freq*obj.time);
                end   
                                
               switch nargin               
               case 3
               dback=mycoef.dc + mycoef.cos*cosx + mycoef.sin*sinx;
               figure; plot(obj.time, dback, 'b',...
                          obj.old_time,obj.data{kk}, 'r','LineWidth',3), ...
                          set(gca,'Box','on','FontName','Arial',...
                         'FontSize',22,'FontWeight','bold','LineWidth',3)
                         set( gcf , 'units' , 'normalized' , 'Position' , [0 0 1 1] );
                          style = hgexport('factorystyle');
                            style.Bounds = 'tight';
                          hgexport( gcf , '.tmpmatlab' , style , 'applystyle' , true );
                         grid on;
                        % title(sprintf(' - data %d',kk));    
                        legend('\fontsize{16} Processed Signal ','\fontsize{16} Original Signal', 'Location', 'southwest');
                case 4            
                if max(plot_list)>numel(mycoef.cos)
                    warning('Wanted max harmonic is bigger than the biggest harmonic of coefs. There are figures can not plotted.');
                    for i=1:numel(plot_list)
                        if max(plot_list)>numel(mycoef.cos);
                            plot_list=plot_list(plot_list~=max(plot_list));
                        else
                            break;
                        end
                    end
                end
                
                b=numel(plot_list);
                if numel(plot_list)>4
                    plot_list=plot_list(1:4);
                    b=4;
                    warning('Only first 4 figure are plotted.');
                end
                dback = zeros( plot_list(end), numel(obj.time));              
                for n=1:plot_list(end)
                dback(n,:)=mycoef.dc + mycoef.cos(1:n)*cosx(1:n,:) + mycoef.sin(1:n)*sinx(1:n,:);
                end    
             
 
          
          figure
          for i=1:4      
              if (b-i)>=0
              switch mod(i,2)
                  case 1
                     subplot(2,2,i), plot(obj.time, dback(plot_list(i),:),'b', ...
                         obj.old_time,obj.data{kk}, 'r','LineWidth',3), ...
%                          axis([obj.time(1) obj.time(end)  (min(dback(plot_list(i),:))-abs(min(dback(plot_list(i),:))*0.15)) ...
%                         (max(dback(plot_list(i),:)) + abs(max(dback(plot_list(i),:))*0.15))]),
%                          set(gca,'units','normalized')
                         set(gca,'Box','on','FontName','Arial',...
                         'FontSize',22,'FontWeight','bold','LineWidth',3)
                         set( gcf , 'units' , 'normalized' , 'Position' , [0 0 1 1] );
                          style = hgexport('factorystyle');
                            style.Bounds = 'tight';
                          hgexport( gcf , '.tmpmatlab' , style , 'applystyle' , true );
                         grid on;
                         title(sprintf('%d harmonic used (blue) - signal @node%d ',plot_list(i),kk));
                  case 0
                     subplot(2,2,i), plot(obj.time, dback(plot_list(i),:),'b',...
                         obj.old_time,obj.data{kk}, 'r','LineWidth',3), ...
%                          axis([obj.time(1) obj.time(end)  (min(dback(plot_list(i),:))-abs(min(dback(plot_list(i),:))*0.15)) ...
%                         (max(dback(plot_list(i),:)) + abs(max(dback(plot_list(i),:))*0.15))]),
                    set(gca,'units','normalized')
                    set(gca,'Box','on','FontName','Arial',...
                    'FontSize',22,'FontWeight','bold','LineWidth',3)
                set( gcf , 'units' , 'normalized' , 'Position' , [0 0 1 1] );
                style = hgexport('factorystyle');
                style.Bounds = 'tight';
                hgexport( gcf , '.tmpmatlab' , style , 'applystyle' , true );
                grid on;
                         title(sprintf('%d harmonic used (blue) - signal @node%d ',plot_list(i),kk));  
              end
              else
                  break;
              end
                end
            end
       
  end %getdataback function
      function freq=get_freq(obj,kk)         
               freq=1/(obj.time(end) - obj.time(1));       
      end %getfreq function
end %methods
end  %Fouriercoef class