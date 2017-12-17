function elif_plot_set (FontSize, LineWidth)
                      set( gca , 'units' , 'normalized' )
                      set( gca , 'Box' ,'on' , 'FontName' , 'Cambria' , ...
                      'FontSize' , FontSize , 'FontWeight' , 'bold' , 'LineWidth' , LineWidth )
                      grid on;    
                       set( gcf , 'units' , 'normalized' , 'Position' , [0 0 1 1] );
                      style = hgexport('factorystyle');
                      style.Bounds = 'tight';
                      hgexport( gcf , '.tmpmatlab' , style , 'applystyle' , true );  
end%elif_plot_set