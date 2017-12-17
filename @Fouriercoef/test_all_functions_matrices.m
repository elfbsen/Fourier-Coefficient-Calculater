 function test_all_functions_matrices
 clear all;
close all;
clear classes;

format long;
 N = 101;
T = 1;
no_periods = 1;
%  try with more than one periods (keep the number odd)

A = 1;
offset = 0;
w = 2 * pi / T;
t = linspace( 0 , no_periods * T , N * no_periods + 1 );

t=t(1:end-1);
%  initial signal
x = A * sin( w * t ) + offset;
x1= A / 4 * cos( w * t) + A / 2 * sin( 4 * w * t ) + offset * 5;
x2 = A / 4 * cos( w * t) + A / 2 * sin( 3.5 * w * t ) + offset * 5;
%D = [t' x'];
D=[ t' x' x1' x2'];
 
        str_interpolation_type = 'spline';
        flag_uniform= 1;
        flag_pointnum_enough =1;
        flag_plotter=1; 
        bigger_pointnum_interpolation = 501; %must be odd
        harmonic_num_for_integralfft = (bigger_pointnum_interpolation-1)/2; 
        fouriercoef_prop_obj=fourier_coef_property (D , str_interpolation_type , flag_uniform , flag_pointnum_enough);
        myfouriercoef_obj=Fouriercoef(D,fouriercoef_prop_obj);
       
        getcoef(myfouriercoef_obj,bigger_pointnum_interpolation,harmonic_num_for_integralfft,flag_plotter);
        dback1=getdataback(myfouriercoef_obj,1,'original',1:1:5);

        derivative_coef(myfouriercoef_obj,1);
        dback2=getdataback(myfouriercoef_obj,1,'derivative');
        title('Derivation');

        integral_coef(myfouriercoef_obj,1);
        dback2=getdataback(myfouriercoef_obj,1,'integral');
        title('Integration');

        timeshift_coef(myfouriercoef_obj, 0.25,1);        
        dback3=getdataback(myfouriercoef_obj,1,'timeshift');
        title('Time Shift');
    end %test_all_functions_matrices