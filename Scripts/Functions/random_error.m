function [stdev_f,tau_f] = random_error(f,n,freq,z,U_bar,delta_min,delta_max,nfilters,use_defaults)
%Estimates random error using the filtering method of Salesky, Chamecki, and Dias (BLM, 2012)
%----------------------------------------------------------------------------------------------------------
% random_error.m 
% Supplementary Material for:
% "Estimating random error in eddy covariance fluxes and other turbulence statistics: the filtering method"
% Scott T. Salesky, Marcelo Chamecki*, and Nelson L. Dias
% Boundary-Layer Meteorology, 2012
% *Corresponding author. The Pennsylvania State University. chamecki@meteo.psu.edu
%----------------------------------------------------------------------------------------------------------
%This is free software released into the public domain.
%
%Anyone is free to copy, modify, publish, use, compile, sell, or
%distribute this software, either in source code form or as a compiled
%binary, for any purpose, commercial or non-commercial, and by any
%means.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
%IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
%OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
%ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%OTHER DEALINGS IN THE SOFTWARE.
%----------------------------------------------------------------------------------------------------------
% Please refer to the Appendix of article for an explanation of the algorithm applied in this code.
% Before calling this routine:
% -> Apply data selection criteria and preprocessing to data
% -> Generate timeseries of instantaneous flux, e.g. f = w'c' for flux, or f = u'u'u' for third-order moment
%----------------------------------------------------------------------------------------------------------
% Input parameters:
% f                 -> Timeseries of the variable of interest (e.g. f = w'c'(t))
% n                 -> Number of points in f
% freq              -> Sampling frequency (Hz)
% z                 -> Measurement height (m)
% U_bar             -> Mean wind speed (m s^-1)
% nfilters          -> Number of filter widths between delta_min and delta_max
% use_defaults      -> Flag for using (1) or not using (0) default values of delta_min, delta_max, and nfilters
%                       -> if use_defaults = 0, user must specify delta_min, delta_max
%                       -> if use_defaults = 1, code uses default values of delta_min = 10.*z/U_bar, delta_max = T/10, nfilters = 50
% delta_min         -> Minimum filter width (s)
% delta_max         -> Maximum filter width (s)
%
% Output parameters:
% stdev_f           -> Error in f, for averaging period T (same units as f)
% tau_f             -> Estimate of integral timescale of f from filtering method (s)
%----------------------------------------------------------------------------------------------------------

% ---------------------------------------------------------------------------------------------------------
% Begin code ----------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------  

%Check whether to use default values of delta_min,delta_max, and nf
if use_defaults == 0    %Use user-specified values
    f_min = delta_min;
    f_max = delta_max;
    nf = nfilters;
elseif use_defaults == 1    %Use default values
    f_min = 10.*z/U_bar;
    f_max = (n/freq)/10.;
    nf = 50;
else
    error('Invalid argument for use_defaults')
end

%Set up log-space array of filter widths as function of delta_t
delta_t = logspace(log10(f_min),log10(f_max),nf);

%Set up array of standard deviation of filtered flux as function of delta_t
sigma_f = zeros(1,nf);

%Calculate streamwise spacing between points
dx = U_bar/freq;

%Begin loop over filter widths -------------------------------------
for k=1:nf

    %Filter f at scale delta_t
    %----------------------------------------------------------------------
    %Wavenumber increment
    dk = (2.*pi)/(dx*n);
    
    %Spatial filter via Taylor's hypothesis
    delta_x = delta_t(k)*U_bar;
    
    %Assemble transfer function of the box filter
    for i=2:(n/2)
        filter(i) = sin((i-1)*dk*delta_x/2.0)/((i-1)*dk*delta_x/2.0);
    end
    
    %Forward FFT
    f_fft = fft(f);
    
    %Multiply FFT(f) by transfer function
    for i=2:(n/2)
        f_fft(i) = f_fft(i)*filter(i);
        f_fft(n-i+2) = f_fft(n-i+2)*filter(i);
    end 
       
    %Inverse FFT
    f_filtered = ifft(f_fft);
    %----------------------------------------------------------------------
    
	%Calculate standard deviation of f, and store
	sigma_f(k) = std(f_filtered,0,1);	

end

%Fit powerlaw to standard deviation of local filtered f
for i=1:nf
    A(i) = delta_t(i).^(-0.5);
end

C_wc = sigma_f/A;

%Evaluate powerlaw fit for delta_t = T to obtain error estimate
stdev_f = C_wc*((n/freq)^(-1/2));

%Calculate variance of *unfiltered* series of f
var_flux = var(f);

%Estimate integral timescale of f from fitted powerlaw coefficient (Eq. 25)
tau_f = (C_wc^2)/(2.*var_flux);
