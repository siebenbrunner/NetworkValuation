%% Auxiliary fair conversion function
% Returns a fair conversion function given the parameter gamma. 
%
% *Inputs*
%
%  dblGamma: parameter for fair conversion, as described in Puig and Siebenbrunner 2018. 

%
%  *Outputs* 
%
%  Conversion function
% 
% Author: Christoph Siebenbrunner and Matias Puig
% Last modified: 03.09.2018
%

function funConversion = fairConversionFun(dblGamma) 
    funConversion = @(x,y,z) fairConversion(x,dblGamma,y,z);
end