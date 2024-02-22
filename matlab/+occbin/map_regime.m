function [regime, regime_start, error_flag]=map_regime(binding_indicator,debug_switch)
% function [regime, regime_start, error_flag]=map_regime(binding_indicator)
% Map regime indicator into information
%
% Inputs:
% - binding_indicator [integer]   [nperiods by 1] vector of regime indices
% - debug_switch      [boolean]   indicator for printing warnings
%
% Outputs:
% - regime          [integer]   [1 by n_transitions] vector of regime number indices
% - regime_start    [integer]   [1 by n_transitions] vectors with period numbers in which regime starts
% - error_flag      [boolean]   1 if regime never leaves 1 or is still there at the end of nperiods
%                               0 otherwise

% Original authors: Luca Guerrieri and Matteo Iacoviello 
% Original file downloaded from:
% https://www.matteoiacoviello.com/research_files/occbin_20140630.zip
% Adapted for Dynare by Dynare Team.
%
% This code is in the public domain and may be used freely.
% However the authors would appreciate acknowledgement of the source by
% citation of any of the following papers:
%
% Luca Guerrieri and Matteo Iacoviello (2015): "OccBin: A toolkit for solving
% dynamic models with occasionally binding constraints easily"
% Journal of Monetary Economics 70, 22-38

error_flag=0;
if isempty(binding_indicator)
    binding_indicator = false;
end
if ischar(binding_indicator)
    % binding indicator represented as a string
    [regime, regime_start] = str2regime(binding_indicator);
    return
end
if any(binding_indicator>1)
    % binding indicator represented a binary (int64)
    binding_indicator = bin2indicator(binding_indicator);
end
% analyze binding_indicator and isolate contiguous periods in the other regime.
regime(1) = binding_indicator(1);
regime_index = 1;
regime_start(1) = 1;
for i=2:length(binding_indicator)
    if binding_indicator(i)~=regime(regime_index)
        regime_index=regime_index+1;
        regime(regime_index) = binding_indicator(i);
        regime_start(regime_index)=i;
    end
end

if (regime(1) == 1 && isscalar(regime_start))
    disp_verbose('map_regime: Binding regime was never left. nperiods needs to be increased.',debug_switch);
    error_flag=1;
end

if (regime(end)==1)
    disp_verbose('map_regime: Constraint(s) are binding at the end of the sample. nperiods needs to be increased.',debug_switch);
    error_flag=1;
end


function binary_indicator = bin2indicator(A)
% function binary_indicator = bin2indicator(A)
% Map regime binary representation into regime indicator
% Inputs:
% - A                 [integer]   scalar whose bits represent regime
%
% Outputs:
% - binding_indicator [integer]   [nperiods by 1] vector of regime indices

binary_indicator = zeros((length(A)-1)*50+length(dec2bin(A(end)))+1,1);
for ka=1:length(A)
    a = dec2bin(A(ka));
    bina = a(end:-1:1);
    for k=1:length(bina)
        binary_indicator(k+50*(ka-1)) = logical(str2double(bina(k)));
    end
end

function [regime, regimestart] = str2regime(regime_string)
% function [regime, regimestart] = str2regime(regime_string)
% Map regime string representation into regime info
%
% Inputs:
% - regime_string  [char array]   string representing regime
%
% Outputs:
% - regime         [array] vector of regimes
% - regimestart    [array] vector of periods where new regime starts

if size(regime_string,1)==1
    regimestart = str2num(regime_string(1,:));
    regime = false(1,length(regimestart));
    for k=length(regimestart)-1:-2:1
        regime(1,k) = true;
    end
else
    regime = logical(str2num(regime_string(1,:)));
    regimestart = str2num(regime_string(2,:));
end

