function [ss,p,pwmd] = dabest3(data,identifiers,p,pwmd,varargin)
%--------------------------------------------------------------------------
% Data inputs have two valid options:
% Option 1:
%   data:           nx1 vector of values (double, int, etc)
%   identifiers:    nx1 cell array of labels
% Option 2:
%   data:           nxm array of values (double, int, etc)
%   identifiers:    1xm cell array of labels (label m corresponds to the
%                   m-th column in data
%--------------------------------------------------------------------------
% Other inputs:
%   p:              panel for comparison plots against label 1 (see panel.m)
%   pwmd:           panel for pair wise multiple comparison (see panel.m)
%   varargin:       input options (see FscatJit2.m)
%
% Outputs:
%   ss:             structure (or cell array of structure) describing the
%                   statistical analysis outputs
%   p:              panel for comparison plots against label 1
%   pwmd:           panel for pair wise multiple comparison
%--------------------------------------------------------------------------
% note that this is changed from the original dabest2.m as provided by the
% toolbox
% 
% Liangyu Tao June 10, 2019
if size(data,2) ~=1
    identifiers = repmat(identifiers,size(data,1),1);
    data = reshape(data,[],1);
    identifiers = reshape(identifiers,[],1);
end


nVarargs = length(varargin);
if nVarargs == 0
    [ss,p,pwmd] = FscatJit2(identifiers,data,p,pwmd);
elseif nVarargs == 1
    [ss,p,pwmd] = FscatJit2(identifiers,data,p,pwmd,varargin{1});
elseif nVarargs == 2
    [ss,p,pwmd] = FscatJit2(identifiers,data,p,pwmd,varargin{1},varargin{2});
elseif nVarargs == 3
    [ss,p,pwmd] = FscatJit2(identifiers,data,p,pwmd,varargin{1},varargin{2},varargin{3});
elseif nVarargs == 4
    [ss,p,pwmd] = FscatJit2(identifiers,data,p,pwmd,varargin{1},varargin{2},varargin{3},varargin{4});
elseif nVarargs == 5
    [ss,p,pwmd] = FscatJit2(identifiers,data,p,pwmd,varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
end

% avr = [];moes = [];
% % if ~isempty(varargin)
% %     if strcmp(varargin{1},'Paired')
% %         [ss] = FscatJit2(identifiers, data,'Y');
% %     
% %     elseif strcmp(varargin{1},'mergeGroups')
% %         [ss,avr,moes] = FscatJit2_mergeGroups(identifiers, data);
% %     end
% % else
%     [ss] = FscatJit2(identifiers, data,varargin);
% % end

end