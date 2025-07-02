function Supdate = setStructDefaults(Sinput,Sdefault,warnExtraInput)
% Compares input struct Sinput to Sdefault. Sets any missing fields in S to
% the defaults in Sdefault.
%
% Similar to 'setstructfields', but does not overwrite preexisting fields
% in Sinput.
%
% Gives warning (by default) if input has fields that differ from fields in
% default, unless warnExtraInput set to false.
%
% Note that fields are case-sensitive
%
% example:
%     Sdefault = struct('a',true,'b','two','c',3,'d',[1 0 0]);
%     Sinput = struct('a',false,'B','too','d',[1 1 1],'ee',[]);
%     Supdate = setStructDefaults(Sinput,Sdefault)
%
% has the output:
%     Warning: Input struct has the following fields not in default struct:
%     B ee 
%     > In setStructDefaults (line 51) 
% 
%     Supdate = 
% 
%       struct with fields:
% 
%          a: 0
%          B: 'too'
%          d: [1 1 1]
%         ee: []
%          b: 'two'
%          c: 3



if nargin < 3
    warnExtraInput = true;
end

Supdate = Sinput;

fieldsDefault = fieldnames(Sdefault);
fieldsInput = fieldnames(Sinput);

fieldsNotSet = setdiff(fieldsDefault,fieldsInput);

for i = 1:length(fieldsNotSet)
    curField = fieldsNotSet{i};
    Supdate.(curField) = Sdefault.(curField);
end

if warnExtraInput
    fieldInputsDiffer = setdiff(fieldsInput,fieldsDefault);

if ~isempty(fieldInputsDiffer) 
    
   warning(['Input struct has the following fields not in default struct:',...
       char(10),strjoin(fieldInputsDiffer,' ')]);
end

end