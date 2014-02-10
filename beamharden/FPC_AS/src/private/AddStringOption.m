function opts = AddStringOption(opts, name, default, list, memo)

%-------------------------------------------------------------------------
% set fields of "opts" which are strings
%
% opts:     input/output structure
% name:     a field to add
% default:  default value of this field
% list:     a list of valid values
% memo:     a description
%
% add a field "name" to opts, name must be in [low, up]
% memo is a string which describes the field "name"
% if "name" exists, check whether it is valid. If it is not valid, set it to the default value
%-------------------------------------------------------------------------

switch default
    case list
    otherwise
        fprintf('\ndefault value: %s, the available options for %s are:\n', default, name);
        for dj = 1: length(list); fprintf('%s \t', list{dj}); end; fprintf('\n');
        error('option error');
end

if isfield(opts, name)
    val = opts.(name);
    switch val
        case list
        otherwise
            fprintf('\ncurrent value: %s, the available options for %s are:\n', val, name);
            for dj = 1: length(list); fprintf('%s \t', list{dj}); end;
            fprintf('\nSet this option to the default value: %s\n\n', default);
            opts.(name) = default;
    end
else
    opts.(name) = default;
end
opts.memo.(name) = memo;


% if isfield(opts, name)
%     val = getfield(opts,name);
%     switch val
%         case list
%         otherwise
%             fprintf('\ncurrent value: %s, the available options for %s are:\n', val, name);
%             for dj = 1: length(list); fprintf('%s \t', list{dj}); end;
%             fprintf('\nSet this option to the default value: %s\n\n', default);
%             opts = setfield(opts, name, default);
%     end
% else
%     opts = setfield(opts, name, default);
% end
% opts.memo = setfield(opts.memo, name, memo);
