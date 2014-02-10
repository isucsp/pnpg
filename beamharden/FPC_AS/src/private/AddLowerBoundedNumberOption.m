function opts = AddLowerBoundedNumberOption(opts, name, default, low, up, memo)
%-------------------------------------------------------------------------
% set fields of "opts" which have real values
%
% opts:     input/output structure
% name:     a field to add
% default:  default value of this field
% low:      low bound
% up:       upper bound
% memo:     a description
%
% add a field "name" to opts, name must be in [low, up]
% memo is a string which describes the field "name"
% if "name" exists, check whether it is valid. If it is not valid, set it to the default value
%-------------------------------------------------------------------------

if low > up; error(sprintf('option.%s, low= %e > up=%e\n',name, low, up)); end
if default < low | default > up;
    error(sprintf('option.%s, default= %g, not in [%g, %g]\n',name, default, low, up));
end
if isfield(opts, name)
    val = opts.(name);
    if val < low | val > up;
        warning(sprintf('option: %s, val= %g, not in [%g, %g], set it to default: %g\n',name, val, low, up, default));
        opts.(name) = default;
    end
else
     opts.(name) = default;
end
opts.memo.(name) = memo;

%fprintf('\\item ''%s'': %s \\\\ default: %g, vaild range: [%g, %g]\n \n', strrep(name,'_','\_'), strrep(memo,'_','\_'), default, low, up);

% if low > up; error(sprintf('option.%s, low= %e > up=%e\n',name, low, up)); end
% if default < low | default > up;
%     error(sprintf('option.%s, default= %g, not in [%g, %g]\n',name, default, low, up));
% end
% if isfield(opts, name)
%     val = getfield(opts,name);
%     if val < low | val > up;
%         warning(sprintf('option: %s, val= %g, not in [%g, %g], set it to default: %g\n',name, val, low, up, default));
%         opts = setfield(opts, name, default);
%     end
% else
%     opts = setfield(opts, name, default);
% end
% opts.memo = setfield(opts.memo, name, memo);
