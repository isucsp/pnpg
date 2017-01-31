function setupPath
%
% Author: Renliang Gu (gurenliang@gmail.com)
%
% 0.5: put path adding to the end
% 0.4: add variable cleaning statements
% 0.3: add the current path
% 0.2: implement by dbstack, so that it can be called from any location

[a,~]=dbstack('-completenames');
a=a(1);
[pathstr,~,~]=fileparts(a.file);

cd('rwt')
tgt=['mdwt.' mexext];
if(~exist(tgt,'file')...
        || isOlder(tgt,{'mdwt.c','mdwt_r.c'}))
    mex -largeArrayDims mdwt.c mdwt_r.c
end
tgt=['midwt.' mexext];
if(~exist(tgt,'file')...
        || isOlder(tgt,{'midwt.c','midwt_r.c'}))
    mex -largeArrayDims midwt.c midwt_r.c
end
cd(pathstr)

addpath([pathstr filesep 'rwt']);
addpath([pathstr filesep 'npg']);
addpath([pathstr filesep 'utils']);
addpath([pathstr filesep 'utils' filesep 'L-BFGS-B-C' filesep 'Matlab']);
addpath([pathstr filesep 'irt' filesep 'nufft']);
addpath([pathstr filesep 'irt' filesep 'systems']);
addpath([pathstr filesep 'irt']);
addpath([pathstr filesep 'irt' filesep 'emission']);
addpath([pathstr filesep 'irt' filesep 'general']);
addpath([pathstr filesep 'irt' filesep 'transmission']);
addpath([pathstr filesep 'irt' filesep 'fbp']);
addpath([pathstr filesep 'irt' filesep 'data']);
addpath([pathstr filesep 'irt' filesep 'utilities']);
addpath([pathstr filesep 'others']);
addpath([pathstr filesep 'others' filesep 'FPC_AS']);
addpath([pathstr filesep 'others' filesep 'FPC_AS' filesep 'src']);
addpath([pathstr filesep 'others' filesep 'FPC_AS' filesep 'prob_gen']);
addpath([pathstr filesep 'others' filesep 'FPC_AS' filesep 'prob_gen' filesep 'classes']);
addpath([pathstr filesep 'others' filesep 'glmnet_matlab']);
addpath([pathstr filesep 'others' filesep 'fpc' filesep 'solvers']);
addpath([pathstr filesep 'others' filesep 'fpc' filesep 'solvers' filesep 'utilities']);
if(exist('others/TFOCS','dir'))
    addpath([pathstr filesep 'others' filesep 'TFOCS']);
end

%slCharacterEncoding('UTF-8');
%disp('INFO: if your editor does not show (α,β) properly as $(\alpha,\beta)$ rendered with tex, please run the following command:');
%disp('    slCharacterEncoding(''UTF-8'');');
end

function o = isOlder(f1, f2)
    if(~iscell(f2))
        dependence{1}=f2;
    else
        dependence=f2;
    end
    file=dir(dependence{1});
    time=datenum(file.date);
    for i=2:length(dependence)
        file=dir(dependence{i});
        time=max(time,datenum(file.date));
    end
    file=dir(f1);
    o=datenum(file.date)<time;
end

