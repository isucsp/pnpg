function getOthersCode

[a,~]=dbstack('-completenames');
keyboard
a=a(1);
[pathstr,~,~]=fileparts(a.file);
currPath=pwd;
cd(pathstr);

% checkout TFOCS code repository from github
cd others
if(exist('TFOCS','dir'))
    system('git pull');
else
    system('git clone git@github.com:cvxr/TFOCS.git');
end
cd ../


cd(currPath);
end % function getOthersCode
