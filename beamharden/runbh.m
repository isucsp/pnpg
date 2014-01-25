
%opt.t3=0;       % set t3 to ignore value of opt.a
opt.numCall=1;
aArray=-6.8:0.2:-6.2;
muLustigArray=logspace(-15,-6,5);
j=3;
opt.muLustig=muLustigArray(j);
opt.skipIe=0;
%string='realCTCircularMaskFullProjASSparseWithResampleMu1000itrRunOnGridforA_2';
%string='castSimCircularMaskFullProjASSparseSparkInputWithResampleMu1000';
%string='twoMaterialsCircularMaskFullProjASWithResampleMu1000itr';
string='castSimCircularMaskFullProj';
if((~isfield(opt,'t3')) || opt.t3~=0 )
    string=[string 'Sparse'];
end
if(spark && ~streq(imageName,'realct'))
    string=[string 'SparkInput'];
end
if(opt.numCall>1)
    string=[string 'WithResampleMu' num2str(opt.numCall) 'Times'];
end
if(opt.skipIe)
    string=[string 'KnownIe'];
end
string=[string num2str(opt.maxItr)];
aArray=-6.5;

if(opt.numCall>1)
    for i=1:length(aArray)
        opt.a=aArray(i);
        out{i}=beamhardenASSparseResampleMu(B,Bt,C,Ct,y,initSig,opt);
        RMSE(i)=1-(out{i}(end).alpha'*opt.trueAlpha/norm(out{i}(end).alpha))^2;
    end
else
    for i=1:length(aArray)
        opt.a=aArray(i);
        out{i}=beamhardenASSparse(B,Bt,C,Ct,y,initSig,opt);
        RMSE(i)=1-(out{i}.alpha'*opt.trueAlpha/norm(out{i}.alpha))^2;
    end
end
if(~ exist('st.mat'))
    save('st.mat','st');
end
clear 'st'
save(string);

return;

for i=length(t3Array)-1 %:-1:6
    for j=3 %1:length(muLustigArray)
        opt.t3=t3Array(i);
        opt.muLustig=muLustigArray(j);
        out(i,j)=beamhardenASSparse(B,Bt,C,Ct,y,initSig,opt);
        RMSE(i,j)=1-(out(i,j).alpha'*trueAlpha/norm(out(i,j).alpha))^2;
        save('realCTRunOnGrid.mat');
    end
end

