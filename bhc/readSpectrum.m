function [epsilon,iota] = readSpectrum(material,kev,rvr)
    % the spectrums are from the SIEMENS X-ray Spectra tool online:
    % https://w9.siemens.com/cms/oemproducts/Home/X-rayToolbox/spektrum/Pages/Default.aspx
    % where we always set Air kerma=1Gy without additional filters
    % and 5% relative voltage ripple.

    filename=lower(material);
    filename=[filename '_' num2str(kev) 'kV_' num2str(rvr) '_1Gy.data'];
    delimiter = ' ';
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, inf,...
        'Delimiter', delimiter, 'MultipleDelimsAsOne', true,...
        'HeaderLines',0, 'ReturnOnError', false);

    fclose(fileID);
    epsilon = dataArray{:, 1};
    iota = dataArray{:, 2};
end

