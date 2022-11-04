% utility script to convert (dis)connection matrices in the output of 
% the LQT toolbox from .mat to .csv; this makes the files compatible with
% the BLDI toolkit
clear
clc

% provide the path to a folder that contains the .mat files; make sure that
% no other mat files are contained in this folder
folder = 'D:\Arbeit_Bern\Projekt_Bayesian_LSM\Disconnection_Prospective';
fileList = struct2cell(dir(fullfile(folder, '*.mat')));



for i=1:size(fileList,2)
    path = strcat(folder, '\', fileList{1,i});

    temp=load(strcat(fileList{2,1}, '\', fileList{1,i}));
    % load opens the mat file as a struct; the next line accesses the first
    % element of this struct and opens it as a matrix
    conn_mat = temp.(subsref(fieldnames(temp),substruct('{}',{1})));

    % replace the extension in the path
    output_filename = strrep(path,'.mat','.csv');
    
    % write as csv with space as delimiter
    writematrix(conn_mat,output_filename,'Delimiter','space')
end
