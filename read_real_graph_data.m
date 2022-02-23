% generate adjacency matrix from .txt real network data

function [A] = read_real_graph_data(realGraphName)

%% read network data
currentFolder = pwd;
dataPath = strcat(currentFolder,'/graphs/',realGraphName,'.txt');
network_data = importdata(dataPath);

text = network_data.textdata;
len = strlength(text{1});
n = str2double(text{1}(3:len));

%% generata adjacency matrix
numEntry = size(network_data.data,1);
connectivity = network_data.data;
A = zeros(n);
for k = 1 : numEntry
    i = connectivity(k, 1) + 1;
    j = connectivity(k, 2) + 1;
    A(i,j) = connectivity(k, 3);
end
A = A + A';
