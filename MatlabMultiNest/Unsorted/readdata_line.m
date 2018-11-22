function data = readdata_line()
%
% read sample data including outliers (from hogg et al., 1008.4686)
% returning relevant data as a cell array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = load('data_line.txt', '-ascii');

% extract loaded data
ID = temp(:,1);
x = temp(:,2);
y = temp(:,3);
sigmay = temp(:,4);
sigmax = temp(:,5);
rhoxy = temp(:,6);

% package relevant data as a cell array
data{1} = x;
data{2} = y;
data{3} = sigmay;

return
