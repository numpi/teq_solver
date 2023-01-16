function data = generate_latex_histogram2(file_name, ns, tex_file)

data = [];
for i = 1 : length(ns)
    n = ns(i);

    timedata = dlmread(sprintf('%s_n=%d.dat', file_name, n));
    data = [ data , timedata ];
end

fh = fopen(sprintf('%s.tex', file_name), 'w');

for j = 1 : size(data, 1)
    fprintf(fh, '\\addplot coordinates { ');
    for i = 1 : length(ns)
        fprintf(fh, '($n=%d$,%.2f) ', ns(i), data(j,i));
    end
    fprintf(fh, ' };\n');
end

fclose(fh);

return;

% nmin_size			decide what to vary either 'nmin' or 'size'
% n					value of the quantities that do not vary
% ind				index in the table of n (if 'nmin' is the row index in the table, if is 'size' it indicates the 1st,2nd or 3rd group of 4 columns)	
	data = load(strcat('../Results/', file_name, '.dat'));
	if strcmp(nmin_size, 'nmin')
		fileID = fopen(strcat('../Results/', file_name,'n=', num2str(n), '.txt'), 'w');
		for j = 1:4
			fprintf(fileID, '\\addplot coordinates {($n_{\\min}=256$,%.2f) ($n_{\\min}=512$,%.2f) ($n_{\\min}=1024$,%.2f) };\n', data(ind, [j,j+4,j+8]));
		end		
	else
		fileID = fopen(strcat('../Results/', file_name,'nmin=', num2str(n), '.txt'), 'w');
		for j = (ind-1)*4+1:ind*4
			fprintf(fileID, '\\addplot coordinates {($n=2048$,%.2f) ($n=4096$,%.2f) ($n=8192$,%.2f) };\n', data(end-2:end, j));
		end
	end	
	fclose(fileID);
end
