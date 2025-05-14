function [] = printpdf (filename)

% give user a chance to close file in case it is open

try
	print ('-dpdf', filename);
catch
	fprint ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n CLOSE %s, then hit enter so that the pdf can be printed\n',filename);
	pause
	print ('-dpdf', filename);
end
	
