function fig2png(fname,fnameS)
%open *.fig & save into png & close

open(fname);
print(fnameS,'-dpng')
close