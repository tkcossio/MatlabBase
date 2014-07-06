function g=AAM_NormalizeAppearance3D(gim,options)
if(options.normalizeg)
	% Normalize appearance data grey values
	color(2)=mean(gim);
	color(1)=std(gim)+eps;
	g=(gim-(color(2)))/(color(1));
else
	g=gim;
end
g(isnan(g))=0;
g(~isfinite(g))=0;