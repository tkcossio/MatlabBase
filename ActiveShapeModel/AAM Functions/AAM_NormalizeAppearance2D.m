function g=AAM_NormalizeAppearance2D(gim)
% Normalize appearance data grey values
color(2)=mean(gim);
color(1)=std(gim);
g=(gim-(color(2)))/(color(1));
