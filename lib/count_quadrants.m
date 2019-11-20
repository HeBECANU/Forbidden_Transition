function out=count_quadrants(txy_in)

% xy = ++,-+,+-,--
out=nan(4,1);
out(1)=sum(txy_in(:,2)>0 & txy_in(:,3)>0);
out(2)=sum(txy_in(:,2)<0 & txy_in(:,3)>0);
out(3)=sum(txy_in(:,2)>0 & txy_in(:,3)<0);
out(4)=sum(txy_in(:,2)<0 & txy_in(:,3)<0);

end