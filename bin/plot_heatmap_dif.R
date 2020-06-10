args=commandArgs(TRUE)
a=read.table(args[1])
name=paste(args[2],"tiff",sep=".")
colnames(a)=c("x","y","Count")
xmin=as.numeric(args[4]);
xmax=as.numeric(args[5]);
attach(a)
r= as.numeric(a[,3])
r.center = (r < quantile(r, as.numeric(args[10])))
detach(a)
a = a[r.center, ]

ymin=as.numeric(args[6]);
ymax=as.numeric(args[7]);
tiff(name,width =2000,height =1500,units = "px",compression = "lzw",res=300)
library(ggplot2)
ggplot(a, aes(x=x, y=y, fill=Count))+geom_tile(width=as.numeric(args[11]))+
geom_hline(yintercept = as.numeric(args[9]),linetype="dotted")+
geom_vline(xintercept = as.numeric(args[8]),linetype="dotted")+
ylim(ymin,ymax)+xlim(xmin,xmax)+labs(title=args[3])+scale_fill_gradient(low='white', high='red')
dev.off()
