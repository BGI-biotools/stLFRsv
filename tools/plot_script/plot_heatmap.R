args=commandArgs(TRUE)
if( file.info(args[1])$size!=0){
	a=read.table(args[1])
	attach(a)
	r= as.numeric(a[,3])
	r.center = (r < quantile(r, as.numeric(args[8])))
	detach(a)
	a = a[r.center, ]
}else{
	 a= data.frame(matrix(ncol = 3, nrow = 0))
}
name=paste(args[2],"tiff",sep=".")
colnames(a)=c("x","y","Count")
cmin=as.numeric(args[4])
cmax=as.numeric(args[5])

tiff(name,width =2000,height =1500,units = "px",compression = "lzw",res=300)
library(ggplot2)
ggplot(a, aes(x=x, y=y, fill=Count))+geom_tile(width=as.numeric(args[9]))+
geom_hline(yintercept = as.numeric(args[6]),linetype="dotted")+
geom_hline(yintercept = as.numeric(args[7]),linetype="dotted")+
geom_vline(xintercept = as.numeric(args[6]),linetype="dotted")+
geom_vline(xintercept = as.numeric(args[7]),linetype="dotted")+
ylim(cmin,cmax)+xlim(cmin,cmax)+labs(title=args[3])+scale_fill_gradient(low='white', high='red')
dev.off()
