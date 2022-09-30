#setwd("./")
pdf('Image/tumor_normal.merge3-3_modify_now.pdf', height=18, width=22)
par(mfrow = c(3, 1), mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=3,yaxs='i')
file1 <- read.table('tumor_normal.bk.0.cpg.distr.num.win50.500')
plot(file1[,2],file1[,3], col='darkred', type='p', pch='',lwd='5',xlab='windows  (50bp)', ylab='ratio of break points', main='CPG',ylim=c(0,12),cex.main=3,cex.lab=3)

par(new=TRUE)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=3,yaxs='i')
file2 <- read.table('tumor_normal.bk.1.cpg.distr.num.win50.500')
plot(file2[,2],file2[,3], col='darkblue', type='p', pch='',lwd=5,xlab='windows  (50bp)', ylab='ratio of break points', main='CPG',ylim=c(0,12),yaxt='n',xaxt='n',cex.main=3,cex.lab=3)

par(new=TRUE)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=3,yaxs='i')
file3 <- read.table('hg19_CpGIsland.txt.len.pct2.500')
plot(file3[,2],file3[,4], col='gray', type='p', pch='',lwd=5,xlab='windows  (50bp)', ylab='ratio of break points', main='CPG',ylim=c(0,12),yaxt='n',xaxt='n',cex.main=3,cex.lab=3)


spl1 <- spline(file1[,2],file1[,3])
lines(spl1, col='darkred', lwd='6')

spl2 <- spline(file2[,2],file2[,3])
lines(spl2, col='darkblue', lwd='6')

spl3 <- spline(file3[,2],file3[,4])
lines(spl3, col='gray', lwd='6')

legend(450, 10, c('tumor','normal', 'expected'),lty=c(1,1),lwd=5,col=c('darkred','darkblue','gray'),cex=2)

#pdf('/ifs5/ST_IM/PMO/F12HPCECSZ0278/HBV-1000/HBV-825-20131022/Tumor_Normal/All_groups/norm10_analysis_5/gte10/tumor_normal/genome_element/tumor_normal.tfbs3-1.pdf',width=24,height=18)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=3,yaxs='i')
file1 <- read.table('tumor_normal.bk.0.tfbs.distr.num.win50.500')
plot(file1[,2],file1[,3], col='darkred', type='p', pch='',lwd='5',xlab='windows  (50bp)', ylab='ratio of break points', main='TFBS',ylim=c(0,12),cex.main=3,cex.lab=3)

par(new=TRUE)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=3,yaxs='i')
file2 <- read.table('tumor_normal.bk.1.tfbs.distr.num.win50.500')
plot(file2[,2],file2[,3], col='darkblue', type='p', pch='',lwd='5',xlab='windows  (50bp)', ylab='ratio of break points', main='TFBS',ylim=c(0,12),yaxt='n',xaxt='n',cex.main=3,cex.lab=3)

par(new=TRUE)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=3,yaxs='i')
file3 <- read.table('/public/home/chshen/auto_Image_Pro/tfbsConsSites.txt.len500')
plot(file3[,2],file3[,3], col='gray', type='p', pch='',lwd=5,xlab='windows  (50bp)', ylab='ratio of break points', main='TFBS',ylim=c(0,12),yaxt='n',xaxt='n',cex.main=3,cex.lab=3)

spl1 <- spline(file1[,2],file1[,3])
lines(spl1, col='darkred', lwd='6')

spl2 <- spline(file2[,2],file2[,3])
lines(spl2, col='darkblue', lwd='6')

spl3 <- spline(file3[,2],file3[,3])
lines(spl3, col='gray', lwd='6')

legend(450, 10, c('tumor','normal', 'expected'),lty=c(1,1),lwd=5,col=c('darkred','darkblue','gray'),cex=2)
#pdf('/ifs5/ST_IM/PMO/F12HPCECSZ0278/HBV-1000/HBV-825-20131022/Tumor_Normal/All_groups/norm10_analysis_5/gte10/tumor_normal/genome_element/tumor_normal.tss3-1.pdf',width=20,height=16)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=2,yaxs='i')
file1 <- read.table('tumor_normal.bk.0.tss.distr.num.win50.500')
plot(file1[,2],file1[,3], col='darkred', type='p', pch='',lwd='5',xlab='windows  (50bp)', ylab='ratio of break points', main='TSS',ylim=c(0,12),cex.main=3,cex.lab=3, cex.axis=3)

par(new=TRUE)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=2,yaxs='i')
file2 <- read.table('tumor_normal.bk.1.tss.distr.num.win50.500')
plot(file2[,2],file2[,3], col='darkblue', type='p', pch='',lwd='5',xlab='windows  (50bp)', ylab='ratio of break points', main='TSS',ylim=c(0,12),yaxt='n',xaxt='n',cex.main=3,cex.lab=3, cex.axis=3)

par(new=TRUE)
par(mar=c(6,10,6,6),mgp=c(5,2,0),bty = 'l',cex.axis=2,yaxs='i')
file3 <- read.table('/public/home/chshen/auto_Image_Pro//switchDbTss.txt.len.pct.500.ts')
plot(file3[,2],file3[,3], col='gray', type='p', pch='',lwd=5,xlab='windows  (50bp)', ylab='ratio of break points', main='TSS',ylim=c(0,12),yaxt='n',xaxt='n',cex.main=3,cex.lab=3, cex.axis=3)

spl1 <- spline(file1[,2],file1[,3])
lines(spl1, col='darkred', lwd='6')

spl2 <- spline(file2[,2],file2[,3])
lines(spl2, col='darkblue', lwd='6')

spl3 <- spline(file3[,2],file3[,3])
lines(spl3, col='gray', lwd='6')

par(new=TRUE)
legend(450, 10, c('tumor','normal', 'expected'),lty=c(1,1),lwd=5,col=c('darkred','darkblue','gray'),cex=2)
