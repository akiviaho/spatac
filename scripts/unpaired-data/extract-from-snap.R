library(Matrix)
library(SnapATAC)
file.list = c("data/atac/cerebrum-study/CEMBA180426_8B.snap","data/atac/cerebrum-study/CEMBA190711_8J.snap", "data/atac/cerebrum-study/CEMBA180430_8B.snap", "data/atac/cerebrum-study/CEMBA190716_8E.snap", "data/atac/cerebrum-study/CEMBA190711_8E.snap", "data/atac/cerebrum-study/CEMBA190716_8J.snap")

sample.list = c('8B1','8J1','8B2','8E2','8E1','8J2')

data = createSnap(file=file.list, sample=sample.list);

data = addGmatToSnap(data)
data = addPmatToSnap(data)

pmat = data@pmat
gmat = data@gmat
meta = data@metaData

head(pmat)
head(gmat)
head(meta)

writeMM(pmat,file="aggregate-section-8-peak-matrix.mtx")
write.csv(data@barcode,'barcodes-aggregate-section-8.csv')
write.csv(data@peak$name,'peaks-aggregate-section-8-peak-matrix.csv')
writeMM(gmat,file="aggregate-section-8-gene-matrix.mtx")
write.csv(colnames(gmat),'genes-aggregate-section-8-gene-matrix.csv')

write.csv(meta,"section-8-metadata.csv")
