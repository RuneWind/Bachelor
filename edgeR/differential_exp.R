library(Rsubread)
library(edgeR)
library(ggplot2)
library(dplyr)
library(viridis)
library(gplots)

setwd("./STAR")

count_matrix = read.delim(file="featureCounts.txt", header = TRUE, sep = "\t", row.names = "Geneid")
count_matrix_trimmed <- count_matrix[,6:150] # Remove columns not containing read counts
View(count_matrix_trimmed)


# Do analysis on full dataset
group = factor(c('Aearl', 'Clfin', 'Mrida', 'Mrida', 'Aearl', 'Rbani', 'Rbani', 'Sster', 'Aearl', 'Banna', 'Volin', 'Banna', 'Banna', 'Ctain', 'Aearl', 'Aearl', 'Clfin', 'Llanc', 'Aearl', 'Aalon', 'Rbani', 'Rling', 'Ctain', 'Ctain', 'Kdike', 'Kdike', 'Borek', 'Borek', 'Ancor', 'Aoost', 'Ccyma', 'Volin', 'Banca', 'Rbani', 'Aoost', 'Ctain', 'Iiona', 'Aoost', 'Clfin', 'Aalon', 'Banna', 'Aearl', 'Banca', 'Sster', 'Aearl', 'Aalon', 'Aalon', 'Aalon', 'Rbani', 'Ctain', 'Aearl', 'Ancor', 'Mrida', 'Clfin', 'Ancor', 'Volin', 'Rling', 'Rling', 'Aearl', 'Ctain', 'Banna', 'Clfin', 'Banna', 'Ancor', 'Aalon', 'Ccyma', 'Clfin', 'Aalon', 'Ccyma', 'Ancor', 'Iiona', 'Volin', 'Aearl', 'Iiona', 'Ancor', 'Borek', 'Iiona', 'Clfin', 'Volin', 'Sster', 'Rbani', 'Aaran', 'Aearl', 'Iiona', 'Kdike', 'Clfin', 'Iiona', 'Llanc', 'Kdike', 'Mrida', 'Aearl', 'Borek', 'Banca', 'Volin', 'Ccyma', 'Banca', 'Ancor', 'Kdike', 'Volin', 'Banna', 'Llanc', 'Sster', 'Llanc', 'Banca', 'Clfin', 'Ccyma', 'Sster', 'Clfin', 'Ancor', 'Aearl', 'Sster', 'Llanc', 'Ctain', 'Banna', 'Sster', 'Aoost', 'Aaran', 'Aearl', 'Llanc', 'Ancor', 'Aoost', 'Ccyma', 'Clfin', 'Aaran', 'Iiona', 'Rbani', 'Llanc', 'Ctain', 'Mrida', 'Ctain', 'Rbani', 'Banna', 'Aalon', 'Aaran', 'Ctain', 'Sster', 'Aoost', 'Rling', 'Aearl', 'Aalon', 'Rbani', 'Aoost', 'Aearl', 'Iiona', 'Iiona'))

# Do analysis without NA's
na = c("Aearl_073","Aearl_098","Aearl_121","Aearl_130", "Aearl_157")
count_matrix_trimmed = select(count_matrix_trimmed, -c(na))
group = factor(c('Aearl', 'Clfin', 'Mrida', 'Mrida', 'Aearl', 'Rbani', 'Rbani', 'Sster', 'Aearl', 'Banna', 'Volin', 'Banna', 'Banna', 'Ctain', 'Aearl', 'Aearl', 'Clfin', 'Llanc', 'Aearl', 'Aalon', 'Rbani', 'Rling', 'Ctain', 'Ctain', 'Kdike', 'Kdike', 'Borek', 'Borek', 'Ancor', 'Aoost', 'Ccyma', 'Volin', 'Banca', 'Rbani', 'Aoost', 'Ctain', 'Iiona', 'Aoost', 'Clfin', 'Aalon', 'Banna', 'Aearl', 'Banca', 'Sster', 'Aearl', 'Aalon', 'Aalon', 'Aalon', 'Rbani', 'Ctain', 'Aearl', 'Ancor', 'Mrida', 'Clfin', 'Ancor', 'Volin', 'Rling', 'Rling', 'Aearl', 'Ctain', 'Banna', 'Clfin', 'Banna', 'Ancor', 'Aalon', 'Ccyma', 'Clfin', 'Aalon', 'Ccyma', 'Ancor', 'Iiona', 'Volin', 'Iiona', 'Ancor', 'Borek', 'Iiona', 'Clfin', 'Volin', 'Sster', 'Rbani', 'Aaran', 'Aearl', 'Iiona', 'Kdike', 'Clfin', 'Iiona', 'Llanc', 'Kdike', 'Mrida', 'Borek', 'Banca', 'Volin', 'Ccyma', 'Banca', 'Ancor', 'Kdike', 'Volin', 'Banna', 'Llanc', 'Sster', 'Llanc', 'Banca', 'Clfin', 'Ccyma', 'Sster', 'Clfin', 'Ancor', 'Sster', 'Llanc', 'Ctain', 'Banna', 'Sster', 'Aoost', 'Aaran', 'Llanc', 'Ancor', 'Aoost', 'Ccyma', 'Clfin', 'Aaran', 'Iiona', 'Rbani', 'Llanc', 'Ctain', 'Mrida', 'Ctain', 'Rbani', 'Banna', 'Aalon', 'Aaran', 'Ctain', 'Sster', 'Aoost', 'Rling', 'Aalon', 'Rbani', 'Aoost', 'Aearl', 'Iiona', 'Iiona'))
ncol(count_matrix_trimmed)

DGE = DGEList(counts=count_matrix_trimmed, group=group)
DGE = calcNormFactors(DGE)
keep = rowSums(cpm(DGE$counts) > 1) > 70
View(keep)
DGE = DGE[keep, , keep.lib.sizes=FALSE]
nrow(DGE$counts)


# Samples for control PCA
controls = select(count_matrix_trimmed, c("Aearl_005", "Aearl_016", "Aearl_019", "Clfin_002", "Clfin_017", "Ctain_014", "Ctain_024", "Aearl_001", "Aearl_015", "Aearl_051"))
colnames(controls) = c("Aearl_2_005", "Aearl_2_016", "Aearl_2_019", "Clfin_2_002", "Clfin_2_017", "Ctain_2_014", "Ctain_2_024", "Aearl_1_001", "Aearl_2_015", "Aearl_3_051")
genotyper = c("Aearl_0529", "Aearl_0529", "Aearl_0529", "Clfin_0320", "Clfin_0320", "Ctain_0636", "Ctain_0636", "Aearl_0749", "Aearl_0749", "Aearl_0749")
View(controls)

# PCA plot
pc_df = prcomp(controls)
percentage = round(pc_df$sdev / sum(pc_df$sdev) * 100, 2)
xlabel = paste("Leading logFC dim 1 (", as.character(percentage[1]), "%)", sep="")
ylabel = paste("Leading logFC dim 2 (", as.character(percentage[2]), "%)", sep="")

plot_values = plotMDS(controls, col=genotyper, pch = 19, plot = FALSE)

plot_data = matrix(plot_values$x)
plot_data = cbind(plot_data, plot_values$y)
plot_data = cbind(plot_data, genotyper)
colnames(plot_data) = c("x", "y", "Sorter")
for (i in 1:length(plot_values$x)){
  plot_data[i,3] = substr(plot_data[i, 3], start = 1, stop = 5)
}

plot_data = as.data.frame(plot_data)
factor.to.numeric <- function(x) as.numeric(levels(x)[x])
plot_data$x = factor.to.numeric(plot_data$x)
plot_data$y = factor.to.numeric(plot_data$y)

hulls <- data.frame()
for(Sort in unique(plot_data$Sorter)){
  Subdf <- plot_data %>% filter(Sorter==Sort)
  hull <- Subdf[chull(Subdf$x, Subdf$y),]
  hulls <- rbind(hulls, hull)
}

colorvector1 = c("#a76d00","#4168ff","#59e069","#420097","#acd44f",
                 "#ff71f1","#00931d","#de006b","#47dada","#fe0e16",
                 "#002552","#9a9900","#70006e","#2d6800","#980016",
                 "#003003","#f6b89c","#4e0726","#ab9173","#504200")

figure = ggplot(plot_data, aes(x=x, y=y, color=Sorter, label=row.names(plot_data)))
figure = figure + geom_point()
figure = figure + theme_bw()
figure = figure + geom_polygon(data = hulls, aes(x=x, y=y, fill=Sorter), alpha = 0.1)
figure = figure + scale_color_manual(values=colorvector1)
figure = figure + scale_fill_manual(values=colorvector1)
figure = figure + labs(x=xlabel, y=ylabel)
figure = figure + geom_text(vjust=1.5)
print(figure)



# Continue analysis ###################
design = model.matrix(~0+group, data=DGE$samples)
colnames(design) = levels(DGE$samples$group)

dispersion = estimateDisp(DGE, design) # Estimates common and tagwise dispersion
plotBCV(dispersion)

# Testing for differential expression
fit = glmQLFit(dispersion, design)
plotQLDisp(fit)

########################### Contrast all at once
con_vector = c()
vector_idx = 0
for (i in 1:19){
  for (j in (i):19){
    vector_idx = vector_idx + 1
    con_name = paste(levels(group)[i], levels(group)[j+1], sep = "-")
    con_vector[vector_idx] = con_name
  }
}
con_vector
contrasts = makeContrasts(contrasts=con_vector, levels=design)
contrasts

qlf = glmQLFTest(fit, contrast=contrasts)
toptags = topTags(qlf)
summary(decideTests(qlf))


######################### Contrast pairwise
levels(group)
de_genes = data.frame()
for (i in 1:(length(levels(group)))){
  print(paste("###############", levels(group)[i], sep=" "))
  for (j in 1:(length(levels(group)))){
    print(colnames(de_genes))
    if(i == j){
      print(paste(paste("Skipping contrast of", levels(group)[i], sep=" "), levels(group)[j], sep = " and "))
    }
    else{
      print(paste(levels(group)[i], levels(group)[j], sep=" - "))
      pw_contrast = makeContrasts(contrasts=c(paste(levels(group)[i], levels(group)[j], sep = "-")), levels=design)
      pw_qlf = glmQLFTest(fit, contrast=pw_contrast)
      con_sum = summary(decideTests(pw_qlf))
      dt = as.data.frame(decideTests(pw_qlf))
      if(i == 1 & j == 2){
        con_matrix = matrix(con_sum, nrow = 3, dimnames = list(c("Down", "Not sig", "Up"), colnames(con_sum)))
        de_genes = dt
      }
      else{
        con_matrix = cbind(con_matrix, con_sum)
        de_genes = cbind(de_genes, dt)
        if (j < i){
          colnames(de_genes)[ncol(de_genes)] = paste(strsplit(colnames(decideTests(pw_qlf)), " ")[[1]][[2]], strsplit(colnames(decideTests(pw_qlf)), " ")[[1]][[1]])
        }
      }
    }
  } 
}
View(de_genes)
no_null = de_genes[apply(de_genes, 1, function(x) !all(x==0)),]
View(no_null[ , order(names(no_null))])
write.csv(no_null, file = "no_null.csv")


# Heatmap of DE-genes
DE.heatmap <- as.matrix(tail(hm_data, n=250))
DE.heatmap = DE.heatmap + 1
View(DE.heatmap)
DE.heatmap = log10(DE.heatmap)
labCol = substr(x=colnames(DE.heatmap), start=1, stop=5)
col.scale <- viridis(14)
heat = heatmap.2(DE.heatmap, col=col.scale, Rowv=FALSE, Colv=FALSE,  dendrogram="none", trace="none", scale="row", labCol=labCol, key=TRUE, offsetRow=0, adjCol=c(0, 0.5), offsetCol=1.4, linecol="black", add.expr=abline(v=c(9.5, 13.5, 30.5, 39.5, 46.5, 51.5, 60.5, 64.5, 70.5, 81.5, 91.5, 100.5, 105.5, 112.5, 117.5, 126.5, 130.5, 138.5, 145.5)))
dev.off()


log10(scale(as.vector(DE.heatmap[(nrow(DE.heatmap)) , ])))

# PCA-plot of heatmap data
hm_data = read.csv("./STAR/countMatrix_DE.csv", sep=";", row.names="X")
hm_pca = prcomp(hm_data)
percentage = round(hm_pca$sdev / sum(hm_pca$sdev) * 100, 2)
xlabel = paste("Leading logFC dim 1 (", as.character(percentage[1]), "%)", sep="")
ylabel = paste("Leading logFC dim 2 (", as.character(percentage[2]), "%)", sep="")

plot_values = plotMDS(hm_data, plot = FALSE)
plot_data = matrix(plot_values$x)
plot_data = cbind(plot_data, plot_values$y)
plot_data = cbind(plot_data, row.names(plot_data))
colnames(plot_data) = c("x", "y", "Sorter")
for (i in 1:length(plot_values$x)){
  plot_data[i,3] = substr(plot_data[i, 3], start = 1, stop = 5)
}

plot_data = as.data.frame(plot_data)
factor.to.numeric <- function(x) as.numeric(levels(x)[x])
plot_data$x = factor.to.numeric(plot_data$x)
plot_data$y = factor.to.numeric(plot_data$y)

hulls <- data.frame()
for(Sort in unique(plot_data$Sorter)){
  Subdf <- plot_data %>% filter(Sorter==Sort)
  hull <- Subdf[chull(Subdf$x, Subdf$y),]
  hulls <- rbind(hulls, hull)
}

colorvector1 = c("#a76d00","#4168ff","#59e069","#420097","#acd44f",
                 "#ff71f1","#00931d","#de006b","#47dada","#fe0e16",
                 "#002552","#9a9900","#70006e","#2d6800","#980016",
                 "#003003","#f6b89c","#4e0726","#ab9173","#504200")

figure = ggplot(plot_data, aes(x=x, y=y, color=Sorter))
figure = figure + geom_point()
figure = figure + theme_bw()
figure = figure + geom_polygon(data = hulls, aes(x=x, y=y, fill=Sorter), alpha = 0.1)
figure = figure + scale_color_manual(values=colorvector1)
figure = figure + scale_fill_manual(values=colorvector1)
figure = figure + labs(x=xlabel, y=ylabel)
print(figure)
