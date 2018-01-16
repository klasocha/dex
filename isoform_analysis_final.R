require(gplots)
require('dendextend')

##Definiujemy kilka pomocniczych funkcji
#nany zamiast błędów
my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#trzykolorowa paleta
colfunc <- colorRampPalette(c("blue", "white", "red"))

#spłaszczanie do thresholdów
thresh <- function(x, threshold = 2.5) {
  x[x < -threshold] <- -threshold
  x[x > threshold] <- threshold
  x
}

#odległość
corr_dist = function(x) as.dist((1-cor(t(x)))/2)


##Wczytywanie danych z plików. Wszystkie pliki z danymi muszą znajdować się wewnątrz data_dir (data_dir/CTRL-AST-3/isoforms.fpkm_tracking)
data_dir = "dex-project/cufflinks"
data_dir_regex = paste0(data_dir, "/*/isoforms.fpkm_tracking")
file_names = sort(Sys.glob(data_dir_regex))[-c(6,7,8,9,24,25,26,27,33)]
dataFiles <- lapply(file_names, read.table) #usuwamy wszystkie NAC i DEX_STR_2


#wczytujemy podstawowe informacje takie jak liczba izoform, ich id, nazwa genów i tworzymy macierz z wynikami fpkms
isoforms_number = nrow(dataFiles[[1]])-1
isoforms_ids = tail(dataFiles[[1]]$V1 , -1)
gene_names = tail(dataFiles[[1]]$V5 , -1)
isoforms_fpkms = matrix(nrow=length(dataFiles), ncol=isoforms_number)

#Z nazw plików wyciągamy treat, tissue i numer badania oraz wartości fpmks
treats = c()
tissues = c()
study_number = c()

for (i in 1:length(file_names)){
  isoforms_fpkms[i,] =  tail(as.numeric(dataFiles[[i]]$V10), -1)
  #wyciągnięcie parametrów z nazwy pliku
  treats = c(treats, strsplit(strsplit(file_names[i], "/")[[1]][3], "-")[[1]][1])
  tissues = c(tissues, strsplit(strsplit(file_names[i], "/")[[1]][3], "-")[[1]][2])
  study_number = c(study_number, strsplit(strsplit(file_names[i], "/")[[1]][3], "-")[[1]][3])
}

#Drobna zmiana formatu danych oraz utworzenie data framów opisujących poszczególne pliki i izoformy
fpkms = data.frame(t(isoforms_fpkms))
rownames(fpkms) <- isoforms_ids
colnames(fpkms) <- 1:length(file_names)
studies = data.frame(treat = treats, tissue = tissues, number = study_number)
isoforms = data.frame(id = isoforms_ids, gene = gene_names)

#sprzątanie niepotrzebnych obiektów
remove(dataFiles, file_names, isoforms_number, isoforms_ids, gene_names, isoforms_fpkms, treats, tissues, study_number, i)

#ANOVA dla izoform
Anova_P_vals = matrix(nrow=nrow(fpkms), ncol=3)

for (i in 1:nrow(fpkms)){
  result=aov(t(as.matrix(fpkms[i,]))~studies$treat*studies$tissue)
  Anova_P_vals[i,] = summary(result)[[1]][["Pr(>F)"]] [1:3]
}
remove(result, i)

#Poprawka FDR, wybranie istotnych, zapisanie wyników ANOVY na potrzeby raportu.
Anova_P_vals_fdr = matrix(p.adjust(Anova_P_vals, "fdr"), ncol =3)

significant_treat = which((Anova_P_vals_fdr[,1] < 0.05) %in% TRUE)
significant_treat_number = sum((Anova_P_vals_fdr[,1] < 0.05) %in% TRUE)
significant_tissue_number = sum((Anova_P_vals_fdr[,2] < 0.05) %in% TRUE)
significant_combination_number = sum((Anova_P_vals_fdr[,3] < 0.05) %in% TRUE)
save(significant_treat_number, significant_tissue_number, significant_combination_number,
     file = "dex-project/report/significant_number.RData")
remove(Anova_P_vals, Anova_P_vals_fdr, significant_treat_number, significant_tissue_number, significant_combination_number)

#T test dla znalezienia istotnych izoform
t_test_pvals = matrix(nrow=length(significant_treat), ncol=length(unique(studies$tissue)))
for (i in 1:length(unique(studies$tissue))){
  tissue = unique(studies$tissue)[i]
  controls = fpkms[significant_treat,studies$tissue == tissue & studies$treat == "CTRL"]
  dexes = fpkms[significant_treat,studies$tissue == tissue & studies$treat == "DEX"]
  
  for (j in 1:length(dexes[,1])){
    #pval = t.test(controls[j,], dexes[j,], var.equal = TRUE)$p.value
    pval = my.t.test.p.value(controls[j,], dexes[j,], var.equal = TRUE)
    t_test_pvals[j,i] = pval
    
  }
}

# Dla każdej tkanki z osobna wprowadzamy poprawkę FDR Wyniki zapisujemy na potrzeby raportu.
t_test_pvals_fdr = t_test_pvals
for (i in 1:length(unique(studies$tissue))){
  t_test_pvals_fdr[,i] = p.adjust(t_test_pvals[,i], "fdr")
  tissue = unique(studies$tissue)[i]
  print(as.character(tissue))
  print(sum(t_test_pvals_fdr[,i] < 0.05 , na.rm=TRUE))
}
save(t_test_pvals_fdr, file = "dex-project/report/t_test_pvals_fdr.RData")


remove(controls, dexes, i, j, pval, tissue)

#Wyplotujemy teraz te, które są istotne dla STR po FDR na poziomie 0.05.
STR_significant = as.matrix(fpkms[significant_treat[(t_test_pvals_fdr[,3] < 0.05) %in% TRUE],  ])

#W nazwie kolumny chcemy zawrzeć tkankę, treat i numer badania. Do id genu dodajemy id tkanki dla uniknięcia duplikatów.
colnames(STR_significant) <- sapply(colnames(STR_significant), function(x) paste0(studies[x,1], "-", studies[x,2], "-", studies[x,3]))
rownames(STR_significant) = paste0(isoforms$gene[significant_treat[(t_test_pvals_fdr[,3] < 0.05) %in% TRUE]], "_", isoforms$id[significant_treat[(t_test_pvals_fdr[,3] < 0.05) %in% TRUE]])

#Tworzymy nowy df na potrzeby plotowania. Będziemy go skalować itp.
STR_significant_plot = STR_significant

#Skalujemy na podstawie mediany i odchylenia kontroli, dla każdej tkanki osobno.
for (i in 1:length(unique(studies$tissue))){
  tissue = unique(studies$tissue)[i]
  med = apply(STR_significant[,studies$tissue == tissue & studies$treat == "CTRL"],1,median)
  sd = apply(STR_significant[,studies$tissue == tissue & studies$treat == "CTRL"],1,sd)
  sd[sd==0] <-1
  STR_significant_plot[,studies$tissue == tissue] <- (STR_significant[,studies$tissue == tissue] - med)/sd
  
}

#Całość skaluję i spłaszczam do thresholdów
gene_median =  apply(STR_significant,1,median)
STR_significant_plot <- t(apply(STR_significant_plot, 1, scale))
STR_significant_plot <- t(apply(STR_significant_plot, 1, thresh, threshold = 2.5))
colnames(STR_significant_plot) <- colnames(STR_significant)

#Wyrzucam geny dla ktorych mediana jest większa niż jeden, tworzę dendrogram i plotuję
hc = hclust(corr_dist(STR_significant_plot[gene_median > 1,]))
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 5 )

heatmap(STR_significant_plot[gene_median > 1,], col = colfunc((1000)), scale = "none", dist = function(x) as.dist((1-cor(t(x)))/2), Colv = NA, Rowv = dend)


#Następnie robię ro samo dla większego progu istotności
STR_significant01 = as.matrix(fpkms[significant_treat[(t_test_pvals_fdr[,3] < 0.1) %in% TRUE],  ])
colnames(STR_significant01) <- sapply(colnames(STR_significant01), function(x) paste0(studies[x,1], "-", studies[x,2], "-", studies[x,3]))
rownames(STR_significant01) = paste0(isoforms$gene[significant_treat[(t_test_pvals_fdr[,3] < 0.1) %in% TRUE]], "_", isoforms$id[significant_treat[(t_test_pvals_fdr[,3] < 0.1) %in% TRUE]])
STR_significant_plot01 = STR_significant01
for (i in 1:length(unique(studies$tissue))){
  tissue = unique(studies$tissue)[i]
  med = apply(STR_significant01[,studies$tissue == tissue & studies$treat == "CTRL"],1,median)
  sd = apply(STR_significant01[,studies$tissue == tissue & studies$treat == "CTRL"],1,sd)
  sd[sd==0] <-1
  STR_significant_plot01[,studies$tissue == tissue] <- (STR_significant01[,studies$tissue == tissue] - med)/sd
}
gene_median01 =  apply(STR_significant01,1,median)
STR_significant_plot01 <- t(apply(STR_significant_plot01, 1, scale))
STR_significant_plot01 <- t(apply(STR_significant_plot01, 1, thresh, threshold = 2.5))
colnames(STR_significant_plot01) <- colnames(STR_significant01)
hc01 = hclust(corr_dist(STR_significant_plot01[gene_median01 > 1,]))
dend01 <- as.dendrogram(hc01)
dend01 <- color_branches(dend01, k = 5 )
dend01_5cl <- color_branches(dend01, k = 5 )
dend01_3cl <- color_branches(dend01, k = 3 )
dend01_2cl <- color_branches(dend01, k = 2 )
htmap01 = heatmap(STR_significant_plot01[gene_median01 > 1,], col = colfunc((1000)), scale = "none", dist = function(x) as.dist((1-cor(t(x)))/2), Colv = NA, Rowv = dend01)

save(STR_significant_plot, gene_median, dend, STR_significant_plot01, gene_median01, dend01, colfunc, dend01_5cl, dend01_3cl, dend01_2cl, hc01, file = "dex-project/report/plot_data.RData")
remove(tissue, i, med, sd, t_test_pvals, dend, dend01, dend01_2cl, dend01_3cl, dend01_5cl, htmap01, 
       STR_significant_plot01, STR_significant_plot, STR_significant, STR_significant01)



#Peak distance analysis
significant_isoforms = isoforms$id[significant_treat[(t_test_pvals_fdr[,3] < 0.1) && (gene_median01 > 1) %in% TRUE]]
save(significant_isoforms, file = "dex-project/peak-distance/01significant-isoforms.RData")

remove(significant_treat, t_test_pvals_fdr)


### ENRICHR analysis
require("enrichR")




get_group_enrichr_data = function(x, clusters, tree) {
  I_Group = sapply(strsplit(names(which(cutree(tree, k = clusters) == x)), "_"), head, 1)
  enrichr(I_Group, databases = c("ENCODE_TF_ChIP-seq_2015", "NCI-Nature_2016", "DrugMatrix", "NCI-60_Cancer_Cell_Lines"))
}

get_group_enrichr_data_all = function(x, clusters, tree) {
  I_Group = sapply(strsplit(names(which(cutree(tree, k = clusters) == x)), "_"), head, 1)
  #NCI-Nature się wysypuje, inna liczba kolumn
  dbs <- setdiff(c(listEnrichrDbs()$libraryName) , c("NCI-Nature_2015"))
  enrichr(I_Group, databases = dbs)
}

enrichr_data <- lapply(1:5, get_group_enrichr_data, clusters = 5, tree = hc)
enrichr_data01 <- lapply(1:5, get_group_enrichr_data, clusters = 5, tree = hc01)

enrichr_data_5cl <- lapply(1:5, get_group_enrichr_data_all, clusters = 5, tree = hc01)
enrichr_data_3cl <- lapply(1:3, get_group_enrichr_data_all, clusters = 3, tree = hc01)
enrichr_data_2cl <- lapply(1:2, get_group_enrichr_data_all, clusters = 2, tree = hc01)

save(enrichr_data, enrichr_data01, enrichr_data_5cl, enrichr_data_3cl, enrichr_data_2cl, file = "dex-project/report/enrichr.RData")


#TMP
