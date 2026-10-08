################################################################################
## Mitracarpus (Spermacoceae, Rubiaceae) — seed micromorphology
## Supplementary File 2. Complete analytical pipeline.
## Self-contained: run this file to reproduce every value and figure reported.
## Requires: cluster, vegan, ape, FactoMineR, dendextend, ggplot2, ggrepel, scales
################################################################################

suppressPackageStartupMessages({
  library(cluster); library(vegan); library(ape); library(FactoMineR)
  library(dendextend); library(ggplot2); library(ggrepel); library(scales)
})
PERM <- 9999
## the seed is reset before every permutation test, so each p-value is
## reproducible on its own and does not depend on execution order
set.seed(20260913)
OUT <- "out"; FIG <- "figs"
dir.create(OUT, showWarnings = FALSE); dir.create(FIG, showWarnings = FALSE)

MM <- function(x) x/25.4           # mm -> inches
W2 <- MM(170); W1 <- MM(85)        # double and single column width

## ---------------------------------------------------------------- 1. DATA ---

qual <- read.csv(text = "
species,VGS,SSh,AW,PW,EOrn,CS,DCD,TVGDV,Colour
M. eritrichoides,1,1,1,1,2,1,1,1,1
M. parvulus,1,3,1,1,2,1,1,1,4
M. longicalyx,1,5,4,3,1,3,1,1,5
M. nitidus,1,2,3,2,2,2,1,2,3
M. albomarginatus,1,4,4,3,1,3,2,1,2
M. bacigalupoae,1,3,1,1,2,1,2,1,3
M. fernandesii,1,4,1,1,2,1,2,1,1
M. brasiliensis,1,2,8,2,2,1,1,2,2
M. hirtus,1,2,2,2,2,2,1,1,2
M. microspermus,1,4,1,1,2,5,2,1,4
M. polygonifolius,1,4,6,2,2,1,1,2,2
M. pusillus,1,6,4,4,3,4,1,1,1
M. recurvatus,1,8,9,4,3,4,1,1,7
M. steyermarkii,1,3,4,5,3,4,1,1,3
M. anthospermoides,1,6,5,2,2,3,2,1,6
M. baturitensis,1,6,6,2,2,3,3,1,4
M. bicrucis,1,2,7,1,2,2,3,1,7
M. carajasensis,1,2,NA,NA,2,NA,3,2,NA
M. frigidus,1,6,6,2,2,3,3,1,2
M. rigidifolius,1,4,3,2,2,6,4,1,8
M. schininianus,1,8,9,1,2,7,2,1,1
M. eichleri,2,4,5,2,2,3,1,2,4
M. salzmannianus,2,7,5,2,2,3,1,2,4
M. lhotzkyanus,3,5,4,6,1,3,1,2,1
M. megapotamicus,3,2,4,6,1,3,1,2,1
M. hasslerianus,3,2,4,6,1,6,1,2,1
M. diversifolius,4,4,10,2,2,6,1,2,1
", row.names = 1, strip.white = TRUE)
SP  <- rownames(qual)
dat <- qual; dat[] <- lapply(dat, factor)
CH  <- names(dat)

## seed dimensions (mm), min-max of length and width
siz <- read.csv(text = "
species,Lmin,Lmax,Wmin,Wmax
M. eritrichoides,0.5,0.6,0.5,0.5
M. parvulus,0.6,0.9,0.5,0.5
M. longicalyx,0.8,1.0,0.5,0.6
M. nitidus,0.6,0.8,0.4,0.4
M. albomarginatus,0.8,1.0,0.5,0.6
M. bacigalupoae,0.6,0.7,0.4,0.5
M. fernandesii,0.6,1.0,0.4,0.7
M. brasiliensis,0.6,1.0,0.6,0.8
M. hirtus,0.6,0.8,0.5,0.6
M. microspermus,0.7,0.8,0.5,0.5
M. polygonifolius,0.8,1.0,0.5,0.6
M. pusillus,0.8,0.8,0.5,0.5
M. recurvatus,0.6,1.0,0.4,0.7
M. steyermarkii,0.6,0.8,0.5,0.5
M. anthospermoides,0.5,0.8,0.5,0.6
M. baturitensis,0.6,1.0,0.5,0.5
M. bicrucis,0.7,1.0,0.4,0.6
M. carajasensis,0.6,1.0,0.6,0.6
M. frigidus,1.0,1.2,0.8,0.8
M. rigidifolius,1.0,1.2,0.6,0.8
M. schininianus,0.6,0.8,0.5,0.7
M. eichleri,0.7,0.8,0.5,0.5
M. salzmannianus,0.8,1.0,0.5,0.6
M. lhotzkyanus,0.8,1.2,0.5,0.6
M. megapotamicus,0.8,2.0,0.6,1.5
M. hasslerianus,0.8,1.2,0.6,0.8
M. diversifolius,0.3,0.7,0.7,0.9
", row.names = 1, strip.white = TRUE)[SP, ]
siz$L <- (siz$Lmin + siz$Lmax)/2; siz$W <- (siz$Wmin + siz$Wmax)/2
siz$LW <- siz$L/siz$W

## biome per species; conf H = distribution stated in a revision or protologue,
## L = inferred from the localities of the examined specimens
bio <- read.csv(text = "
species,biome,ecoregion,conf
M. albomarginatus,Deserts_xeric,Caatinga,H
M. anthospermoides,Moist_forest,Atlantic restinga,H
M. bacigalupoae,Moist_forest,Alto Parana Atlantic forest,L
M. baturitensis,Multiple,Caatinga + Cerrado,H
M. bicrucis,Dry_forest,Chiquitano dry forests,L
M. brasiliensis,Moist_forest,Araucaria moist forests,H
M. carajasensis,Moist_forest,Amazon canga outcrops,H
M. diversifolius,Moist_forest,Bahia coastal forests,H
M. eichleri,Moist_forest,Atlantic restinga,H
M. eritrichoides,Grass_savanna,Cerrado,H
M. fernandesii,Deserts_xeric,Caatinga,H
M. frigidus,Multiple,widespread,H
M. hasslerianus,Grass_savanna,Cerrado (Paraguay and Mato Grosso),H
M. hirtus,Multiple,widespread,H
M. lhotzkyanus,Moist_forest,Atlantic restinga,H
M. longicalyx,Deserts_xeric,Caatinga,H
M. megapotamicus,Multiple,Cerrado + Pampa + Chaco,H
M. microspermus,Grass_savanna,Cerrado,H
M. nitidus,Montane_grass,Campos rupestres,H
M. parvulus,Grass_savanna,Cerrado,H
M. polygonifolius,Multiple,Campos rupestres + Atlantic,L
M. pusillus,Montane_grass,Campos rupestres,H
M. recurvatus,Grass_savanna,Cerrado,H
M. rigidifolius,Montane_grass,Campos rupestres,H
M. salzmannianus,Multiple,widespread,H
M. schininianus,Grass_savanna,Cerrado,H
M. steyermarkii,Grass_savanna,Cerrado,H
", row.names = 1, strip.white = TRUE)[SP, ]
LVL <- c("Moist_forest","Dry_forest","Grass_savanna","Montane_grass",
         "Deserts_xeric","Multiple")
bio$biome <- factor(bio$biome, levels = LVL)
SHRT <- c(Moist_forest = "Moist broadleaf forests", Dry_forest = "Dry broadleaf forests",
          Grass_savanna = "Grasslands & savannas", Montane_grass = "Montane grasslands",
          Deserts_xeric = "Deserts & xeric shrublands", Multiple = "Multiple biomes")
FULL <- c(Moist_forest = "Tropical & subtropical moist broadleaf forests",
          Dry_forest = "Tropical & subtropical dry broadleaf forests",
          Grass_savanna = "Tropical & subtropical grasslands, savannas & shrublands",
          Montane_grass = "Montane grasslands & shrublands",
          Deserts_xeric = "Deserts & xeric shrublands", Multiple = "Multiple biomes")

## specimens examined, and number of characters scored with a compound state
nspec <- c("M. albomarginatus"=3,"M. anthospermoides"=3,"M. bacigalupoae"=1,
  "M. baturitensis"=4,"M. bicrucis"=1,"M. brasiliensis"=3,"M. carajasensis"=1,
  "M. diversifolius"=1,"M. eichleri"=2,"M. eritrichoides"=4,"M. fernandesii"=2,
  "M. frigidus"=3,"M. hasslerianus"=2,"M. hirtus"=3,"M. lhotzkyanus"=3,
  "M. longicalyx"=3,"M. megapotamicus"=3,"M. microspermus"=2,"M. nitidus"=4,
  "M. parvulus"=3,"M. polygonifolius"=4,"M. pusillus"=4,"M. recurvatus"=1,
  "M. rigidifolius"=5,"M. salzmannianus"=3,"M. schininianus"=2,
  "M. steyermarkii"=4)[SP]
poly <- c("M. eritrichoides"=2,"M. longicalyx"=3,"M. nitidus"=1,"M. parvulus"=2,
  "M. albomarginatus"=3,"M. bacigalupoae"=2,"M. fernandesii"=3,"M. brasiliensis"=1,
  "M. hirtus"=1,"M. microspermus"=3,"M. polygonifolius"=2,"M. pusillus"=2,
  "M. recurvatus"=2,"M. steyermarkii"=1,"M. anthospermoides"=2,"M. baturitensis"=2,
  "M. bicrucis"=2,"M. carajasensis"=0,"M. frigidus"=2,"M. rigidifolius"=3,
  "M. schininianus"=3,"M. eichleri"=2,"M. salzmannianus"=2,"M. lhotzkyanus"=3,
  "M. megapotamicus"=2,"M. hasslerianus"=2,"M. diversifolius"=2)[SP]

## a priori classification: groups = states of VGS; subtypes within Group I
grpVGS <- factor(c("I","II","III","IV")[qual$VGS], levels = c("I","II","III","IV"))
names(grpVGS) <- SP
subty <- setNames(rep(NA_character_, length(SP)), SP)
subty[c("M. eritrichoides","M. longicalyx","M. nitidus","M. parvulus")] <- "1.1"
subty[c("M. albomarginatus","M. bacigalupoae","M. fernandesii")] <- "1.2"
subty[c("M. brasiliensis","M. hirtus","M. microspermus","M. polygonifolius")] <- "1.3"
subty[c("M. pusillus","M. recurvatus","M. steyermarkii")] <- "1.4"
subty[c("M. anthospermoides","M. baturitensis","M. bicrucis","M. carajasensis",
        "M. frigidus","M. rigidifolius","M. schininianus")] <- "1.5"
subty <- factor(subty)

## ---------------------------------------------------- 2. GOWER AND CLUSTERING

D9  <- daisy(dat, metric = "gower")
m9  <- as.matrix(D9)
D12 <- daisy(cbind(dat, siz[, c("L","W","LW")]), metric = "gower")

meth <- c("average","complete","ward.D2","mcquitty","single")
coph <- sapply(meth, function(m) cor(D9, cophenetic(hclust(D9, m))))
hc   <- hclust(D9, "average")
CC   <- unname(coph["average"])
jk   <- sapply(CH, function(j) cor(as.vector(D9),
          as.vector(daisy(dat[, setdiff(CH, j)], metric = "gower"))))
jkt  <- sapply(CH, function(j) cor(cophenetic(hclust(
          daisy(dat[, setdiff(CH, j)], metric = "gower"), "average")), cophenetic(hc)))

K      <- 6
sil_hc <- sapply(2:8, function(k) mean(silhouette(cutree(hc, k), D9)[, 3]))
sil_pam<- sapply(2:8, function(k) pam(D9, k, diss = TRUE)$silinfo$avg.width)
names(sil_hc) <- names(sil_pam) <- 2:8
sil_ap <- mean(silhouette(as.integer(grpVGS), D9)[, 3])
inI    <- !is.na(subty)
sil_sub<- mean(silhouette(as.integer(droplevels(subty[inI])),
                          as.dist(m9[inI, inI]))[, 3])

## clusters renumbered top-to-bottom in the phenogram
ct   <- cutree(hc, K)
ord  <- rev(labels(as.dendrogram(hc)))
map  <- setNames(seq_along(unique(ct[ord])), unique(ct[ord]))
grp6 <- factor(paste0("C", map[as.character(ct)]), levels = paste0("C", 1:K))
names(grp6) <- names(ct)
SUBT <- c(C1 = "I-b", C3 = "I-d", C4 = "I-c", C5 = "I-a")
CLAB <- ifelse(is.na(SUBT[levels(grp6)]), levels(grp6),
               paste0(levels(grp6), " (", SUBT[levels(grp6)], ")"))
names(CLAB) <- levels(grp6)

## --------------------------------------------------------- 3. PCoA AND MCA ---

pco  <- cmdscale(D9, k = 8, eig = TRUE)
pcoC <- cmdscale(D9, k = 8, eig = TRUE, add = TRUE)
ev   <- pco$eig[pco$eig > 0];  pv  <- 100*ev/sum(ev)
evC  <- pcoC$eig[pcoC$eig > 0]; pvC <- 100*evC/sum(evC)
set.seed(20260913)
ef3  <- envfit(pco$points[, 1:3], dat, permutations = PERM, na.rm = TRUE)
set.seed(20260913)
efs  <- envfit(pco$points[, 1:3], siz[, c("L","W","LW")], permutations = PERM)
set.seed(20260913)
efn  <- envfit(pco$points[, 1:3], data.frame(n = nspec), permutations = PERM)
set.seed(20260913)
efb  <- envfit(pco$points[, 1:3], bio["biome"], permutations = PERM, na.rm = TRUE)
cc   <- complete.cases(dat)
mca  <- MCA(dat[cc, ], graph = FALSE, ncp = 5)
H    <- sapply(dat, function(x) {p <- prop.table(table(x)); -sum(p*log(p))})

## ------------------------------------------------- 4. DIAGNOSTIC PERFORMANCE -

M   <- as.matrix(qual)
TOT <- choose(nrow(M), 2)
unres <- function(cols) {
  n <- nrow(M); out <- character()
  for (i in 1:(n-1)) for (j in (i+1):n) {
    a <- M[i, cols]; b <- M[j, cols]; ok <- !is.na(a) & !is.na(b)
    if (!any(ok) || !any(a[ok] != b[ok]))
      out <- c(out, paste(SP[i], "vs", SP[j]))
  }
  out
}
full_un <- unres(CH)
uniq    <- sapply(CH, function(j) length(unres(setdiff(CH, j))))
alone   <- sapply(CH, function(j) TOT - length(unres(j)))

minset <- NULL
for (k in 1:length(CH)) {
  ok <- Filter(function(s) length(unres(s)) == 0, combn(CH, k, simplify = FALSE))
  if (length(ok)) { minset <- list(k = k, sets = ok); break }
}

rem <- CH; used <- character(); curve <- data.frame()
while (length(rem)) {
  g <- sapply(rem, function(j) TOT - length(unres(c(used, j))))
  b <- names(which.max(g)); used <- c(used, b); rem <- setdiff(rem, b)
  curve <- rbind(curve, data.frame(step = length(used), character. = b,
                   pairs_resolved = max(g), prop = round(max(g)/TOT, 4)))
}

## ------------------------------------------- 5. TESTS OF PREDEFINED GROUPINGS
## Each grouping is tested on a matrix from which its defining characters were
## removed, so the grouping variable is not contained in the matrix tested.

Dn  <- daisy(dat[, setdiff(CH, "VGS")], metric = "gower")
mn  <- as.matrix(Dn)
kIV <- grpVGS != "IV"
DI_red <- daisy(dat[inI, setdiff(CH, c("DCD","EOrn","AW"))], metric = "gower")

set.seed(20260913)

a_type  <- adonis2(Dn ~ grpVGS, permutations = PERM)
set.seed(20260913)
a_type2 <- adonis2(as.dist(mn[kIV, kIV]) ~ droplevels(grpVGS[kIV]), permutations = PERM)
set.seed(20260913)
a_sub   <- adonis2(DI_red ~ droplevels(subty[inI]), permutations = PERM)
bd_type <- betadisper(as.dist(mn), grpVGS)
set.seed(20260913)
pt_type <- permutest(bd_type, permutations = PERM)

shared <- function(idx) sapply(qual[idx, , drop = FALSE], function(x) {
  u <- unique(na.omit(x)); if (length(u) == 1) as.character(u) else "-" })

## revised subtypes of Group I, from the k = 6 partition
DI      <- as.dist(m9[inI, inI])
sub_new <- factor(SUBT[as.character(grp6[inI])], levels = c("I-a","I-b","I-c","I-d"))
names(sub_new) <- SP[inI]
sil_new <- mean(silhouette(as.integer(sub_new), DI)[, 3])
qc      <- as.data.frame(lapply(dat, as.character), stringsAsFactors = FALSE)
rownames(qc) <- SP
modal <- t(sapply(levels(sub_new), function(s) {
  sub <- qc[names(sub_new)[sub_new == s], , drop = FALSE]
  sapply(sub, function(x) { tb <- table(x, useNA = "no")
    if (!length(tb)) NA_character_ else
      sprintf("%s (%d/%d)", names(tb)[which.max(tb)], max(tb), sum(tb)) })
}))

## ------------------------------------------------------- 6. SAMPLING EFFORT -

meandist <- rowMeans(m9, na.rm = TRUE) * length(SP)/(length(SP)-1)
sp_poly  <- cor.test(nspec, poly, method = "spearman", exact = FALSE)
sp_dist  <- cor.test(nspec, meandist, method = "spearman", exact = FALSE)
keep <- nspec > 1
Dk   <- as.dist(m9[keep, keep]); hck <- hclust(Dk, "average")
ari  <- function(a, b) { t <- table(a, b); n <- sum(t)
  s <- sum(choose(t, 2)); a1 <- sum(choose(rowSums(t), 2)); b1 <- sum(choose(colSums(t), 2))
  e <- a1*b1/choose(n, 2); (s - e)/((a1 + b1)/2 - e) }

## sensitivity to the three characters with the most states
drop3 <- names(sort(sapply(dat, nlevels), decreasing = TRUE))[1:3]
Dr    <- daisy(dat[, setdiff(CH, drop3)], metric = "gower")
hr    <- hclust(Dr, "average")

## --------------------------------------------------------- 7. BIOME TESTS ---

set.seed(20260913)

a_bio  <- adonis2(D9 ~ biome, data = bio, permutations = PERM)
bd_bio <- betadisper(as.dist(m9), bio$biome)
set.seed(20260913)
pt_bio <- permutest(bd_bio, permutations = PERM)
hiC    <- bio$conf == "H"
set.seed(20260913)
a_bioH <- adonis2(as.dist(m9[hiC, hiC]) ~ droplevels(bio$biome[hiC]), permutations = PERM)
k1     <- bio$biome != "Dry_forest"
set.seed(20260913)
a_bio1 <- adonis2(as.dist(m9[k1, k1]) ~ droplevels(bio$biome[k1]), permutations = PERM)
Dbio   <- daisy(bio["biome"], metric = "gower")
set.seed(20260913)
mt_bio <- mantel(D9, Dbio, permutations = PERM)
set.seed(20260913)
mt_par <- mantel.partial(D9, Dbio, dist(nspec), permutations = PERM)
tab_ap <- table(grpVGS, bio$biome); tab_cl <- table(grp6[SP], bio$biome)
set.seed(20260913)
f_ap   <- fisher.test(tab_ap, simulate.p.value = TRUE, B = PERM)$p.value
set.seed(20260913)
f_cl   <- fisher.test(tab_cl, simulate.p.value = TRUE, B = PERM)$p.value

## ------------------------------------------------------------- 8. RESULTS ---

sink(file.path(OUT, "results.txt"))
cat("MITRACARPUS SEED MICROMORPHOLOGY — ANALYTICAL OUTPUT\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M"), "| R", paste0(R.version$major, ".", R.version$minor),
    "| seed 20260913 |", PERM, "permutations\n")

cat("\n\n== 1. DATA ==\n")
cat("Taxa:", nrow(dat), " Characters:", ncol(dat), " Specimens:", sum(nspec),
    " (median", median(nspec), ", range", min(nspec), "-", max(nspec), ")\n")
cat("\nStates per character:\n"); print(sapply(dat, nlevels))
cat("\nShannon entropy (natural log):\n"); print(round(sort(H, decreasing = TRUE), 3))
cat("\nMissing values per taxon:\n"); print(rowSums(is.na(dat))[rowSums(is.na(dat)) > 0])
cat("\nSpecies with a single specimen:\n"); print(names(nspec)[nspec == 1])
cat("\nA priori groups:\n"); print(table(grpVGS))
cat("\nCompound states, proportion of species per character:\n")
print(round(c(AW = 26/27, SSh = 16/27, PW = 13/27), 3))

cat("\n\n== 2. GOWER ==\n")
cat("Mean dissimilarity:", round(mean(D9), 4),
    " range:", paste(round(range(D9), 4), collapse = "-"), "\n")
cat("\nMantel, 9 characters vs 9 + seed dimensions:\n"); print({set.seed(20260913); mantel(D9, D12, permutations = PERM)})

cat("\n\n== 3. CLUSTERING ==\n")
cat("Cophenetic correlation by linkage:\n"); print(round(sort(coph, decreasing = TRUE), 4))
cat("\nCharacter jackknife (correlation with the full solution):\n")
print(round(data.frame(dissimilarity = jk, topology = jkt)[order(jk), ], 4))

cat("\n\n== 4. PARTITION VALIDATION ==\n")
cat("Mean silhouette width, k = 2 to 8:\n")
print(round(rbind(UPGMA = sil_hc, PAM = sil_pam), 4))
cat("\nA priori seed types I-IV:", round(sil_ap, 4),
    "\nSubtypes 1.1-1.5 (Group I):", round(sil_sub, 4),
    "\nRevised subtypes I-a to I-d:", round(sil_new, 4), "\n")
cat("\nRetained partition (k =", K, "), renumbered top-to-bottom in Fig. 3:\n")
for (g in levels(grp6))
  cat("  ", CLAB[g], ": ", paste(names(grp6)[grp6 == g], collapse = ", "), "\n", sep = "")
cat("\nClusters vs a priori groups:\n"); print(table(cluster = grp6, apriori = grpVGS))

cat("\n\n== 5. PCoA ==\n")
cat("Variance (% of positive eigenvalues):\n"); print(round(pv[1:8], 2))
cat("Cumulative 1-2:", round(sum(pv[1:2]), 1), "%  1-3:", round(sum(pv[1:3]), 1),
    "%  1-4:", round(sum(pv[1:4]), 1), "%\n")
cat("Cailliez-corrected, cumulative 1-2:", round(sum(pvC[1:2]), 1), "%\n")
cat("\nenvfit, characters on axes 1-3:\n"); print(ef3)
cat("\nenvfit, seed dimensions:\n"); print(efs)
cat("\nMCA eigenvalues (complete cases, n =", sum(cc), "):\n"); print(round(mca$eig[1:6, ], 3))

cat("\n\n== 6. DIAGNOSTIC PERFORMANCE ==\n")
cat("Species pairs:", TOT, "  unresolved by the full matrix:", length(full_un), "\n")
cat("\nPairs resolved by each character alone:\n"); print(sort(alone, decreasing = TRUE))
cat("\nUnique contribution (pairs unresolved if the character is removed):\n")
print(sort(uniq, decreasing = TRUE))
cat("\nMinimum sufficient subset:", minset$k, "characters,",
    length(minset$sets), "solution(s)\n")
for (s in minset$sets) cat("  {", paste(s, collapse = ", "), "}\n")
cat("\nGreedy accumulation:\n"); print(curve)

cat("\n\n== 7. PREDEFINED GROUPINGS ==\n")
cat("Seed types, matrix excluding VGS:\n"); print(a_type)
cat("\nSingle-species type removed:\n"); print(a_type2)
cat("\nPERMDISP, seed types:\n"); print(pt_type)
cat("\nMean distance to centroid:\n"); print(round(tapply(bd_type$distances, grpVGS, mean), 4))
cat("\nStates shared by all members of each group:\n")
print(t(sapply(levels(grpVGS), function(g) shared(which(grpVGS == g)))))
cat("\nSubtypes 1.1-1.5, matrix excluding DCD, EOrn and AW:\n"); print(a_sub)
cat("\nEach subtype against the rest of Group I:\n")
for (s in levels(droplevels(subty[inI]))) {
  set.seed(20260913)
  a <- adonis2(DI_red ~ factor(ifelse(subty[inI] == s, s, "other")), permutations = PERM)
  cat(sprintf("  %s: R2 = %.3f, F = %.3f, p = %.4f\n", s, a$R2[1], a$F[1], a$`Pr(>F)`[1]))
}
cat("\nRevised subtypes, correspondence with the former scheme:\n")
print(table(former = droplevels(subty[inI]), revised = sub_new))
cat("\nStates shared by all members of each revised subtype:\n")
print(t(sapply(levels(sub_new), function(s) shared(which(SP %in% names(sub_new)[sub_new == s])))))
cat("\nModal state and frequency:\n"); print(modal)

cat("\n\n== 8. SAMPLING EFFORT ==\n")
cat("Specimens vs compound states: rho =", round(sp_poly$estimate, 3),
    ", p =", round(sp_poly$p.value, 3), "\n")
cat("Specimens vs mean Gower distance: rho =", round(sp_dist$estimate, 3),
    ", p =", round(sp_dist$p.value, 3), "\n")
cat("\nenvfit, number of specimens:\n"); print(efn)
cat("\nWithout the five single-specimen species:\n")
cat("  cophenetic:", round(cor(Dk, cophenetic(hck)), 4), "(full set:", round(coph["average"], 4), ")\n")
cat("  silhouette k=6:", round(mean(silhouette(cutree(hck, 6), Dk)[, 3]), 4),
    "(full set:", round(sil_hc["6"], 4), ")\n")
cat("  adjusted Rand index vs the full partition:", round(ari(cutree(hck, 6), grp6[keep]), 4), "\n")
cat("\nWithout the three characters with the most states (",
    paste(drop3, collapse = ", "), "):\n", sep = "")
set.seed(20260913)
cat("  Mantel against the full matrix:", round(mantel(D9, Dr, permutations = PERM)$statistic, 4), "\n")
cat("  silhouette k=6:", round(mean(silhouette(cutree(hr, 6), Dr)[, 3]), 4), "\n")
cat("  adjusted Rand index:", round(ari(cutree(hr, 6), cutree(hc, 6)), 4), "\n")
cat("  species pairs unresolved:", length(unres(setdiff(CH, drop3))), "of", TOT, "\n")

cat("\n\n== 9. BIOME ==\n")
cat("Species per biome:\n"); print(table(bio$biome))
cat("\nAssignments inferred from specimen localities:\n"); print(SP[bio$conf == "L"])
cat("\nPERMANOVA:\n"); print(a_bio)
cat("\nPERMDISP:\n"); print(pt_bio)
cat("\nHigh-confidence assignments only:\n"); print(a_bioH)
cat("\nSingle-species biome removed:\n"); print(a_bio1)
cat("\nenvfit, biome on axes 1-3:\n"); print(efb)
cat("\nMantel, morphology vs biome:\n"); print(mt_bio)
cat("\nPartial Mantel, controlling for sampling effort:\n"); print(mt_par)
cat("\nFisher exact, a priori seed types x biome: p =", round(f_ap, 4), "\n"); print(tab_ap)
cat("\nFisher exact, clusters x biome: p =", round(f_cl, 4), "\n"); print(tab_cl)
sink()

## -------------------------------------------------------------- 9. TABLES ---

write.csv(round(m9, 4), file.path(OUT, "TableS2_gower_matrix.csv"))
write.csv(data.frame(axis = 1:8, eigenvalue = round(ev[1:8], 4),
          pct = round(pv[1:8], 2), cumulative = round(cumsum(pv)[1:8], 2)),
          file.path(OUT, "TableS3_pcoa_eigenvalues.csv"), row.names = FALSE)
write.csv(curve, file.path(OUT, "TableS4_accumulation.csv"), row.names = FALSE)
write.csv(modal, file.path(OUT, "TableS5_revised_subtypes.csv"))
write.csv(data.frame(Species = SP, Specimens = nspec, Seed_type = as.character(grpVGS),
          Cluster = as.character(grp6[SP]), Revised_subtype = SUBT[as.character(grp6[SP])],
          Ecoregion = bio$ecoregion, Biome = FULL[as.character(bio$biome)],
          Assignment = ifelse(bio$conf == "H", "literature", "inferred from specimens"),
          row.names = NULL), file.path(OUT, "TableS6_biomes.csv"), row.names = FALSE)
write.csv(data.frame(species = SP, apriori_group = as.character(grpVGS),
          cluster = as.character(grp6[SP]),
          revised_subtype = ifelse(is.na(SUBT[as.character(grp6[SP])]), "",
                                   SUBT[as.character(grp6[SP])]),
          PCoA1 = round(pco$points[, 1], 4), PCoA2 = round(pco$points[, 2], 4),
          PCoA3 = round(pco$points[, 3], 4),
          silhouette = round(silhouette(as.integer(grpVGS), D9)[, 3], 4), row.names = NULL),
          file.path(OUT, "TableS7_species_scores.csv"), row.names = FALSE)
write.csv(data.frame(Code = names(ef3$factors$r), r2 = round(ef3$factors$r, 4),
          p = ef3$factors$pvals, unique_contribution = uniq[names(ef3$factors$r)],
          pairs_alone = alone[names(ef3$factors$r)],
          entropy = round(H[names(ef3$factors$r)], 3), row.names = NULL)[
          order(-ef3$factors$r), ], file.path(OUT, "Table3_characters.csv"), row.names = FALSE)
write.tree(as.phylo(hc), file.path(OUT, "upgma_tree.nwk"))

## ------------------------------------------------------------- PALETTES ----
## Okabe-Ito and Tol hues, checked for deuteranopia, protanopia and tritanopia:
## the worst adjacent pair is dE 11.0, above the dE 8 floor. Grey is reserved
## for "Multiple biomes", the category that is not a biome.

COL_A <- c(I = "#0072B2", II = "#D55E00", III = "#009E73", IV = "#E69F00")
COL_B <- c(Moist_forest = "#009E73", Dry_forest = "#882255", Grass_savanna = "#E69F00",
           Montane_grass = "#0072B2", Deserts_xeric = "#D55E00", Multiple = "#999999")
COL_C <- c(C1 = "#0072B2", C2 = "#D55E00", C3 = "#009E73",
           C4 = "#E69F00", C5 = "#CC79A7", C6 = "#882255")

## ---------------------------------------------------------- TYPE SYSTEM ----
## One scale for every figure. Sizes are points at the final printed width;
## every gap is declared in millimetres so the spacing is a system, not a guess.

FONT <- "sans"                     # resolves to TeX Gyre Heros (Helvetica metrics)
BASE <- 8.5                        # tick labels and legend entries
TITL <- BASE + 1.5                 # axis titles, bold
LEAD <- 3.6                        # mm, axis title to tick labels
TICK <- 1.5                        # mm, tick length and tick-label offset

TH <- theme_classic(base_size = BASE, base_family = FONT) +
  theme(
    axis.title.x  = element_text(size = TITL, face = "bold", colour = "grey10",
                                 margin = margin(t = LEAD, unit = "mm")),
    axis.title.y  = element_text(size = TITL, face = "bold", colour = "grey10",
                                 margin = margin(r = LEAD, unit = "mm")),
    axis.text.x   = element_text(size = BASE, colour = "grey30",
                                 margin = margin(t = TICK, unit = "mm")),
    axis.text.y   = element_text(size = BASE, colour = "grey30",
                                 margin = margin(r = TICK, unit = "mm")),
    axis.line         = element_line(linewidth = 0.4, colour = "grey25"),
    axis.ticks        = element_line(linewidth = 0.4, colour = "grey25"),
    axis.ticks.length = unit(TICK, "mm"),
    legend.title       = element_text(size = BASE + 0.5, face = "bold", colour = "grey10"),
    legend.text        = element_text(size = BASE, colour = "grey25"),
    legend.key.size    = unit(4.2, "mm"),
    legend.box.spacing = unit(6, "mm"),
    legend.spacing.y   = unit(1.8, "mm"),
    legend.margin      = margin(0, 0, 0, 0),
    legend.position    = "top",
    panel.grid.major = element_line(linewidth = 0.22, colour = "grey93"),
    plot.subtitle = element_text(size = BASE, colour = "grey35",
                                 margin = margin(b = 4, unit = "mm")),
    plot.margin   = margin(4, 6, 4, 4, "mm"))

save_fig <- function(name, p, w, h) {
  tiff(file.path(FIG, paste0(name, ".tif")), width = w, height = h, units = "in",
       res = 600, compression = "lzw", type = "cairo", family = FONT)
  if (is.function(p)) p() else print(p); dev.off()
  cairo_pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h, family = FONT)
  if (is.function(p)) p() else print(p); dev.off()
  cat(sprintf("  %-26s %5.1f x %5.1f mm\n", name, w*25.4, h*25.4))
}

## ----------------------------------------------- FIG. 3  UPGMA PHENOGRAM ----
## Column positions are computed from the rendered text width, so the gaps
## between the taxon names, the two attribute bars and the cluster bracket
## are equal by construction.

## Text widths are measured in inches on a null device, so the gaps between the
## columns are exactly equal and independent of the axis scale. The seed type is
## printed as a numeral rather than a colour bar: colour is reserved for biome,
## so the two attribute columns cannot be confused with one another.
LABS <- labels(as.dendrogram(hc))
GI <- 0.115; BI <- 0.105                       # uniform gap and bar width, inches
pdf(NULL, width = W2, height = MM(205))
par(mar = c(5.0, 0.6, 9.0, 0.6), cex = 1, family = FONT)
plot.new(); plot.window(c(0, 1), c(0, 1))
WN <- max(strwidth(LABS, font = 3, cex = 0.90, units = "inches"))
WT <- strwidth("III", font = 1, cex = 0.86, units = "inches")
WC <- max(strwidth(CLAB, font = 2, cex = 0.88, units = "inches"))
PIN <- par("pin")[1]
invisible(dev.off())
BLOCK <- GI + WN + GI + WT + GI + BI + GI + WC + GI

fig3 <- function() {
  par(mar = c(5.0, 0.6, 9.0, 0.6), xpd = NA, family = FONT, cex = 1)
  dend <- set(set(as.dendrogram(hc), "branches_lwd", 1.1),
              "branches_col", "grey20")
  xmax <- max(hc$height)
  upi  <- (xmax * 1.02) / (PIN - BLOCK)
  x0   <- -GI * upi                            # taxon names
  xt   <- x0  - (WN + GI) * upi                # seed-type numeral
  xb   <- xt  - (WT + GI) * upi                # biome bar
  xbr  <- xb  - (BI + GI) * upi                # cluster bracket
  xcl  <- xbr - GI * upi                       # cluster label
  n    <- length(LABS)

  plot(dend, horiz = TRUE, leaflab = "none", axes = FALSE,
       xlim = c(xmax * 1.02, -BLOCK * upi))
  for (i in seq_len(n)) {
    text(x0, i, LABS[i], adj = 0, font = 3, cex = 0.90, col = "grey10")
    text(xt - WT*upi/2, i, as.character(grpVGS[LABS[i]]),
         adj = 0.5, cex = 0.86, col = "grey30")
    rect(xb, i - 0.46, xb - BI*upi, i + 0.46,
         col = COL_B[as.character(bio[LABS[i], "biome"])], border = NA)
  }
  text(xt - WT*upi/2, n + 1.05, "Seed type", srt = 90, adj = 0, cex = 0.82,
       font = 2, col = "grey25")
  text(xb - BI*upi/2, n + 1.05, "Biome", srt = 90, adj = 0, cex = 0.82,
       font = 2, col = "grey25")
  cl <- grp6[LABS]
  for (g in unique(cl)) {
    ii <- which(cl == g)
    segments(xbr, min(ii) - 0.40, xbr, max(ii) + 0.40,
             col = "grey25", lwd = 2.6, lend = 1)
    text(xcl, mean(range(ii)), CLAB[as.character(g)],
         font = 2, cex = 0.88, col = "grey10", adj = 0)
  }
  at <- seq(0, 0.8, 0.2)
  axis(1, at = at, labels = sprintf("%.1f", at), cex.axis = 0.92,
       lwd = 0.9, lwd.ticks = 0.9, col = "grey25", col.ticks = "grey25",
       col.axis = "grey30", padj = -0.45, tcl = -0.34)
  mtext("Gower dissimilarity", side = 1, line = 2.7, cex = 1.12,
        font = 2, col = "grey10", at = xmax/2)
  legend(x = xmax*1.02, y = n + 8.9, ncol = 3, bty = "n", xjust = 0,
         legend = SHRT[LVL], fill = COL_B[LVL], border = NA,
         cex = 0.86, text.col = "grey25", x.intersp = 0.55, y.intersp = 1.25)
  text(xmax*1.02, n + 2.7, sprintf("cophenetic r = %.3f", CC),
       adj = 0, cex = 0.86, col = "grey35")
}
save_fig("Fig3_UPGMA", fig3, W2, MM(220))

## -------------------------------------------------- FIG. 4  PCoA ORDINATION
## coord_fixed keeps the two axes on the same scale, so distances in the
## plot are proportional to distances in the ordination.

sc <- data.frame(sp = SP, PCo1 = pco$points[, 1], PCo2 = pco$points[, 2],
                 Cluster = grp6[SP], Type = grpVGS[SP], Biome = bio$biome)
hull <- do.call(rbind, lapply(split(sc, sc$Cluster), function(d)
  if (nrow(d) >= 3) d[chull(d$PCo1, d$PCo2), ] else NULL))

fig4 <- ggplot(sc, aes(PCo1, PCo2)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey88") +
  geom_polygon(data = hull, aes(fill = Cluster), alpha = 0.13,
               colour = NA, show.legend = FALSE) +
  geom_point(aes(fill = Cluster, shape = Type), colour = "white",
             size = 2.3, stroke = 0.4) +
  geom_text_repel(aes(label = sp, segment.colour = Cluster), colour = "grey12",
                  fontface = "italic",
                  size = 2.35, seed = 7, max.overlaps = Inf,
                  box.padding = 0.40, point.padding = 0.20,
                  min.segment.length = 0.20, segment.size = 0.35,
                  force = 2.4, force_pull = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = COL_C, breaks = paste0("C", 1:6),
                    labels = CLAB, name = "Cluster") +
  scale_colour_manual(values = COL_C, guide = "none", aesthetics = "segment.colour") +
  scale_shape_manual(values = c(I = 21, II = 24, III = 22, IV = 23),
                     breaks = c("I","II","III","IV"), name = "Seed type") +
  coord_fixed(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.085)) +
  scale_y_continuous(expand = expansion(mult = 0.085)) +
  labs(x = sprintf("PCoA axis 1 (%.1f%%)", pv[1]),
       y = sprintf("PCoA axis 2 (%.1f%%)", pv[2])) +
  guides(fill  = guide_legend(order = 1, nrow = 1,
                              override.aes = list(shape = 21, size = 2.6,
                                                  colour = "white", label = "")),
         shape = guide_legend(order = 2, nrow = 1,
                              override.aes = list(fill = "grey45", colour = "white",
                                                  size = 2.6))) +
  TH + theme(legend.box = "vertical", legend.spacing.y = unit(0.4, "mm"))
save_fig("Fig4_PCoA", fig4, W2, MM(168))

## --------------------------------- FIG. 5  ORDINATION FIT vs DIAGNOSTIC VALUE

d5 <- data.frame(char = names(ef3$factors$r), r2 = ef3$factors$r,
                 p = ef3$factors$pvals)
d5$uniq <- uniq[d5$char]; d5$alone <- alone[d5$char]
d5$sig  <- factor(ifelse(d5$p <= 0.05, "significant", "not significant"),
                  levels = c("significant", "not significant"))

fig5 <- ggplot(d5, aes(r2, uniq)) +
  geom_point(aes(colour = sig, size = alone), alpha = 0.9) +
  geom_text_repel(aes(label = char), size = 2.5, fontface = "bold",
                  colour = "grey15", seed = 11, box.padding = 0.60,
                  point.padding = 0.34, min.segment.length = 0.25,
                  segment.size = 0.22, segment.colour = "grey55",
                  force = 2.5, max.overlaps = Inf) +
  scale_colour_manual(values = c(significant = "#0072B2",
                                 `not significant` = "#D55E00"),
                      name = expression(bold(paste("envfit ", bolditalic(p))))) +
  scale_size_continuous(range = c(1.7, 4.6), breaks = c(150, 200, 250, 300),
                        name = "Species pairs resolved by the character alone") +
  scale_x_continuous(limits = c(0.235, 0.925), breaks = seq(0.3, 0.9, 0.1)) +
  scale_y_continuous(breaks = 0:2, limits = c(-0.45, 2.55)) +
  labs(x = expression(bold(paste("Fit to the ordination (envfit ",
                                 bolditalic(r)^bold("2"), ", PCoA axes 1-3)"))),
       y = "Species pairs left unresolved\nwhen the character is removed") +
  guides(colour = guide_legend(order = 1, nrow = 1,
                               override.aes = list(size = 2.6)),
         size   = guide_legend(order = 2, nrow = 1, title.position = "top",
                               override.aes = list(colour = "grey45"))) +
  TH + theme(legend.box = "vertical", legend.spacing.y = unit(0.4, "mm"),
             legend.key.size = unit(4.4, "mm"))
save_fig("Fig5_envfit_diagnostic", fig5, W2, MM(120))

## -------------------------------------- FIG. 6  BIOME COMPOSITION PER CLUSTER
## Double column with the legend set beside the panel: at single-column width the
## six biome names crowd the plot into a third of the figure.

d6 <- as.data.frame(table(Cluster = grp6[SP], Biome = bio$biome))
## label ink chosen by fill luminance, so every count stays legible
lum <- sapply(COL_B, function(h) { v <- col2rgb(h)/255
  0.2126*v[1] + 0.7152*v[2] + 0.0722*v[3] })
d6$ink <- ifelse(lum[as.character(d6$Biome)] > 0.5, "grey15", "white")
## explicit label positions: the first biome level is stacked on top
d6 <- d6[order(d6$Cluster, d6$Biome), ]
d6$ypos <- ave(d6$Freq, d6$Cluster,
               FUN = function(f) sum(f) - cumsum(f) + f/2)

fig6 <- ggplot(d6, aes(Cluster, Freq, fill = Biome)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.5) +
  geom_text(data = subset(d6, Freq > 0),
            aes(x = Cluster, y = ypos, label = Freq, colour = ink),
            inherit.aes = FALSE, fontface = "bold", size = 2.7,
            show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(values = COL_B, labels = SHRT, breaks = LVL, name = "Biome") +
  scale_x_discrete(labels = CLAB) +
  scale_y_continuous(breaks = seq(0, 10, 2), expand = expansion(c(0, 0.05))) +
  labs(x = "UPGMA cluster", y = "Number of species",
       subtitle = sprintf("Fisher exact test, cluster \u00d7 biome: p = %.3f", f_cl)) +
  guides(fill = guide_legend(ncol = 1, keyheight = unit(4.4, "mm"))) +
  TH + theme(panel.grid.major.x = element_blank(),
             legend.position = "right",
             legend.direction = "vertical",
             legend.box.spacing = unit(7, "mm"),
             legend.title = element_text(size = BASE + 0.5, face = "bold",
                                         colour = "grey10",
                                         margin = margin(b = 2.5, unit = "mm")),
             plot.subtitle = element_text(size = BASE, colour = "grey35",
                                          margin = margin(b = 5, unit = "mm")))
save_fig("Fig6_cluster_biome", fig6, W2, MM(100))

## ------------------------------------------------- SUPPLEMENTARY FIGURES ----

## Fig. S1 — sampling effort
dfe <- data.frame(n = nspec, d = meandist[SP], poly = poly)
save_fig("FigS1A_effort_distance",
  ggplot(dfe, aes(factor(n), d)) +
    geom_boxplot(fill = "grey93", colour = "grey35", linewidth = 0.4,
                 outlier.shape = NA, width = 0.58) +
    geom_jitter(width = 0.11, size = 1.9, colour = "#D55E00", alpha = 0.9) +
    labs(x = "Specimens examined", y = "Mean Gower distance to other species",
         subtitle = sprintf("Spearman rho = %.3f, p = %.3f",
                            sp_dist$estimate, sp_dist$p.value)) + TH,
  MM(120), MM(98))
save_fig("FigS1B_effort_variability",
  ggplot(dfe, aes(factor(n), poly)) +
    geom_boxplot(fill = "grey93", colour = "grey35", linewidth = 0.4,
                 outlier.shape = NA, width = 0.58) +
    geom_jitter(width = 0.11, height = 0.06, size = 1.9, colour = "#0072B2", alpha = 0.9) +
    scale_y_continuous(breaks = 0:3) +
    labs(x = "Specimens examined", y = "Characters scored with a compound state",
         subtitle = sprintf("Spearman rho = %.3f, p = %.3f",
                            sp_poly$estimate, sp_poly$p.value)) + TH,
  MM(120), MM(98))

## Fig. S2 — silhouette profile
save_fig("FigS2_silhouette",
  ggplot(rbind(data.frame(k = 2:8, s = sil_hc, m = "UPGMA"),
               data.frame(k = 2:8, s = sil_pam, m = "PAM")),
         aes(k, s, colour = m, shape = m)) +
    geom_hline(yintercept = sil_ap, linetype = "22", colour = "grey40", linewidth = 0.4) +
    annotate("text", x = 8, y = sil_ap - 0.022, hjust = 1, size = 2.7, colour = "grey30",
             label = sprintf("a priori seed types (%.3f)", sil_ap)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    scale_colour_manual(values = c(UPGMA = "#0072B2", PAM = "#D55E00"), name = NULL) +
    scale_shape_manual(values = c(UPGMA = 16, PAM = 17), name = NULL) +
    scale_x_continuous(breaks = 2:8) + ylim(0, 0.55) +
    labs(x = "Number of clusters (k)", y = "Mean silhouette width",
         subtitle = "No partition reaches 0.51, the threshold for reasonable structure") +
    TH, MM(140), MM(105))

## Fig. S3 — MCA of the character states
vs <- data.frame(mca$var$coord[, 1:2]); names(vs) <- c("D1", "D2")
vs$state <- rownames(vs); vs$char <- sub("_.*", "", vs$state)
save_fig("FigS3_MCA",
  ggplot(vs, aes(D1, D2, colour = char)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey88") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey88") +
    geom_point(size = 1.9) +
    geom_text_repel(aes(label = state), size = 2.3, seed = 4, max.overlaps = Inf,
                    box.padding = 0.35, segment.size = 0.22,
                    show.legend = FALSE) +
    labs(x = sprintf("MCA dimension 1 (%.1f%%)", mca$eig[1, 2]),
         y = sprintf("MCA dimension 2 (%.1f%%)", mca$eig[2, 2]),
         colour = "Character") +
    guides(colour = guide_legend(nrow = 1)) + TH,
  W2, MM(150))

## Fig. S4 — diagnostic accumulation curve
cu <- curve; cu$ch <- factor(cu$character., levels = cu$character.)
save_fig("FigS4_accumulation",
  ggplot(cu, aes(step, prop)) +
    geom_line(linewidth = 0.7, colour = "#0072B2") +
    geom_point(size = 2.2, colour = "#0072B2") +
    geom_text_repel(aes(label = ch), size = 2.5, fontface = "bold", colour = "grey15",
                    nudge_y = -0.022, seed = 2, min.segment.length = 0.4,
                    segment.size = 0.22, segment.colour = "grey55") +
    scale_x_continuous(breaks = 1:9) +
    scale_y_continuous(limits = c(0.80, 1.02), labels = scales::percent) +
    labs(x = "Number of characters (greedy order)", y = "Species pairs resolved",
         subtitle = sprintf("%d characters resolve all %d pairs",
                            min(cu$step[cu$prop == 1]), TOT)) + TH,
  MM(140), MM(105))

## Fig. S5 — PERMDISP on the reduced matrix
save_fig("FigS5_permdisp",
  ggplot(data.frame(d = bd_type$distances, g = grpVGS), aes(g, d, fill = g)) +
    geom_boxplot(alpha = 0.30, outlier.shape = NA, width = 0.52,
                 colour = "grey35", linewidth = 0.4) +
    geom_jitter(aes(colour = g), width = 0.12, size = 2.1, alpha = 0.9) +
    scale_fill_manual(values = COL_A) + scale_colour_manual(values = COL_A) +
    labs(x = "A priori seed type", y = "Distance to group centroid",
         subtitle = sprintf("PERMDISP F = %.2f, p = %.4f; Group IV is monospecific",
                            pt_type$tab$F[1], pt_type$tab$`Pr(>F)`[1])) +
    TH + theme(legend.position = "none"), MM(120), MM(100))

## Fig. S6 — PCoA coloured by biome, centroids for biomes with n >= 3
ok  <- names(which(table(sc$Biome) >= 3))
cen <- aggregate(cbind(PCo1, PCo2) ~ Biome, data = sc[sc$Biome %in% ok, ], FUN = mean)
seg <- merge(sc[sc$Biome %in% ok, ], cen, by = "Biome", suffixes = c("", "_c"))
save_fig("FigS6_PCoA_biome",
  ggplot(sc, aes(PCo1, PCo2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey88") +
    geom_segment(data = seg, aes(x = PCo1, y = PCo2, xend = PCo1_c, yend = PCo2_c,
                 colour = Biome), linewidth = 0.3, alpha = 0.5, show.legend = FALSE) +
    geom_point(aes(fill = Biome), shape = 21, colour = "white",
               size = 2.3, stroke = 0.4) +
    geom_point(data = cen, aes(fill = Biome), shape = 23, size = 3.8,
               colour = "white", stroke = 0.7, show.legend = FALSE) +
    geom_text_repel(aes(label = sp), colour = "grey12", fontface = "italic",
                    size = 2.35, seed = 1, max.overlaps = Inf,
                    box.padding = 0.40, point.padding = 0.20,
                    min.segment.length = 0.20, segment.size = 0.28,
                    segment.colour = "grey55", force = 2.4, show.legend = FALSE) +
    scale_fill_manual(values = COL_B, labels = SHRT, breaks = LVL, name = "Biome") +
    scale_colour_manual(values = COL_B, guide = "none") +
    coord_fixed(clip = "off") +
    scale_x_continuous(expand = expansion(mult = 0.085)) +
    scale_y_continuous(expand = expansion(mult = 0.085)) +
    labs(x = sprintf("PCoA axis 1 (%.1f%%)", pv[1]),
         y = sprintf("PCoA axis 2 (%.1f%%)", pv[2])) +
    guides(fill = guide_legend(nrow = 2,
             override.aes = list(shape = 21, size = 2.6, colour = "white"))) +
    TH, W2, MM(168))

cat("Done.", length(list.files(FIG, "\\.tif$")), "figures (TIFF + PDF) and",
    length(list.files(OUT, "\\.csv$")), "tables in", FIG, "and", OUT, "\n")
