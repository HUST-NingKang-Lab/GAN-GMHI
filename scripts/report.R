H_signature <- data.frame(subset(all_matrix, all_matrix$PH_fold >= 1.4 & all_matrix$PH_diff >=8))
NH_signature <- data.frame(subset(all_matrix, all_matrix$PNH_fold >= 1.4 & all_matrix$PH_diff <= -8))

alpha_gmhi <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}

H_shannon <- apply((H_signature[,-c(4348:4350)]/100), 2, alpha_gmhi)
NH_shannon <- apply((NH_signature[,-c(4348:4350)]/100), 2, alpha_gmhi)

H_sig_count <- apply(H_signature[,-c(4348:4350)], 2, function(i) (sum(i > 0)))
NH_sig_count <- apply(NH_signature[,-c(4348:4350)], 2, function(i) (sum(i > 0)))

constant <- data.frame(cbind(H_sig_count,NH_sig_count))

HC1 <- constant[with(constant, order(-H_sig_count, NH_sig_count)), ]
H_constant <- length(rownames(H_signature))

NHC1 <- constant[with(constant, order(H_sig_count, -NH_sig_count)), ]
NH_constant <- length(rownames(NH_signature))

H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)

GMHI <- data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
Healthy_GMHI <- data.frame(GMHI[grep('Healthy', row.names(GMHI)),])
Healthy_GMHI$Phenotype<-"Healthy"

Nonhealthy_GMHI <- data.frame(GMHI[-grep('Healthy', row.names(GMHI)),])
Nonhealthy_GMHI$Phenotype<-"Nonhealthy"

colnames(Healthy_GMHI)[1] <- "GMHI"
colnames(Nonhealthy_GMHI)[1] <- "GMHI"

GMHI_20 <- data.frame(rbind(Healthy_GMHI, Nonhealthy_GMHI))
Healthy_accuracy <- sum(Healthy_GMHI$GMHI>0)*100/2636
Nonhealthy_accuracy <- sum(Nonhealthy_GMHI$GMHI<0)*100/1711
total_accuracy <- (Healthy_accuracy+Nonhealthy_accuracy)/2
report <- cbind(nrow(H_signature),nrow(NH_signature),Healthy_accuracy,Nonhealthy_accuracy,total_accuracy)
report