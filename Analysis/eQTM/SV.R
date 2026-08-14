library(data.table)
re<-fread("re_hierarchical_fdr_ACAT_4.txt")
re <- re[grepl("SV", re$id),]
re$id_first <- sub(":.*", "", re$id)
re_filtered <- re[
  , if (.N > 1 && uniqueN(id) > 1) .SD,
  by = id_first
]
re_filtered<-re_filtered[,c(1:5,25)]

### mixed significance: same SV segment, same transcipt, different coordinates ---> different effects
re_mixed <- re_filtered[
  , if (
    .N > 1 &&
    any(!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05) &&
    any(!is.na(E_meth_p_FDR) & E_meth_p_FDR > 0.05)
  ) .SD,
  by = .(id_first, id2)
]


length(unique(
  re_mixed$id_first
))/length(unique(re_filtered$id_first[!is.na(re_filtered$E_meth_p_FDR) & re_filtered$E_meth_p_FDR<0.05])) # 0.5357143

fwrite(re_mixed, file = "re_sv_mixedsignificance_fdr_ACAT_4.txt", sep = "\t", quote = FALSE, row.names = FALSE)
