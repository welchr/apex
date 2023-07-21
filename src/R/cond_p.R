loop_gene_tryCatch <- function (x, rlist=NULL, meta=NULL) {
  if (missing(rlist) || is.null(rlist)) {
    rlist = res_list
  }
  if (missing(meta) || is.null(meta)) {
    meta = ssm
  }
  print(x)
  return(
    tryCatch(
      loop_gene(gene=names(rlist)[x],c_list=rlist[[x]],heter=T,meta=meta),
      error=function(e) {
        warning(e)
        data.frame("gene"=names(rlist)[x])
      }
    )
  )
}

loop_gene=function(gene,c_list,heter=T,meta=NULL){
  if (missing(meta) || is.null(meta)) {
    meta = ssm
  }
  meta_object = meta$getSuffStats(gene)
  meta_V <- meta_object$XtX
  meta_U <- meta_object$Xty
  meta_S <- meta_object$yty
  meta_n <- meta_object$n_samples
  meta_m <- meta_object$n_covariates
  gene_sf_het <- getStudySuffStats(gene, meta)

  all_out=data.frame()
  for (i in 1:length(c_list)){
    not_c=c_list[i]
    c_condi = c_list[c_list != not_c]
    c_idx=match(c(c_condi),rownames(meta_V))
    out=data.frame()
    for (k_idx in seq(1:nrow(meta_V))[-c_idx]) {
      meta_V_ck = meta_V[c_idx,k_idx]
      meta_V_cc = meta_V[c_idx,c_idx]
      meta_V_kk = meta_V[k_idx,k_idx]

      RSQ = as.numeric((t(meta_V_ck) %*% solve(meta_V_cc) %*% meta_V_ck)/meta_V_kk)
      ok = !any(is.na(meta_V_ck))

      if(ok & (RSQ <= 0.8) & heter) {
        k_out=heter_cond_p(gene_sf_het,c_idx,k_idx)
      } else if (ok & (RSQ <= 0.8) & !heter){
        k_out=cond_p(meta_V,meta_U,meta_S,meta_n,meta_m,c_idx,k_idx)
      } else {
        U_star = 0
        V_star = 0
        S_star = 0
        beta_star = 0
        se_star = 0
        df_star = 0
        pval = 1
        k_out=data.frame(U_star,V_star,S_star,beta_star,se_star,
                         df_star,pval,row.names = rownames(meta_V)[k_idx])
      }
      out=rbind(out,k_out)
    }
    out$not_c = not_c
    out$gene = gene
    all_out=rbind(all_out,out)
  }
  return(all_out)
}

cond_p=function(V,U,S,n,m,c_idx,k_idx){
  V_ck = V[c_idx,k_idx]
  V_cc = V[c_idx,c_idx]
  V_kk = V[k_idx,k_idx]
  U_star = U[k_idx]-t(V_ck) %*% solve(V_cc) %*% U[c_idx]
  V_star = V_kk-t(V_ck) %*% solve(V_cc) %*% V_ck
  S_star = S-t(U[c_idx]) %*% solve(V_cc) %*%U[c_idx]
  beta_star <- U_star/V_star
  df_star=n - m -length(c_idx) - 1
  se_star <- sqrt( (S_star/V_star -beta_star^2  )/df_star )
  pval = pf( (beta_star/se_star)^2, 1, df_star, low =FALSE )
  out=data.frame(U_star,V_star,S_star,beta_star,se_star,
                 df_star,pval,row.names = rownames(V)[k_idx])
  return(out)
}

heter_cond_p=function(obj,c_idx,k_idx){
  gene_in_study=unlist(lapply(1:length(obj),function(x){!is.na(obj[[x]]$n_covariates)}))
  obj=obj[gene_in_study]

  tmp=lapply(obj,function(x){heter_cond_p_within_study(x,c_idx,k_idx)})
  tmp=do.call(rbind,tmp)
  colnames(tmp)=c('U_i','V_i','S_i','beta_i','se_i','df_i','weight_i','p_i')
  U_star = sum(tmp[,'U_i']*tmp[,'weight_i'])
  V_star = sum(tmp[,'V_i']*tmp[,'weight_i'])
  S_star = sum(tmp[,'S_i']*tmp[,'weight_i'])
  beta_star <- U_star/V_star
  df_star = sum(tmp[,'df_i'])
  se_star <- sqrt( (S_star/V_star -beta_star^2  )/df_star )
  pval = pf( (beta_star/se_star)^2, 1, df_star, low =FALSE )
  return(data.frame(U_star,V_star,S_star,beta_star,se_star,df_star,pval,row.names = rownames(tmp)[1]))
}

heter_cond_p_within_study=function(x,c_idx,k_idx){
  V <- x$XtX
  U <- x$Xty
  S <- x$yty
  n <- x$n_samples
  m <- x$n_covariates

  V_ck = V[c_idx,k_idx]
  V_cc = V[c_idx,c_idx]
  V_kk = V[k_idx,k_idx]

  U_star = U[k_idx]-t(V_ck) %*% solve(V_cc) %*% U[c_idx]
  V_star = V_kk-t(V_ck) %*% solve(V_cc) %*% V_ck
  S_star = S-t(U[c_idx])%*% solve(V_cc)%*%U[c_idx]
  beta_star <- U_star/V_star
  df_star=n - m -length(c_idx) - 1
  #weight_star = df_star/S_star
  #weight_star = df_star/S
  weight_star = df_star/(S_star - (U_star^2)/(V_star))
  se_star <- sqrt( (S_star/V_star -beta_star^2  )/df_star )
  pval = pf( (beta_star/se_star)^2, 1, df_star, low =FALSE )
  return(data.frame(U_star,V_star,S_star,beta_star,se_star,df_star,weight_star,pval,row.names = rownames(V)[k_idx]))
}

getStudySuffStats <- function(gene, object){

   use_studies <- which(sapply(object$study_data, function(x) gene %in% x$genes))

  lapply(object$study_data[use_studies], function(x){
    x$getSuffStats(gene)
  })
}

