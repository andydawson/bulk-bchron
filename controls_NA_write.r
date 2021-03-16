# go through files 

fnames = list.files('Cores-new', '*.csv', recursive=TRUE)

fname_str = sapply(fnames, function(x) paste0('Cores-new/', x))

for (fname in fname_str){
  
  geochron <- try(read.table(fname, sep=',', header=TRUE), silent=TRUE)
  if(is(geochron, "try-error")) {
    next
  }
  
  if (any(is.na(geochron$limitolder) | any(is.na(geochron$limityounger)))) {
    print(fname)
    print(geochron)
    
    na_rc = which((is.na(geochron$limityounger)| is.na(geochron$limitolder)) & (geochron$type == "Radiocarbon"))
    na_ct = which((is.na(geochron$limityounger)| is.na(geochron$limitolder)) & (geochron$type == "Core top"))
    na_teph = which((is.na(geochron$limityounger)| is.na(geochron$limitolder)) & (geochron$type == "Tephra"))
    
    if (length(na_rc) > 0) {
      new_age = data.frame(age=geochron$age[na_rc])
      model_p = cbind(new_age, data.frame(predict.gam(mod, new_age, type='response', se.fit=TRUE)))
      younger = geochron$age[na_rc] - (crit.t * model_p$se.fit)
      older   = geochron$age[na_rc]  + (crit.t * model_p$se.fit)

      geochron$limityounger[na_rc] =  younger
      geochron$limitolder[na_rc]   =  older
      
      geochron$error[na_rc] = crit.t * model_p$se.fit

    }
    
    if (length(na_ct) > 0) {
      geochron$limityounger[na_ct] = ifelse((geochron$age[na_ct] - ct_error) > -70, geochron$age[na_ct] - ct_error, -70)
      geochron$limitolder[na_ct]   = geochron$age[na_ct] + ct_error
      
      geochron$error[na_ct] = ct_error
      
    }
    
    if (length(na_teph) > 0) {
      geochron$limityounger[na_teph] = geochron$age[na_teph] - tephra_error
      geochron$limitolder[na_teph]   = geochron$age[na_teph] + tephra_error
      
      geochron$error[na_teph] = tephra_error
      
    }
    
    
  }
  
}
