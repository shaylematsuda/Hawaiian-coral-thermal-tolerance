suppressPackageStartupMessages(library(lmerTest))
out <- 'output/reviews/pam_t0_audit'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
x <- readLines('ACE21_physiology.Rmd')
y <- rep('', length(x)); active <- FALSE
for (i in seq_along(x)) {
  if (startsWith(x[i], '```{r')) { active <- !grepl('eval=FALSE', x[i], fixed=TRUE); next }
  if (startsWith(x[i], '```')) { active <- FALSE; next }
  if (active && i <= 3776) y[i] <- x[i]
}
e <- new.env(parent=globalenv())
exprs <- parse(text=y, keep.source=TRUE)
for (expr in exprs) tryCatch(invisible(eval(expr, e)), error=function(err) NULL)
stopifnot(exists('PAM_stats_T0_model', e))
fit <- e$PAM_stats_T0_model
write.csv(as.data.frame(anova(fit)), file.path(out, 'pam_t0_anova.csv'))
old <- e$PAM_stats_T0
old$Treatment[old$Frag.ID=='208'] <- 'Ambient'
old$Treatment[old$Frag.ID=='755'] <- 'High'
old$Treatment[old$Frag.ID=='N2'] <- 'Ambient'
old$Treatment[old$Frag.ID=='N35'] <- 'Ambient'
old$Treatment[old$Frag.ID=='N84'] <- 'High'
old$Treatment[old$Frag.ID=='1108'] <- 'High'
old_fit <- lmer(Yield_cu ~ Treatment * Species + (1|Tank.num), old)
write.csv(as.data.frame(anova(old_fit)), file.path(out, 'pam_t0_anova_before_correction.csv'))
saveRDS(fit, file.path(out, 'pam_t0_model.rds'))
write.csv(as.data.frame(model.frame(fit)), file.path(out, 'pam_t0_input.csv'), row.names=FALSE)
print(anova(fit))
