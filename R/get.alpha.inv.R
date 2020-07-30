get.alpha.inv <-
function(chi, W, l, rho, alpha, overdisperse="no", phi=1, od.alpha=0) {
	t = if (rho < 1) { qnorm(alpha) } else { qnorm(1-alpha) }
	if (overdisperse %in% c("ql", "quasi-likelihood")) { t = t*sqrt(phi) }
	c = W*chi/l
	# k = c/(c-t^2)
	# r = k + sqrt(k^2-1)
	oldWarn = getOption("warn")
	options(warn=-1)
	if (overdisperse %in% c("nb", "negative binomial")) {
		A = 1 + c*od.alpha
		r = (c + t*sqrt(2*c*A - A^2*t^2))/(c-A*t^2)
	} else {
		r = (c + t*sqrt(2*c-t^2))/(c-t^2) 
	}
	options(warn=oldWarn)
	if (any(r < 0 | is.nan(r))) { r[r<0] = NaN }
	return(r)
}

#get a cutoff r(ρ) which yields desired minimum spec/sens for testing a particular copy-number ratio ρ at a particular exon with some depth and length.
#qnorm is inverse operation of pnorm, probability distribution function.
