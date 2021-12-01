
# some global vars
pij <- NA
sparse.lambda <- NA
ll.list <- c()

STENSL <- function(
    C,
    metadata,
    EM_iterations,
    COVERAGE,
	l.range=c(0.01, 10)
	# l.range=c(0.01, 0.1, 1, 10, 100)
) {
    max.lambda <- NA
	results <- list()
	for (li in 1:length(l.range)) {
		sparse.lambda <<- l.range[li]
		ll.list <<- c()
        pij <<- NA
		set.seed(0)
		em_result <- FEAST(
			C=C,
			metadata=metadata,
			EM_iterations=EM_iterations,
			# rarefy_sink=T,
			COVERAGE=COVERAGE,
			different_sources_flag=0,
			method='stensl'
		)
		cat('\n')

		results[[li]] <- list(
			lambda=sparse.lambda,
			ll=ll.list,
			alpha=unlist(em_result$proportions_mat)
		)
	}

	max.ll <- results[[1]]$ll[length(results[[1]]$ll)]
	max.ind <- 1
	for (ii in 1:length(results)) {
		blob <- results[[ii]]
		print(paste(ii, blob$ll[length(blob$ll)], unlist(blob$alpha)[51], blob$lambda))
		if ( max.ll < blob$ll[length(blob$ll)]) {
			max.ll <- blob$ll[length(blob$ll)]
			max.ind <- ii
		}
	}
	print(paste('Choosing:', max.ind, max.ll))
    max.lambda <<- l.range[max.ind]

	feast_result <- list(
		alpha=results[[max.ind]]$alpha,
		iterations=results[[max.ind]]$iterations,
        max.lambda=max.lambda
	)

    return(feast_result)
}