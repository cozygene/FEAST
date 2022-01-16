
#' @import glmnet
#' @import CVXR
#' @export
STENSL <- function(
    C,
    metadata,
    EM_iterations,
    COVERAGE,
	l.range=c(0.01, 0.1, 1, 10, 100),
	method='stensl'
) {
    max.lambda <- NA
	results <- list()
	for (li in 1:length(l.range)) {
		sparse.lambda <- l.range[li]
        pij <<- NA
		set.seed(0)
		em_result <- FEAST(
			C=C,
			metadata=metadata,
			EM_iterations=EM_iterations,
			# rarefy_sink=T,
			COVERAGE=COVERAGE,
			different_sources_flag=0,
			method=method,
			options=list(sparse.lambda=sparse.lambda, pij=NA)
		)
		cat('\n')

		results[[li]] <- list(
			lambda=sparse.lambda,
			em=em_result
		)
	}

	max.ll <- results[[1]]$ll[length(results[[1]]$ll)]
	max.ind <- 1
	for (ii in 1:length(results)) {
		blob <- results[[ii]]
		print(paste(ii,
			blob$ll[length(blob$ll)],
			unlist(blob$em$proportions_mat)[51],
			blob$lambda))
		if ( max.ll < blob$ll[length(blob$ll)]) {
			max.ll <- blob$ll[length(blob$ll)]
			max.ind <- ii
		}
	}
	print(paste('Choosing:', max.ind, max.ll))
    max.lambda <<- l.range[max.ind]

	feast_result <- list(
		proportions_mat=results[[max.ind]]$em$proportions_mat,
        max.lambda=max.lambda
	)

    return(feast_result)
}