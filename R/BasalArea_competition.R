##############################################################
# Function to convert tree size distribution in BA_comp metric
# Will Vieira
# March, 29, 2023
##############################################################


size_to_BAind <- function(
	mesh # size distribution frequency in milimiters
){
	# size to individual Basal area
  indBA_vec <- pi * ( mesh / 2 * 1e-3 )^2
  
	# return basal area in m²/ha
	return( indBA_vec )
}


size_to_BAplot <- function(
	N,
	plot_size
){
	# get BAind
	BAind <- size_to_BAind(N$meshpts)
	sum(BAind * N$Nvec * 1e4/plot_size)
}


#' Function to compute BA competition in function of population size vectors
#' N_intra: size distribution for focal species; output of `init_pop` function
#' N_inter: size distribution of competition species
#' If `N_inter` is NULL, BA_comp is from intraspecific competition. If `N_inter` is defined, then BA_comp is from interspecific competition
size_to_BAcomp <- function(
	N_intra,
	N_inter = NULL,
	plot_size
){
	# Define weather intra or inter competition is used
	if(is.null(N_inter)) {
		N_focal <- N_intra
	}else{
		N_focal <- N_inter
	}

	# size to individual basal area in m2
	baind = size_to_BAind(N_focal$meshpts)
	
	# basal in m2/ha corrected by the number of individuals
	ba_N_plot = baind * N_focal$Nvec * 1e4/plot_size

	sapply(
		N_intra$meshpts,
		function(x)
			sum(ba_N_plot[N_focal$meshpts > x])
	)
}


# Deprecated function
dbh_to_meshpts <- function(
	dbh,
	mesh
){
	# empty mshpts to count number of individuals per size class
	mesh_dbh <- rep(0, length(mesh))

	# approximate each individual size to its class
	dbh_class <- cut(dbh, breaks = mesh, labels = F)

	# count number of ind dbh at each class
	dbh_count <- table(dbh_class)

	# fill meshpts with the count of each class
	mesh_dbh[as.numeric(names(dbh_count))] <- dbh_count

	return( mesh_dbh)
}
