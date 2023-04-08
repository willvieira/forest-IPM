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


BAind_to_BAplot <- function(
	BAind,
	N,
	plot_size
){
	sum(BAind * N * 1e4/plot_size)
}


size_to_BAcomp <- function(
	mesh,
	N,
	plot_size
){
	# size to individual basal area in m2
	baind = size_to_BAind(mesh)
	
	# basal in m2/ha corrected by the number of individuals
	ba_N_plot = baind * N * 1e4/plot_size
	
	sapply(
		mesh,
		function(x)
			sum(ba_N_plot[mesh > x])
	)
}


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
