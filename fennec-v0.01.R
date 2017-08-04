##########################################################################################
## Requirements/Settings

library(jsonlite)
library(phyloseq)
dbversion = "1.0"
baseurl = "http://fennec.molecular.eco/api/"

##########################################################################################
## get list of traits
list_traits <- fromJSON("http://fennec.molecular.eco/api/listing/traits?dbversion=1.0&limit=1000")
# list_orgs <- fromJSON("http://fennec.molecular.eco/api/listing/organisms?dbversion=1.0&limit=10")
# list_sample <- fromJSON("http://fennec.molecular.eco/api/details/traits?dbversion=1.0&trait_type_id=1&fennec_ids[]=122277&fennec_ids[]=667818")
# name2fennid <- fromJSON("http://fennec.molecular.eco/api/mapping/byOrganismName?dbversion=1.0&ids%5B%5D=Bellis%20perennis&ids%5B%5D=Brassica%20napus")
# http://fennec.molecular.eco/api/details/traitsOfOrganisms?dbversion=1.0&fennec_ids[]=122277&fennec_ids[]=667818

##########################################################################################
## mapping scientific names to fennec IDs
name2fennecID <- function(names){

	# prepare the scientific names for the json call	
	names = unlist(lapply(example_names, URLencode))
	names = unlist(lapply(names, function(x) paste("&ids%5B%5D=",x, sep="") ))
	
	# prepare the call
	orgCall = paste(names, collapse="")
	apiCall = paste(baseurl, "mapping/byOrganismName?dbversion=", dbversion, orgCall, sep="")
	
	# perform the call and return the result
	idlist <- fromJSON(apiCall)
	return(idlist)	
}

##########################################################################################
## mapping trait names to fennec IDs
trait2traitID <- function(traits){
	traits  = unlist(lapply(traits, function (x) list_traits[list_traits$name == x,"trait_type_id"]))

	## wahts the API for traits by name?
	#traits = unlist(lapply(example_names, URLencode))

}

##########################################################################################
## FENNEC mapping of a OTU table in vegan format
fennec4vegan <- function(mat){
	taxa = rownames(mat) 
	fennecIDs 	= as.character(unlist(name2fennecID(taxa)))
	##traitIDs 	= as.character(unlist(trait2traitID(traits)))
	
	idCallC = unlist(lapply(fennecIDs, function(x) paste("&fennec_ids[]=",x, sep="") ))
	idCall  = paste(idCallC, collapse="")
	
	# get all traits
	apiCallB = paste(baseurl, "details/traitsOfOrganisms?dbversion=", dbversion, idCall, sep="")

	traitlist <- fromJSON(apiCallB)
	#traitlist= traitlist[12]
	fmat = matrix(ncol=length(traitlist),nrow=length(fennecIDs))
	colnames(fmat) <- names(traitlist)
	rownames(fmat) <- fennecIDs
	
	#traitdetails=list()
#	for (i in 1:length(names(traitlist))){
	for (i in 1:length(names(traitlist))){

		apiCall = paste(baseurl, "details/traits?dbversion=", dbversion, "&trait_type_id=",names(traitlist)[i], idCall, sep="")
		cat("\n", apiCall)
		tmptraits <- fromJSON(apiCall)
		#cat(unlist(tmptraits))
		#traitdetails[[names(traitlist)[i]]] <- fromJSON(apiCall)
		
		for (j in 1:length(tmptraits$values)){
			
			
				if (tmptraits$trait_format == "categorical_free"){
					for (k in 1:length(tmptraits$values[[names(tmptraits$values)[j]]])){
						value 	= names(tmptraits$values)[j]
						x 		= as.character(tmptraits$values[[names(tmptraits$values)[j]]][k])			
						fmat[x, names(traitlist)[i]] = value					
					}			
				}else{
					if (length(tmptraits$values[[names(tmptraits$values)[j]]])>0){
						value 	= as.character(tmptraits$values[[names(tmptraits$values)[j]]][1]) ### taking only the first value atm
						x 		= names(tmptraits$values)[j]
						fmat[x, names(traitlist)[i]] = value	
					}		
			}
			
#			if (length(tmptraits$values[[names(tmptraits$values)[j]]])>0){ # needs to sort out what to do with more than 1 value
				#fmat[x, names(traitlist)[i]] = value			
			#}else{
				#fmat[x, names(traitlist)[i]] = "DD"
#			}
#			} # end k
		
		} # end j
	} # end i

#	idCall = unlist(lapply(fennecIDs, function(x) paste("&fennec_ids[]=",x, sep="") ))
#	idCall = paste(idCall, collapse="")
#	apiCall = paste(baseurl, "details/traits?dbversion=", dbversion, "&trait_type_id=1&trait_type_id=2", idCall, sep="")
#	cat(apiCall)	
#	traitlist <- fromJSON(apiCall)
	
	#lapply(rownames(mat), function(x) fennecIDs[[x]])
	
	# nice up the names
	colnames(fmat) <- lapply(colnames(fmat), function(x) gsub(" ",".",as.character(traitlist[[x]]["trait_type"])))
	rownames(fmat) <- taxa
	
	return(fmat)	
}


##########################################################################################
## FENNEC mapping of a OTU table in phyloseq format
fennec4phyloseq <- function(phyloseq){
	fennec4vegan(t(otu_table(phyloseq)))
}


##########################################################################################
## Example dataset
example_names 		= c("Bellis perennis","Brassica napus","Ranunculus acris")
example_otutable 	= matrix(nrow=length(example_names),ncol=10,data=sample(c(rep(0,25),1:10),length(example_names)*10,replace=T))
rownames(example_otutable) = example_names
example_traits 		= c("Plant Growth Habit","Leaf Area")

f = fennec4vegan(example_otutable)
f
