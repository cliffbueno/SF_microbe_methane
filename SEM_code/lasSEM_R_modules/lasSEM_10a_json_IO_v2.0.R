
# make lasSEM objects JSON
export_lasSEM_as_json = function(lasSEM_object, path = ''){
    
    # create nested json
    lasSEM_list <-list(lasSEM_object$filt_models, lasSEM_object$sem_models)
    names(lasSEM_list) <- c('filt_models', 'sem_models')
    json <- toJSON(lasSEM_list)
    
    # get filename and save
    name <- deparse(substitute(CH4ac_mod.4_sem))
    outfile <- paste0(path, name, ".json")
    write(json, outfile) 
}
################################################################################
# JSON write, read for lasSEM objects 
# (resid. filtered models; SEM models, R2 filtered)

# for notebook running speed
# take lasSEM output for model branches
# pass to composite branch testing w/o re-running large lasSEM model sets


# a) writes lasSEM as json -- e.g. -- export_lasSEM_as_json(CH4ac_mod.4_sem)
# b) reads json to lasSEM -- e.g. -- import_lasSEM_from_json('CH4ac_mod.4_sem.json')


library('jsonlite')

################################################################################
## a)  make lasSEM objects JSON & write file

# TODO: pass in features?

export_lasSEM_as_json = function(lasSEM_object, path = ''){

    # create nested json
    lasSEM_list <-list(lasSEM_object$filt_models, lasSEM_object$sem_models)
    names(lasSEM_list) <- c('filt_models', 'sem_models')
    json <- toJSON(lasSEM_list)

    # get filename and save
    name <- deparse(substitute(CH4ac_mod.4_sem))
    outfile <- paste0(path, name, ".json")
    write(json, outfile) 

}

# demonstrate function
# js_test <- export_lasSEM_as_json(CH4ac_mod.4_sem)
# js_test


################################################################################
## b) get lasSEM objects from JSON

# TODO: add path?

import_lasSEM_from_json = function(lasSEM_json){
    lasSEM_object <- fromJSON(lasSEM_json)
    names(lasSEM_object) <- c('filt_models', 'sem_models')
    lasSEM_object
}

# demonstrate function
#js_read <- import_lasSEM_from_json(js_test)
#js_read <- import_lasSEM_from_json('CH4ac_mod.4_sem.json')
#js_read