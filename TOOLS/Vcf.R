library('stringr')
library('itertools')
library('foreach')

vcfFile <- setRefClass(
    "vcfFile",
    field = list(
        reference = 'character',
        samples = 'character',
        Items = 'list'
    )
)

vcfRecord <- setRefClass(
    "vcfRecord",
    field = c(
		'samples',
		'CHROM',
		'POS',
		'ID', 
		'REF',
		'ALT',
		'QUAL',
		'FILTER',
		'INFO',
		'FORMAT',
		'SAMPLE',
		'start',   #Maf format 1-based
		'end',     #Maf format 1-based end inclusive
		'type'
    ),
    method = list(
        initialize = function(record, samples){   #这里的等号=不可以改成<-
            Record <- as.vector(str_split(record, "\\s+", simplify = T))
            .self$samples <- samples
            .self$CHROM <- Record[1]
            .self$POS <- as.numeric(Record[2])
            .self$ID <- Record[3]
            .self$REF <- Record[4]
            .self$ALT <- as.vector(str_split(Record[5], ',', simplify = T))
            .self$QUAL <- Record[6]
            .self$FILTER <- Record[7]
            .self$INFO <- getInfo(Record[8])
            .self$FORMAT <- as.vector(str_split(Record[9], ':', simplify = T))
            .self$SAMPLE <- getSample(Record[10:length(Record)])
            .self$start <- .self$POS
            .self$end <- .self$POS + length(.self$REF) - 1
            .self$type <- getType()
        },
        getInfo = function(info){
            tempINFO <- list()
            for (eachItem in as.vector(str_split(info, ':', simplify = T))){
                if (!str_detect(eachItem, '=')){
                    next
                }
                tempEachItem <- as.vector(str_split(eachItem, '=', simplify = T))
                tempINFO[tempEachItem[1]] <- tempEachItem[2]
            }
            return(tempINFO)
        },
        getSample = function(format){
            if (length(format) != length(.self$samples)){
                stop('ERROR: NOT EQUAL -- number of samples vs FORMAT!')
            }
            tempSample <- list()
            for (eachFormat in as.list(enumerate(format))){   #enumerate返回的迭代器for循环无法操作，只能先变成list，迭代器可以和foreach很好的配合
                if (!.self$samples[eachFormat$index] %in% tempSample){
                    tempSample[[.self$samples[eachFormat$index]]] <- list()   #务必注意对于list []和[[]] 的显著区别, 也可以考虑用$
                }else{
                    stop('ERROR: SAME sample name for different FORMAT!')
                }
                eachCall <- as.vector(str_split(eachFormat$value, ':', simplify = T))
                if (length(.self$FORMAT) != length(eachCall)){
                    stop('ERROR: FORMAT terms deficiency')
                }
                for (eachCallItem in as.list(enumerate(eachCall))){
                    if (!.self$FORMAT[eachCallItem$index] %in% tempSample[[.self$samples[eachFormat$index]]]){
                        tempSample[[.self$samples[eachFormat$index]]][[.self$FORMAT[eachCallItem$index]]] <- eachCallItem$value
                    }else{
                        stop('ERROR: duplicated FORMAT term!')
                    }
                }

            }
            return(tempSample)
        },
        getType = function(){
            if (length(.self$REF) == 1){
                for (eachAlt in .self$ALT){
                    if (length(eachAlt) > 1){
                        return('INS')
                    }
                }
                return('SNP')
            }else{
                for (eachAlt in .self$ALT){
                    if (length(eachAlt) >= length(.self$REF)){
                        return('INS')
                    }else{
                        return('DEL')
                    }
                }
            }
        }
    )
)


vcfReader <- function(vcf){
    varSet <- vcfFile$new()
    flag <- FALSE
    con <- file(vcf, open = 'r')
    foreach(eachLine = ireadLines(con, warn = FALSE), .combine = c) %do% {
        if (startsWith(eachLine, '##reference')){
            varSet$reference <- str_replace(str_trim(eachLine), '##reference=file://', '')
        }
        if (startsWith(eachLine, '#CHROM')){
            temp <- as.vector(str_split(str_trim(eachLine), "\\s+", simplify = T))
            varSet$samples <- temp[which(temp == 'FORMAT') + 1 ：length(temp)]
            flag <- TRUE
        }
        if (!startsWith(eachLine, '#') && flag){
            varSet$Items <- c(varSet$Items, vcfRecord$new(str_trim(eachLine), varSet$samples))
        }
        NULL   #让foreach不返回任何东西
    }
    close(con)
    return(varSet)
}