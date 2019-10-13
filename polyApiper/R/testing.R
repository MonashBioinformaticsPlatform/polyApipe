
#' Confident-correlation gene-set enrichment
#'
#' Effect size is the Pearson correlation between the indicator variable for the gene-set and the statistic.
#'
#' @export
enrichment <- function(stats, orgdb, confident_sign, min_size=2, ad=FALSE) {
    gene_term <- AnnotationDbi::select(
            orgdb, keys=names(stats), keytype="ENSEMBL", columns="GOALL") %>%
        select(name=ENSEMBL, GOID=GOALL) %>%
        filter(!is.na(GOID)) %>%
        distinct(name, GOID) # Duplicate entries, by evidence type?

    go_info <- AnnotationDbi::select(
        GO.db, keys=unique(gene_term$GOID), keytype="GOID",
        columns=c("GOID","TERM","ONTOLOGY"))

    #all_m <- mean(stats)
    #all_s2 <- sum((stats-all_m)^2)

    # R^2 from using set as a predictor in a linear model
    if (ad) {
        median_all <- median(stats)
        ad_all <- sum(abs(stats-median_all))
        doit <- function(set_names) {
            l <- names(stats) %in% set_names
            am <- median(stats[l])
            bm <- median(stats[!l])
            r <- 1-(sum(abs(stats[l]-am))+sum(abs(stats[!l]-bm)))/ad_all
            data.frame(
                effect=r,
                diff=am-bm)
        }
    } else {
        doit <- function(set_names) {
            l <- names(stats) %in% set_names
            r <- cor(stats, l)
            if (is.na(r))
                r <- 0
            #se <- sqrt( (1-r*r)/(length(stats)-2) )

            am <- mean(stats[l])
            bm <- mean(stats[!l])
            data.frame(
                effect=r,
                #se=se,
                #df=length(stats)-2,
                diff=am-bm)
        }
    }

    result <- gene_term %>%
        mutate(stat=stats[name]) %>%
        group_by(GOID) %>%
        # A gene set containing all genes can be discarded:
        filter(length(name) >= min_size, length(name)+min_size <= length(stats)) %>%
        summarize(data=list(doit(name)), size=length(name), names=list(name)) %>%
        unnest(data) %>%
        arrange(-abs(effect)) %>%
        left_join(go_info, by="GOID")

    seen <- c()
    novel <- rep(0,nrow(result))
    up <- rep(0,nrow(result))
    down <- rep(0,nrow(result))
    novel_up <- rep(0,nrow(result))
    novel_down <- rep(0,nrow(result))
    for(i in seq_len(nrow(result))) {
        this_names <- result$names[[i]]
        this_novel <- setdiff(this_names, seen)
        novel[i] <- length(this_novel)
        up[i] <- sum(confident_sign[this_names] > 0)
        down[i] <- sum(confident_sign[this_names] < 0)
        novel_up[i] <- sum(confident_sign[this_novel] > 0)
        novel_down[i] <- sum(confident_sign[this_novel] < 0)

        seen <- c(seen, this_novel)
    }

    result$up <- up
    result$down <- down
    result$novel <- novel
    result$novel_up <- novel_up
    result$novel_down <- novel_down

    result %>%
        select(r=effect, diff, size, up, down, novel, novel_up, novel_down,
            id=GOID, term=TERM, ontology=ONTOLOGY, gene_ids=names)

    # On reflection, this doesn't seem to make much sense:
    #
    #top <- normal_confects(result$effect, result$se, result$df)
    #top$table <- bind_cols(top$table,
    #    select(result, diff, size, goid=GOID, term=TERM, ontology=ONTOLOGY, gene_ids=names)[top$table$index,,drop=F])
    #top$effect_desc <- "correlation"
    #top
}


#' @export
do_test <- function(out_dir, weitrix, organism, design, coef, confect=FALSE, rank=FALSE, ad=FALSE, fdr=0.05) {
    weitrix <- load_banquet(weitrix)
    organism <- load_banquet(organism)

    #TODO: do this better
    #TODO: comp can be reused when testing found components
    comp <- weitrix_components(weitrix, p=0, design=design)
    cal_weitrix <- weitrix_calibrate_trend(weitrix, comp)
    fit <- lmFit(weitrix_elist(cal_weitrix), design)
    top <- limma_confects(fit, coef, full=TRUE, fdr=fdr)

    ensure_dir(out_dir)
    saveRDS(top, file.path(out_dir, "confects.rds"))
    write_csv(top$table, file.path(out_dir, "confects.csv"))

    for(signed in c(F,T)) {
        if (confect)
            stats <- top$table$confect
        else if (rank)
            stats <- (nrow(top$table)+1-top$table$rank) * sign(top$table$effect)
        else
            stats <- top$table$effect
        stats[is.na(stats)] <- 0.0
        names(stats) <- top$table$name
        if (!signed) stats <- abs(stats)
        
        confident_sign <- ifelse(is.na(top$table$confect),0,sign(top$table$effect))
        names(confident_sign) <- top$table$name
        
        rich <- enrichment(stats, organism$orgdb, ad=ad, confident_sign=confident_sign)
        
        suffix <- if (signed) "_signed" else "_unsigned"
        
        saveRDS(rich, file.path(out_dir, paste0("enrichment",suffix,".rds")))
        
        df <- mutate(rich, gene_ids=map_chr(gene_ids, paste, collapse="/"))
        write_csv(df, file.path(out_dir, paste0("enrichment",suffix,".csv")))
        
        filter(df, novel_up+novel_down > size/5) %>%
            write_csv(file.path(out_dir, paste0("enrichment",suffix,"_brief.csv")))
    }
}

#' @export
do_components_tests <- function(out_dir, input, organism, confect=FALSE, rank=FALSE, ad=FALSE, fdr=0.05) {
    input <- load_banquet(input)
    organism <- load_banquet(organism)

    comp_seq <- input$comp_seq
    for(p in seq_along(comp_seq)) {
        message(p)
        comp_names <- comp_seq[[p]]$ind_components
        for(comp in seq_len(p)) {
            comp_name <- paste0("C",comp)
            message(comp_name)
            do_test(file.path(out_dir, paste0("p",p), comp_name),
                input$weitrix, organism, comp_seq[[p]]$col, comp_name,
                confect=confect, rank=rank, ad=ad, fdr=fdr)
        }
    }
}







