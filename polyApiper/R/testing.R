
#' @export
get_gene_sets <- function(orgdb=NULL, go=FALSE, msigdb=FALSE, custom=NULL) {
    # TODO: fixme

    gene_term <- tibble(name=character(0), term_name=character(0))
    term_info <- tibble(term_name=character(0), desc=character(0), collection=character(0))
        
    if (!is.null(custom)) {
        gene_term <- tibble(
                name=map(custom, geneIds),
                term_name=names(custom)) %>%
            pmap_df(tibble) %>%
            bind_rows(gene_term)
        
        term_info <- tibble(
                term_name=names(custom),
                desc=map_chr(custom, description),
                collection="custom") %>%
            bind_rows(term_info)
    }
    
    if (msigdb) {
        entrez <- AnnotationDbi::select(
                orgdb, keys=AnnotationDbi::keys(orgdb, "ENSEMBL"), keytype="ENSEMBL", columns="ENTREZID") %>%
            filter(!is.na(ENTREZID)) %>%
            transmute(name=ENSEMBL, entrez_gene=as.integer(ENTREZID))
        organism <- filter(metadata(orgdb), name=="ORGANISM")$value
        msig <- msigdbr(species=organism) %>%
            inner_join(entrez, "entrez_gene")
        
        # Omit GO terms if they are going to be included separately
        # (MSigDB only has a subset of GO terms)
        if (go)
           msig <- filter(msig, gs_cat != "C5")
        
        gene_term <- msig %>%
            select(name, term_name=gs_id) %>%
            bind_rows(gene_term)
        
        term_info <- distinct(msig, gs_id, gs_name, gs_cat, gs_subcat) %>%
            transmute(term_name=gs_id, desc=gs_name, collection=paste0(gs_cat,":",gs_subcat)) %>%
            bind_rows(term_info)

    } 
    
    if (go) {
        gene_term_go <- AnnotationDbi::select(
                orgdb, keys=AnnotationDbi::keys(orgdb, "ENSEMBL"), keytype="ENSEMBL", columns="GOALL") %>%
            select(name=ENSEMBL, term_name=GOALL) %>%
            filter(!is.na(term_name)) %>%
            distinct(name, term_name) # Duplicate entries, by evidence type?
        
        gene_term <- bind_rows(gene_term_go, gene_term)
        
        term_info <- AnnotationDbi::select(
            GO.db, keys=unique(gene_term_go$term_name), keytype="GOID",
            columns=c("GOID","TERM","ONTOLOGY")) %>%
            select(term_name=GOID, desc=TERM, collection=ONTOLOGY) %>%
            bind_rows(term_info)
    }
    
    list(gene_term=gene_term, term_info=term_info)
}


#' Confident-correlation gene-set enrichment
#'
#' Effect size is the Pearson correlation between the indicator variable for the gene-set and the statistic.
#'
enrichment <- function(stats, confident_sign, gene_sets, min_size=2, ad=FALSE) {
       
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

    result <- gene_sets$gene_term %>%
        filter(name %in% names(stats)) %>%
        mutate(stat=stats[name]) %>%
        group_by(term_name) %>%
        # A gene set containing all genes can be discarded:
        filter(length(name) >= min_size, length(name)+min_size <= length(stats)) %>%
        summarize(data=list(doit(name)), size=length(name), names=list(name)) %>%
        unnest(data) %>%
        arrange(-abs(effect)) %>%
        left_join(gene_sets$term_info, by="term_name")

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
            term_name, desc, collection, gene_ids=names)

    # On reflection, this doesn't seem to make much sense:
    #
    #top <- normal_confects(result$effect, result$se, result$df)
    #top$table <- bind_cols(top$table,
    #    select(result, diff, size, goid=GOID, term=TERM, ontology=ONTOLOGY, gene_ids=names)[top$table$index,,drop=F])
    #top$effect_desc <- "correlation"
    #top
}


enrichment_fgsea <- function(stats, confident_sign, gene_sets, min_size=10, max_size=500) {
    pathways <- gene_sets$gene_term %>%
        filter(name %in% names(stats)) %>%
        group_by(term_name) %>%
        filter(n() >= min_size, n() <= max_size) %>%
        ungroup() %>%
        { split(.$name, .$term_name) }

    message(length(pathways), " pathways")
    iters <- max(1e5, length(pathways)*100)
    message(iters, " samples will be taken")

    result <- fgsea(pathways, stats, iters)
    
    result %>%
        select(pval,padj,NES,ES,size, term_name=pathway) %>%
        arrange(pval) %>%
        left_join(gene_sets$term_info, "term_name") %>%
        mutate(gene_ids=pathways[term_name])
}


enrichment_wilcox <- function(stats, confident_sign, gene_sets, min_size=2, max_prop=1/3) {
    max_size <- length(stats)*max_prop

    doit <- function(names) {
        mask <- names(stats) %in% names
        bg_mask <- !mask
        tibble(
            size=length(names),
            p=wilcoxGST(mask, stats, type="f"),
            pos=mean(confident_sign[mask]>0),
            neg=mean(confident_sign[mask]<0),
            bg_pos=mean(confident_sign[bg_mask]>0),
            bg_neg=mean(confident_sign[bg_mask]<0),
            gene_ids=list(names))
    }


    gene_sets$gene_term %>%
        filter(name %in% names(stats)) %>%
        group_by(term_name) %>%
        filter(n() >= min_size, n() <= max_size) %>%
        summarize(data=list(doit(name))) %>%
        unnest(cols="data") %>%
        arrange(p) %>%
        mutate(padj = p.adjust(p, "fdr")) %>%
        left_join(gene_sets$term_info, "term_name") %>%
        select(p,padj,size,pos,neg,bg_pos,bg_neg,term_name,desc,collection,gene_ids)
}


#' @export
do_test <- function(out_dir, fit, coef, gene_sets=NULL, fdr=0.05, trend=FALSE, ...) {    
    top <- limma_confects(fit, coef, full=TRUE, fdr=fdr, trend=trend)

    ensure_dir(out_dir)
    saveRDS(top, file.path(out_dir, "confects.rds"))
    write_csv(top$table[,-c(1,2)], file.path(out_dir, "confects.csv"))
    
    if (is.null(gene_sets))
        return()

    stats <- nrow(top$table)-top$table$rank
    names(stats) <- top$table$name
        
    confident_sign <- ifelse(is.na(top$table$confect),0,sign(top$table$effect))
    names(confident_sign) <- top$table$name
        
    rich <- enrichment_wilcox(stats, confident_sign=confident_sign, gene_sets=gene_sets, ...)
        
    saveRDS(rich, file.path(out_dir, "enrichment.rds"))
        
    df <- mutate(rich, gene_ids=map_chr(gene_ids, paste, collapse="/"))
    write_csv(df, file.path(out_dir, "enrichment.csv"))

    #confect=TRUE, rank=FALSE,   
    #for(signed in c(F)) {
    #    if (confect)
    #        stats <- top$table$confect
    #    else if (rank)
    #        stats <- (nrow(top$table)+1-top$table$rank) * sign(top$table$effect)
    #    else
    #        stats <- top$table$effect
    #    stats[is.na(stats)] <- 0.0
    #    names(stats) <- top$table$name
    #    if (!signed) stats <- abs(stats)
    #    
    #    confident_sign <- ifelse(is.na(top$table$confect),0,sign(top$table$effect))
    #    names(confident_sign) <- top$table$name
    #    
    #    rich <- enrichment_wilcox(stats, confident_sign=confident_sign, gene_sets=gene_sets, ...)
    #    
    #    suffix <- if (signed) "_signed" else "_unsigned"
    #    
    #    saveRDS(rich, file.path(out_dir, paste0("enrichment",suffix,".rds")))
    #    
    #    df <- mutate(rich, gene_ids=map_chr(gene_ids, paste, collapse="/"))
    #    write_csv(df, file.path(out_dir, paste0("enrichment",suffix,".csv")))
    #    
    #    #filter(df, novel_up+novel_down > size/5) %>%
    #    #    write_csv(file.path(out_dir, paste0("enrichment",suffix,"_brief.csv")))
    #}
}


#' @export
do_tests <- function(out_dir, weitrix, comp, coefs=NULL, ...) {
    cal_weitrix <- weitrix_calibrate_trend(weitrix, comp)
    fit <- lmFit(weitrix_elist(cal_weitrix), comp$col)

    if (is.null(coefs))
        coefs <- colnames(comp$col)[comp$ind_components]

    for(coef in coefs) {
        message(coef)
        do_test(file.path(out_dir, coef),
            fit, coef, ...)
    }
}

#' @export
do_components_tests <- function(out_dir, input, ...) {
    input <- load_banquet(input)
    
    comp_seq <- input$comp_seq
    for(p in seq_along(comp_seq)) {
        message(paste0("p",p))
        do_tests(file.path(out_dir, paste0("p",p)),
            input$weitrix, comp_seq[[p]], ...)
    }
}







