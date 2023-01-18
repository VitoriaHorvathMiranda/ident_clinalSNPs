## SNP Calling

Os SNPs foram identificados com o
[snape-pooled](https://github.com/EmanueleRaineri/snape-pooled). Esse
programa considera a qualidade e a profundidade do sequenciamento e
devolve a frequência das bases de todas as posições e a “probabilidade”
dessa frequência estar correta. Esse arquivo foi filtrado para incluir
apenas as posições cuja frequência indicada tinha “probabilidade” maior
ou igual a 0.9 e para incluir apenas os cromossomos 2L, 2R, 3L, 3R e X:

    awk '(($9>=0.9) && ($1==\"2L\" || $1==\"2R\" || $1==\"3L\" || $1==\"3R\" || $1==\"X\"))' {input} > {output}

O que gerou essa tabela para cada população:

    ##     chrom position ref_base ref_count alt_count ref_qual alt_qual bases   prob
    ##  1:    2L     5037        A         0         1       70       70     G 0.9998
    ##  2:    2L     5357        C         4         1       70       70    CG 0.9277
    ##  3:    2L     5372        T         2         3       70       66    AT 1.0000
    ##  4:    2L     5390        T         3         2       70       70    TA 1.0000
    ##  5:    2L     5403        C         4         1       67       70    CG 0.9277
    ##  6:    2L     5762        T         3         4       66       70    CT 1.0000
    ##  7:    2L     5867        A         4         1       70       70    AG 0.9277
    ##  8:    2L     5889        T         3         1       66       70    TC 0.9452
    ##  9:    2L     5892        A         3         1       70       70    AG 0.9452
    ## 10:    2L     5904        C         3         1       70       70    CA 0.9452
    ##          p(1)   freq
    ##  1: 9.480e-01 0.9899
    ##  2: 5.748e-13 0.1913
    ##  3: 4.522e-05 0.5997
    ##  4: 9.032e-09 0.4001
    ##  5: 9.115e-12 0.1914
    ##  6: 7.146e-07 0.5717
    ##  7: 5.748e-13 0.1913
    ##  8: 3.463e-08 0.2419
    ##  9: 2.184e-09 0.2418
    ## 10: 2.184e-09 0.2418

## Effective Number of Chromosomes

Antes de criar os modelos para os snps clinais, calculei o número
efetivo de cromossomos (NE) para cada posição da tabela anterior e para
cada população. Para isso usei o
[script\_test.R](https://github.com/VitoriaHorvathMiranda/ident_clinalSNPs/blob/main/script_test.R)
no Snakemake. Para fazer isso, obtive a profundidade de leitura de cada
posição com o samtools:

    samtools depth -aa -q 20 -Q 30 {input} -o {output}

Calculei o número efetivo de cromossomos da seguinte forma:

    merged <- merge(snapePool, samtoolsDepth, by = c("chrom", "position"), all.x = TRUE)

    #gets chrom number
    chrom_sampled <- xargs$ChromNumber

    #calcula o NE 
    merged[, NE := ((1/depth) + (1/chrom_sampled))^-1]
    merged[, latitude := xargs$Latitude]
    merged[, population := xargs$PopName]

O número de cromossomos amostrados em cada população estava na tabela de
metadados. No final desse script tive a seguinte tabela para cada
população:

    ##     chrom position   freq depth        NE latitude population
    ##  1:    2L     5037 0.9899     1 0.9937888    33.39      ESC97
    ##  2:    2L     5357 0.1913     4 3.9024390    33.39      ESC97
    ##  3:    2L     5372 0.5997     5 4.8484848    33.39      ESC97
    ##  4:    2L     5390 0.4001     5 4.8484848    33.39      ESC97
    ##  5:    2L     5403 0.1914     6 5.7831325    33.39      ESC97
    ##  6:    2L     5762 0.5717     8 7.6190476    33.39      ESC97
    ##  7:    2L     5867 0.1913     4 3.9024390    33.39      ESC97
    ##  8:    2L     5889 0.2419     3 2.9447853    33.39      ESC97
    ##  9:    2L     5892 0.2418     3 2.9447853    33.39      ESC97
    ## 10:    2L     5904 0.2418     3 2.9447853    33.39      ESC97

Em seguida usei o script
[join\_time\_pops.R](https://github.com/VitoriaHorvathMiranda/ident_clinalSNPs/blob/main/join_time_pops_script.R)
para juntar todas as populações de um mesmo período em uma única tabela,
então no final formei três tabelas, uma com as populações de 1997, uma
com as populações de 2009/2010 e uma última com as populações de 2017.

## Glm

Para identificar os SNPs clinais era preciso correr um GLM para cada SNP
de cada período. Mas antes disso, filtrei os SNPs que tinham
profundidade de leitura igual a zero. Pode parecer estranho que o
snape-pooled tenha identificado snps com profundidade de leitura igual a
zero, mas isso ocorreu porque quando calculei a profundidade com o
`samtools depth` usei um filtro de qualidade de leitura e qualidade de
mapeamento. Esse filtro considerou a profundidade como zero para todas
as posições que não alcançavam o valor mínimo de qualidade. Mesmo assim
o snape-pooled conseguiu chamar esses SNPs, mas achei melhor filtrá-los.
Também exclui os SNPs que foram chamados apenas em duas ou menos
populações, e precisei aninhar cada SNP na tabela.

    nested_snps <- joined_pops[depth>0,][,group_nest_dt(.SD, position2)][, 
                                           n_pops := purrr::map_dbl(data, nrow)][
                                           n_pops>2,][,
                                           !c("n_pops"), with = FALSE]

    ##     position2              data
    ##  1:   2L:5317 <data.table[3x5]>
    ##  2:   2L:5357 <data.table[3x5]>
    ##  3:   2L:5372 <data.table[7x5]>
    ##  4:   2L:5390 <data.table[7x5]>
    ##  5:   2L:5403 <data.table[7x5]>
    ##  6:   2L:5465 <data.table[6x5]>
    ##  7:   2L:5598 <data.table[4x5]>
    ##  8:   2L:5762 <data.table[7x5]>
    ##  9:   2L:5867 <data.table[7x5]>
    ## 10:   2L:5889 <data.table[7x5]>

Corri o seguinte glm para cada SNP:

    #runs glm for each snp
    nested_snps[, models := purrr::map(data, ~ glm(freq~latitude, 
                                                   weights = NE,
                                                   data = .x,
                                                   family = binomial()))]
    head(nested_snps)

    ##    position2              data    models
    ## 1:   2L:5317 <data.table[3x5]> <glm[30]>
    ## 2:   2L:5357 <data.table[3x5]> <glm[30]>
    ## 3:   2L:5372 <data.table[7x5]> <glm[30]>
    ## 4:   2L:5390 <data.table[7x5]> <glm[30]>
    ## 5:   2L:5403 <data.table[7x5]> <glm[30]>
    ## 6:   2L:5465 <data.table[6x5]> <glm[30]>

olhando para um modelo:

    pluck(nested_snps[[3]]) %>% pluck(1)

    ## 
    ## Call:  glm(formula = freq ~ latitude, family = binomial(), data = .x, 
    ##     weights = NE)
    ## 
    ## Coefficients:
    ## (Intercept)     latitude  
    ##    -4.35989      0.04843  
    ## 
    ## Degrees of Freedom: 2 Total (i.e. Null);  1 Residual
    ## Null Deviance:       0.731 
    ## Residual Deviance: 0.02412   AIC: 13.21

    nested_snps[, summary := purrr::map(models, 
                                            ~summary(.x))]
    pluck(nested_snps[[4]]) %>% pluck(1)

    ## 
    ## Call:
    ## glm(formula = freq ~ latitude, family = binomial(), data = .x, 
    ##     weights = NE)
    ## 
    ## Deviance Residuals: 
    ##        1         2         3  
    ##  0.02854  -0.11004   0.10583  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept) -4.35989    2.49140   -1.75   0.0801 .
    ## latitude     0.04843    0.06133    0.79   0.4297  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 0.731028  on 2  degrees of freedom
    ## Residual deviance: 0.024123  on 1  degrees of freedom
    ## AIC: 13.212
    ## 
    ## Number of Fisher Scoring iterations: 4

Então extrai o coeficiente de inclinação e o p-valor associado ao
coeficiente de inclinação de cada glm:

    #gets lat coefficient
    nested_snps[, coefficients := purrr::map_dbl(models, ~coef(.x) %>% pluck("latitude"))]

    #gets p-value
    nested_snps[, p_value := purrr::map_dbl(models, 
                                             ~summary(.x) %>% 
                                               pluck("coefficients") %>% pluck(8))]
    head(nested_snps, n=3)

    ##    position2              data    models           summary coefficients
    ## 1:   2L:5317 <data.table[3x5]> <glm[30]> <summary.glm[17]>   0.04843086
    ## 2:   2L:5357 <data.table[3x5]> <glm[30]> <summary.glm[17]>  -0.04014949
    ## 3:   2L:5372 <data.table[7x5]> <glm[30]> <summary.glm[17]>  -0.01348586
    ##      p_value
    ## 1: 0.4297093
    ## 2: 0.6503347
    ## 3: 0.6314901

Nesta última parte eu estou bastante insegura em relação ao p-valor que
extraí, é esse mesmo? Toda essa última etapa foi feita no script
[glm\_script.test](https://github.com/VitoriaHorvathMiranda/ident_clinalSNPs/blob/main/glm_script.R)

## P-value cutoff

No mesmo script de antes
([glm\_script.R](https://github.com/VitoriaHorvathMiranda/ident_clinalSNPs/blob/main/glm_script.R)),
calculei o p-valor para um FDR de 10% (foi um valor arbitrário para
testar, acho que temos que pensar melhor nisso). Isso foi feito com o
pacote [qvalue](https://github.com/StoreyLab/qvalue), criado pelos
autores do paper [Storey and Tibshirani,
2003](https://www.pnas.org/doi/abs/10.1073/pnas.1530509100).

    p_values <- nested_snps$p_value
    qobj <- qvalue(p = p_values)
    p_value_cutoff <- max(qobj$pvalues[qobj$qvalues <= 0.1])

Em seguida filtrei todos os SNPs cujo p-valor do glm era menor que o
p-valor do cut-off:

    #gets snps below that p-value
    significant_snps <- nested_snps[p_value <= p_value_cutoff,]

Obtive a seguinte tabela:

    ##     position2 coefficients      p_value
    ##  1:   2L:5465  -0.08990523 0.0070812282
    ##  2:  2L:42709   0.11863900 0.0053168082
    ##  3:  2L:55977  -0.40931142 0.0046357487
    ##  4:  2L:56815   0.17012102 0.0062187496
    ##  5:  2L:84886   0.10610724 0.0043126608
    ##  6: 2L:136625   0.46104193 0.0047295115
    ##  7: 2L:143327  -0.21965328 0.0072895720
    ##  8: 2L:144125  -0.11816785 0.0070629710
    ##  9: 2L:144142  -0.11122934 0.0003831638
    ## 10: 2L:144188  -0.09107743 0.0043110069

Eu só fiz esses dois últimos passos para os SNPs de 1997, porque ainda
estava testando o script e ainda não tentei rodar ele com o Snakemake.
