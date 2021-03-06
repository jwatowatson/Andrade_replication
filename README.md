# Andrade_replication

The code in this repository provides a re-analysis of RNA seq data from [Andrade et al](https://www.nature.com/articles/s41591-020-1084-0) and some simulations to test sensitivity of the underlying algorithm that estimates the hours post invasion (hpi) from in vivo transcriptomics data. The RMarkdown script *hpi_estimation.Rmd* downloads the data; does some cleaning of gene names; fits the Lemieux et al MLE model; does some gene expression comparisons based on ranks; and then does some simulations that show pathological behaviour of the Lemieux et al algorithm in asynchronous infections.


Key messages:

* Comparisons of developmental stages in asymptomatic and symptomatic patients are confounded by selection bias (fever in uncomplicated malaria influences treatment seeking behaviour; fever synchronises the infection; fever can trigger cytoadherence; patients are selected on the basis of parasitaemia). A key difference is in the synchronicity of the parasite stages across infected red cells: asymptomatics are probably very asynchronous, whereas symptomatics are probably highly synchronous (especially if you select patients on the basis of high parasitaemia, see discussion below).
* Andrade et al have not explored the influence of gametocyte gene expression: gametocyte proportions in asymptomatics will be much higher than in symptomatics
* The maximum likelihood method from Lemieux et al used to estimate hpi basically breaks down when you simulate asynchronous infections (give wildly biased results)


## Synchronicity of infection

Brigadier Sir Neil Hamilton Fairley did some pretty cool experiments at the Cairns experimental station (*Sidelights on malaria in man obtained by subinoculation experiments*, 1947). Vivax doesn't sequester and he showed a log-linear increase in parasite densities (no fluctuations in density). Falciparum does sequester so we see a sine-wave increase in parasite densities:


![Sine-wave increase in parasite densities from Hamilton Fairley's experiments](Fairley1947.png)

The sine-wave pattern is explained by two things:

* Sequestration of mature stages (large rings onwards)
* Relative synchronicity of infection (this is helped by fever): if the infection was asynchronous we would not see the sine-wave (cancelled out by the variance in ages)

A more recent (but still 30 years old) paper by Nick White and colleagues presents a simple mathematical model for understanding the effects of synchronicity in falciparum infections (*The effects of multiplication and synchronicity on the vascular distribution of parasites in falciparum malaria*. 1992). They put together a mathematical model of peripheral parasite density as a function of the mean age (hours post invasion) and the synchronicity (standard deviation of the age disribution). Parasites ages are assumed to be normally distributed (modulo 48) with varying mean age and varying standard deviation. They have some really interesting data from 4 patients who were treated with cipro (completely ineffective). These patients were then given rescue quinine treatment when the parasitaemia didn't go down. You can see the amazing sine-wave parasitaemia curves in the 4 patients in Figure 8 of the paper:

![](White_Chapman_Watt_Fig8.png)



## Simulations

We simulate transcriptomics data using the reference dataset from Bozdech et al. The simulation algorithm does the following:

* Input: a mean age $\mu\in[0,48]$; a standard deviation in age $\sigma>0$; a vector of circulation probabilities for each hour 0 until 48.
* We calculate the proportion of parasites ages for each hourly block from 0 to 48 under the assumption that the parasite ages are normally distributed with mean $\mu$ and variance $\sigma^2$, modulo 48.
* We then sample parasites from the expected hourly blocks and use the reference transcriptome to construct an aggregate trancriptional profile
* We do the same but after multiplying the expected proportions per hour by the circulating probabilities
* We compute the mle using the method from Lemieux et al on both simulated gene expression datasets






