Results
===


First approach
---

![Image](https://github.com/pilarcormo/Correlations/blob/master/Results/Rplot_First_brewer.png?raw=false)



To create the permutated genomes I did:

```
10.times do
	start_pop.each do |perm|
		70.times do
			perm = PMeth.adjacent_swap(perm) 
```

So, 700 adjacent swaps change the  correctly ordered genome. In every case, the number of homozygous SNP is kept constant while I changed the number of heterozygous SNP as shown in the titles. 

In every case, the increase in the number of divisions decrease the correlation between the ratios (decrease the fitness)


I am tented to say that the correlation is in general better when we have a higher number of heterozygous SNPs (1/500). I kept the number of homozygous SNP more or less constant for all the genome sizes in these examples. 

Effect of the contig size on the correlation 
----


![Image](https://github.com/pilarcormo/Correlations/blob/master/Results/Rplot_Contig_brewer.png?raw=false)


To improve the quality of the plots and compare them more effectively, the permutation of the genome used in all these cases was created by:

```
2.times do
	start_pop.each do |perm|
		70.times do
			perm = PMeth.adjacent_swap(perm) 
```

So, only 140 adjacent swaps change the correctly ordered genome (instead of 700). Doing this, I avoid having specific genomes with bad correlation probably due to the randomness of the permutation (as seen above).

Also, to create these genomes I changed the number of homozygous and heterozygous SNPs so I have the same number of both types of SNP. 

We observe a clear trend of improvement in the correlation (red colour) when increasing the size of the genome for high number of divisions. If we focus on 250 divisions, for instance, this improvement occurs regardless the number of contigs. The value of r is much better for lower number of divisions for all the genomes.

There is no significant effect of the number of contigs on the fitness score. Again, we observe an increase in the fitness score for all the contig sizes for larger genomes. The trend is similar for the 3 contig sizes. Therefore, the contig size does no affect the fitness. 



Effect of the SNP density 
----
For the following results, I used 700 contigs to build the genomes. 

###1. Same number of homozygous and heterozygous SNPs. 
![Image](https://github.com/pilarcormo/Correlations/blob/master/Results/Rplot.SNP_density_brewer.png?raw=false)


###2. Increase ht SNP density while keeping the number of hm SNPs constant

![Image](https://github.com/pilarcormo/Correlations/blob/master/Results/Rplot.Ht_Hmconstant.png?raw=false)

It seems that a high SNP frequency has a possitive effect in the fitness score. 


Conclusions 
----

- In general, the correlation between the ratios worsens when using more than 250 divisions.
- The contig size has no significant effect on the correlation.
- When the number of divisions increases, the best correlations are observed in larger genomes. 
- In general, the SNP frequency has no effect on the correlation, although the highest frequency used gave better correlation results. 


