Correlations
============

The [distance metrics](https://github.com/pilarcormo/small_genomes_SNPs/tree/master/Results/Rplots.%20Distances) used to evaluate the performance of the genetic algorithm turn out to be not discriminating enough for the purpose of the experiment. 

However, we want to test if the fitness score calculated in the count ratio method can be used as a measure of the efficiency of the algorithm to order contigs. 

The count ratio method, for each genome division, defines the hm and ht SNPs and calculate a ratio hm/ht. Then, it compares all the ratios obtained with the expected ratios by using the absolute Pearson's correlation coefficient that is _"a measure of the linear correlation (dependence) between two variables X_ (expected ratios) _and Y_ (permutation ratios), _giving a value between +1 and −1 inclusive, where 1 is total positive correlation, 0 is no correlation, and −1 is total negative correlation"._ In this case, as it is the aboslute correlation, we only have values oscillating between 1 and 0 (correlation or not).

Therefore, we are going to measure the influence of the genome divisions and the genome length on this fitness score (r). Instead of using the algorithm, I compare the correctly ordered genome with a permutation obtained from it by using the [PMeth.adjacent_swap method](https://github.com/pilarcormo/pmeth/blob/master/lib/pmeth.rb). I used 70 adjacent swaps and 10 populations to generate the example permutation. 
