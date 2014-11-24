#!/usr/bin/env ruby

#encoding: utf-8

class FitnessScore
	require 'rinruby'
	require_relative 'lib/snp_dist'

	### Count/ratio method ################################################################################
		# Input 0: Array of SNP positions
		# Input 1: Number of breaks (divisions) in the genome to count the number of SNPs in
		# Input 2: The length of the genome
		# Output: Array of number of SNPs in each genome division
		def self.count(snp_pos, div, genome_length)
			myr = RinRuby.new(echo = false)
			myr.assign 'snp_pos', snp_pos
			myr.assign 'div', div.to_i
			myr.assign 'l', genome_length
			myr.eval 'breaks <- c(0)
			for(i in 1:div){
			  breaks <- c(breaks,(l/div)*i)
			}'
			myr.eval 'counts <- hist(snp_pos, breaks=breaks, plot=FALSE)$counts'
			counts = myr.pull 'counts'
			myr.quit
			return counts
		end

		# Input 0: Array of homozygous SNP positions
		# Input 1: Array of heterozygous SNP positions
		# Input 2: Number of breaks (divisions) in the genome to count the number of SNPs in
		# Input 3: The length of the genome
		# Output: Array of ratios (floats) of homozygous to heterozygous SNPs for each division of the genome permutation
		def self.ratio(hm, ht, div, genome_length)
			hm_count = FitnessScore::count(hm, div, genome_length)
			ht_count = FitnessScore::count(ht, div, genome_length)
			x = 0
			ratios = []
			div.times do
				count_ratio = ((hm_count[x] + 1).to_f / (ht_count[x] + 1).to_f) # a measure of ratio
				ratios << count_ratio
				x+=1
			end
			return ratios
		end

		# Input 0: Array of homozygous SNP positions
		# Input 1: Array of heterozygous SNP positions
		# Input 2: Number of breaks (divisions) in the genome to count the number of SNPs in
		# Input 3: The length of the genome
		# Input 4: Array of expected ratios (floats) of homozygous to heterozygous SNPs for each division of the genome
		# Output: Float between 0.0 and 1.0 where closely matching inputs are closer to 1.0 (pearson correlation)
		def self.count_ratio(hm, ht, div, genome_length, expected_ratios)
			ratios = ratio(hm, ht, div, genome_length)
			myr = RinRuby.new(echo = false)
			myr.assign 'x', expected_ratios
			myr.assign 'y', ratios
			myr.eval 'score <- abs(cor(x,y))'
			fitness_score = myr.pull 'score'
			myr.quit
			return fitness_score
		end
end
