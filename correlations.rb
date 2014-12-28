#!/usr/bin/env ruby

#encoding: utf-

require_relative 'lib/fitness_score'
require_relative 'lib/reform_ratio'
require 'pmeth'
require 'csv'
dataset = ARGV[0]
name = ARGV[1]


## Files ################################################
vcf_file = "genomes/#{dataset}/snps.vcf"
fasta_file = "genomes/#{dataset}/frags.fasta"
fasta_shuffle = "genomes/#{dataset}/frags_shuffled.fasta"
#############################################################

fasta = ReformRatio::fasta_array(fasta_file)

snp_data = ReformRatio.get_snp_data(vcf_file)

ht, hm = ReformRatio.perm_pos(fasta, snp_data)

fasta_shuffle = ReformRatio::fasta_array(fasta_shuffle)

ht2, hm2 =ReformRatio.perm_pos(fasta_shuffle, snp_data)

genome_length = ReformRatio::genome_length(fasta_file)

# CSV.open("genomes/#{dataset}/#{name}.csv", "ab") do |csv|
# 	csv << ["Genome length", "Divisions", "Fitness"]
# end 
list = (526..1000)
fitness = []
list.each do |div|
	if div % 2 == 0 
		expected_ratios = FitnessScore.ratio(hm, ht, div, genome_length) # Array of expected ratios (floats) of homozygous to heterozygous SNPs for each division of the genome
		ratios = FitnessScore.ratio(ht2, hm2, div, genome_length)
		myr = RinRuby.new(echo = false)
		myr.assign 'x', expected_ratios
		myr.assign 'y', ratios
		myr.eval 'score <- abs(cor(x,y))'
		fitness_score = myr.pull 'score'
		myr.quit
		puts "Calculating fitness for #{div} divisions"
		CSV.open("genomes/#{dataset}/#{name}.csv", "ab") do |csv|
			csv << [genome_length, div, fitness_score]
		end
	end
end
