#!/usr/bin/env ruby
#encoding: utf-8

require 'PMeth'
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'

require_relative 'count_ratio'
require_relative 'lib/reform_ratio'
require 'PMeth'
require_relative 'swaps'
require 'pp'
dataset = ARGV[0]
size = ARGV[1].to_i # size of each population of permuations that are progressively further from correct
pop_num = ARGV[2].to_i # number of populations
swap_num = ARGV[3].to_i # number on adjacent swaps performed on permutations between each population



## Files ################################################
vcf_file = "genomes/#{dataset}/snps.vcf"
fasta_file = "genomes/#{dataset}/frags.fasta"
fasta_shuffle = "genomes/#{dataset}/frags_shuffle.fasta"
#############################################################

fasta = ReformRatio::fasta_array(fasta_file)



# class Swap

# def self.adjacent_swaps (dataset, size, pop_num, swap_num)
fasta = ReformRatio::fasta_array("genomes/#{dataset}/frags.fasta") # correct permutation
snp_data = ReformRatio::get_snp_data("genomes/#{dataset}/snps.vcf")
genome_length = ReformRatio::genome_length("genomes/#{dataset}/frags.fasta")


start_pop = []
# size.times do
start_pop << fasta
# end

puts start_pop

pop_num.times do
	adj_pop = []
	start_pop.each do |perm|
		swap_num.times do
			perm = PMeth.adjacent_swap(perm) 
		end
		# puts "adjacent_swaps: another #{swap_num} for pop: #{x}"
		adj_pop << perm # need this population to be the next starting population

	end
	start_pop = adj_pop

end
puts start_pop

# end
# end

