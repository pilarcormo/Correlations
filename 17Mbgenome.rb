#encoding utf-8
require_relative 'lib/model_genome'
require_relative 'lib/write_it'

require_relative 'lib/reform_ratio'
name = ARGV[0]
require 'PMeth'

# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "/Correlations/genomes/#{name}"))

# Create the lists of homozygous and heterozygous SNPs
fasta_file = "TAIR10_chr4.fasta"
genome_length = ReformRatio::genome_length(fasta_file)

hm_r = "hm <- rnorm(18000, 10000000, 1000)" # Causative SNP at/near 1000
ht_r = "ht <- runif(18000, 1, #{genome_length})" # Genome length is same as arabidopisis chromosome 4
hm, ht = ModelGenome::get_snps(hm_r, ht_r)
snp_pos = [hm, ht].flatten

puts "There are #{hm.length} homozygous SNPs"
puts "There are #{ht.length} heterozygous SNPs"
puts "Is there a SNP at the centre of the distribution? -- #{snp_pos.include?(10000000)}"

arabidopsis_c4 = ModelGenome::fasta_to_char_array(fasta_file)
puts "Fragmenting arabidopsis_c4..."
contig_size = 18000# 10-20kb
frags = ModelGenome::get_frags(arabidopsis_c4, contig_size)
puts "Done!"
puts "Arabidopsis chr4 length: #{arabidopsis_c4.length}/ kb"
puts "Fragmented seq   length: #{frags.join.length} = close enough? You decide."
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"

# Get the positions of the SNPs on fragments
pos_on_frags, snp_pos_all = ModelGenome::pos_each_frag(snp_pos, frags)

fastaformat_array = ModelGenome::fasta_array(frags)
fastaformat_array_shuf = fastaformat_array.shuffle # shuffle it to show that the order doesn't need to be conserved when working out density later on

vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

#Shuffle the original order of fragments by using the adjacent_swap method. 
start_pop = []
start_pop << fastaformat_array
10.times do
	adj_pop = []
	start_pop.each do |perm|
		70.times do
			perm = PMeth.adjacent_swap(perm) 
		end
		adj_pop << perm # need this population to be the next starting population
	end
	start_pop = adj_pop
end

WriteIt::write_data("genomes/#{name}/frags.fasta", fastaformat_array)
WriteIt::write_data("genomes/#{name}/snps.vcf", vcf)
WriteIt::write_data("genomes/#{name}/frags_shuffled.fasta", fastaformat_array_shuf)
WriteIt::write_txt("genomes/#{name}/info", [hm_r, ht_r, "Contig size = #{contig_size}"])
WriteIt::write_txt("genomes/#{name}/hm_snps", hm)
WriteIt::write_txt("genomes/#{name}/ht_snps", ht)